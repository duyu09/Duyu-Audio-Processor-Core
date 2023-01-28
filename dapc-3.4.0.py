# Python 3.x source code of Duyu-PCM-Processor-Core software system.
# Copyright (c) 2020~2022 DuYu (No.202103180009), Faculty of Computer Science & Technology, Qilu University of Technology.
# Copyright (c) 2010 Bill Cox. (Sonic Library (sonic.dll, Built with TDM-GCC 4.9.2 64-bit Release by Duyu.))
# Copyright (c) 2000~2021 the FFmpeg developers(FFmpeg (ffmpeg.exe, Built with GCC 11.2.0 (Rev5, Built by MSYS2 project))
# Version 3.4.0 __August 28, 2022.

# This source code file should be process with PyInstaller(version>=5.1, architecture=x86x64, OS=WinAll).
# FFmpeg.exe and sonic.dll should be encapsuled in EXE file.
# FFmpeg.exe, sonic.dll & this file should be in the same directory if you run the source code alone.

# Import Librarys.
import ctypes
import os.path
import subprocess
import sys
from os import popen
from sys import argv
from sys import stderr
from sys import exit
from numpy import append as npappend
from numpy import zeros as npzeros
from numpy import float64 as npfloat64
from numpy import float32 as npfloat32
from numpy import int16 as npint16
from numpy import int32 as npint32
from numpy import array as nparray
from numpy import linspace as nplinspace
from numpy import insert as npinsert
from numpy import vstack as npvstack
from numpy import divide as npdivide
from numpy import add as npadd
from numpy import log10 as nplog10
from numpy import clip as npclip
from numpy import abs as npabs
from numpy import interp as npinterp
from numpy import arange as nparange
from numpy.fft import rfft as nprfft
from numpy.ctypeslib import ndpointer
from scipy.signal import medfilt
from scipy.io.wavfile import read as wavRead
from scipy.io.wavfile import write as wavWrite
from argparse import ArgumentParser
from warnings import filterwarnings

# Global const variables declare
NAME_AND_VERSION = "Duyu PCM-WAV Audio Processor Core v3.4.0"
BUILDER = "Built and packaged with PyInstaller 5.1"
BUILD_TIME = "__August 28 2022."
COPYRIGHT_STATEMENT = "Copyright (c) 2020~2022 DuYu (No.202103180009),\
Faculty of Computer Science & Technology, Qilu University of Technology."


# APIs

def getTypeFactor(gtf_dataType):
    typeFactor = 2 ** 15
    if gtf_dataType == npint16:
        typeFactor = 2 ** 15
    elif gtf_dataType == npint32:
        typeFactor = 2 ** 31
    elif gtf_dataType == npfloat32:
        typeFactor = 1
    return typeFactor


def reverb_funDefault(musicDataArr, volumeDb, offsetSample):
    k = 10.0 ** (-abs(volumeDb) / 20.0)
    return npappend(npzeros(offsetSample, npfloat64), musicDataArr * k)[:len(musicDataArr)]


def readPcmWavData(inputWavFileName_Original):
    print("[STEPS]Reading WAV data...")
    SampleRate_Original = 0
    MusicData = []
    try:
        SampleRate_Original, MusicData = wavRead(inputWavFileName_Original)
    except Exception:
        print("[ERROR]Read WAV file failed.", file=stderr)
        exit(1)
    readWavData_left = []
    readWavData_right = []
    for item in MusicData:
        readWavData_left.append(item[0])
        readWavData_right.append(item[1])
    Larr, Rarr = nparray(readWavData_left), nparray(readWavData_right)
    if not Larr.dtype == Rarr.dtype:
        raise Exception
    else:
        readWavData_dataType = Larr.dtype
    readWavData_left = Larr * 1.0
    readWavData_right = Rarr * 1.0
    return SampleRate_Original, readWavData_left, readWavData_right, readWavData_dataType


def writePcmWavData(write_outputWavFileName, write_left, write_right, write_SampleRate, write_dataType,
                    sampleRate_Factor=0, isFilter=False, FilterDistance=7):
    print("[STEPS]Writing WAV data...")
    if not sampleRate_Factor == 0:
        write_left = write_left[range(0, len(write_left) - 1, sampleRate_Factor)]
        write_right = write_right[range(0, len(write_right) - 1, sampleRate_Factor)]
    else:
        sampleRate_Factor = 1
    if isFilter:
        print("[STEPS]Filtering...")
        medfilt(write_left, FilterDistance)
        medfilt(write_right, FilterDistance)
    MixedData = npvstack((write_left.astype(write_dataType), write_right.astype(write_dataType))).T
    try:
        wavWrite(write_outputWavFileName, int(write_SampleRate / sampleRate_Factor), MixedData)
    except Exception:
        print("[ERROR]Read WAV file failed.", file=stderr)
        exit(1)


def reverb(reverb_SampleRate, reverb_left, reverb_right, echoFunction=reverb_funDefault, reverb_numberOfEcho=3,
           reverb_maxVolumeDb=2, reverb_offsetSecond=0.22, CONSTANT_REVERB_RANGE=12.0):
    print("[STEPS]Reverb Processing...")
    offsetSample = int(reverb_offsetSecond * reverb_SampleRate)
    reverb_sum = 0.0
    factorArr = nplinspace(reverb_maxVolumeDb, reverb_maxVolumeDb * CONSTANT_REVERB_RANGE, reverb_numberOfEcho)
    for a in factorArr:
        reverb_sum += 10.0 ** (-abs(a) / 20.0)
    reverb_left = npdivide(reverb_left, (reverb_sum + 1.0))
    reverb_right = npdivide(reverb_right, (reverb_sum + 1.0))
    c = 1
    for b in factorArr:
        reverb_left += echoFunction(musicDataArr=reverb_left, volumeDb=b, offsetSample=offsetSample * c)
        reverb_right += echoFunction(musicDataArr=reverb_right, volumeDb=b, offsetSample=offsetSample * c)
        c = c + 1
    return reverb_left, reverb_right


def mixer(mixer_left, mixer_right, mixer_left_leftRate=1.0, mixer_left_rightRate=-1.0, mixer_right_leftRate=-1.0,
          mixer_right_rightRate=1.0):
    print("[STEPS]Mixing...")
    left_out = npadd(mixer_left * mixer_left_leftRate, mixer_right * mixer_left_rightRate)
    right_out = npadd(mixer_left * mixer_right_leftRate, mixer_right * mixer_right_rightRate)
    return left_out, right_out


# mode='Factor' or 'DB'
def gain(gain_left, gain_right, leftFactor=1.0, rightFactor=1.0, leftDB=0.0, rightDB=0.0, mode='Factor'):
    print("[STEPS]Gaining...")
    if mode.upper() == 'DB':
        leftFactor = 10.0 ** (-abs(leftDB) / 20.0)
        rightFactor = 10.0 ** (-abs(rightDB) / 20.0)
    return gain_left * leftFactor, gain_right * rightFactor


# This function uses Sonic library.
# Copyright (c) 2010 Bill Cox.
# File sonic.dll and its source code is licensed under the Apache 2.0 license.
def pitch(pitch_left, pitch_right, pitch_dataType, input_sampleRate, output_sampleRate, pitch_pitchFactor=0.8,
          pitch_speedFactor=1.0, lib_file=r"./sonic.dll"):
    print("[STEPS]Pitching...")
    pitch_left, pitch_right, current_datatype = changeBitRate(pitch_left, pitch_right, pitch_dataType, npint16)
    wav_len = len(pitch_left)  # length of input
    out_len = round(wav_len / pitch_pitchFactor)  # length of output.
    # Call the compiled C library(sonic.dll).
    lib = ctypes.cdll.LoadLibrary
    sonic_lib = lib(lib_file)
    wav_speech_change = sonic_lib.wavChangeSpeed
    wav_speech_change.argtypes = [ndpointer(ctypes.c_short), ndpointer(
        ctypes.c_short), ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_float]
    wav_speech_change.restypes = None
    # ----left pitch----
    resultLeft = npzeros([out_len], dtype=current_datatype)  # must int16.
    wav_speech_change(pitch_left, resultLeft, 1, input_sampleRate, pitch_left.shape[0], pitch_pitchFactor)
    # ----right pitch----
    wav_len = len(pitch_right)  # length of input
    out_len = round(wav_len / pitch_pitchFactor)  # length of output.
    # Call the compiled C library(sonic.dll).
    resultRight = npzeros([out_len], dtype=current_datatype)  # must int16.
    wav_speech_change(pitch_right, resultRight, 1, input_sampleRate, pitch_left.shape[0], pitch_pitchFactor)
    current_sampleRate, resultLeft, resultRight = resample(int(input_sampleRate / pitch_pitchFactor), output_sampleRate,
                                                           resultLeft, resultRight)
    if not pitch_speedFactor == 1.0:
        print("[STEPS]Change speed...")
        resultLeft, resultRight, current_datatype = changeBitRate(resultLeft, resultRight, resultLeft.dtype, npint16)
        # ----left speed----
        wav_len = len(resultLeft)  # length of input
        out_len = round(wav_len / pitch_speedFactor)  # length of output.
        out_resultLeft = npzeros([out_len], dtype=current_datatype)  # must int16.
        wav_speech_change(resultLeft, out_resultLeft, 1, current_sampleRate, resultLeft.shape[0], pitch_speedFactor)
        # ----right speed----
        wav_len = len(resultRight)  # length of input
        out_len = round(wav_len / pitch_speedFactor)  # length of output.
        out_resultRight = npzeros([out_len], dtype=current_datatype)  # must int16.
        wav_speech_change(resultRight, out_resultRight, 1, current_sampleRate, resultRight.shape[0], pitch_speedFactor)
    else:
        out_resultLeft, out_resultRight = resultLeft, resultRight
    return current_sampleRate, out_resultLeft, out_resultRight


def resample(original_sampleRate, current_sampleRate, resample_left, resample_right):
    print("[STEPS]Resample...")
    factor = original_sampleRate * 1.0 / current_sampleRate
    resample_left = npinterp(nparange(0, len(resample_left), factor), nparange(0, len(resample_left)), resample_left)
    resample_right = npinterp(nparange(0, len(resample_right), factor), nparange(0, len(resample_right)),
                              resample_right)
    return current_sampleRate, resample_left, resample_right


# showMode=0:only print copyright information.
# showMode=1:only print usage.
# showMode=2:print both.
def showInformation(showMode):
    if showMode == 0 or showMode == 2:
        print(NAME_AND_VERSION)
        print(COPYRIGHT_STATEMENT)
        print(r"Copyright (c) 2010 Bill Cox. (Sonic Library)")
    if showMode == 1 or showMode == 3:
        print()
        # Writing usage there.


# UNIT = second or sample.
def trim(sampleRate, trim_start, trim_end, trim_left, trim_right, trim_UNIT="second"):
    print("[STEPS]Trim...")
    startSample = trim_start
    endSample = trim_end
    if trim_UNIT.lower() == "second":
        startSample = trim_start * sampleRate
        endSample = trim_end * sampleRate
    try:
        return trim_left[startSample:endSample + 1], trim_right[startSample:endSample + 1]
    except Exception:
        print("[ERROR]Can not trim.", file=stderr)
        exit(1)


# Mode can be 'append' or 'mix'.
def joint(joint_left01, joint_left02, joint_right01, joint_right02, sampleRate01, sampleRate02, out_sampleRate,
          joint_dataType01, joint_dataType02, out_dataType, joint_MODE="append"):
    print("[STEPS]Jointing...")
    if joint_MODE.lower() == "mix":
        joint_left01 = joint_left01 / 2.0
        joint_left02 = joint_left02 / 2.0
        joint_right01 = joint_right01 / 2.0
        joint_right02 = joint_right02 / 2.0
    currentSampleRate01, joint_left01, joint_right01 = resample(sampleRate01, out_sampleRate, joint_left01,
                                                                joint_right01)
    joint_left01, joint_right01, currentDataType01 = changeBitRate(joint_left01, joint_right01, joint_dataType01,
                                                                   out_dataType)
    currentSampleRate02, joint_left02, joint_right02 = resample(sampleRate02, out_sampleRate, joint_left02,
                                                                joint_right02)
    joint_left02, joint_right02, currentDataType02 = changeBitRate(joint_left02, joint_right02, joint_dataType02,
                                                                   out_dataType)
    if joint_MODE.lower() == "mix":
        out_left = joint_left01 + joint_left02
        out_right = joint_right01 + joint_right02
    else:
        out_left = npappend(joint_left01, joint_left02)
        out_right = npappend(joint_right01, joint_right02)
    return out_sampleRate, out_left, out_right, out_dataType


def changeBitRate(cbr_left, cbr_right, original_dataType, current_dataType):
    print("[STEPS]Changing bit rate...")
    cbr_left = (cbr_left.astype(npfloat64) / getTypeFactor(original_dataType)) * getTypeFactor(current_dataType)
    cbr_right = (cbr_right.astype(npfloat64) / getTypeFactor(original_dataType)) * getTypeFactor(current_dataType)
    return cbr_left.astype(current_dataType), cbr_right.astype(current_dataType), current_dataType


# UNIT="second" or "sample", MODE="insert" or "override".
def addSilence(sampleRate, addSilence_start, addSilence_length, addSilence_left, addSilence_right, addSilence_dataType,
               addSilence_UNIT="second", addSilence_MODE="insert"):
    print("[STEPS]Adding silence...")
    startSample = addSilence_start
    lengthSample = addSilence_length
    if addSilence_UNIT.lower() == "second":
        startSample = addSilence_start * sampleRate
        lengthSample = addSilence_length * sampleRate
    try:
        if addSilence_MODE.lower() == "insert":
            addSilence_left = npinsert(addSilence_left, startSample, npzeros(lengthSample, addSilence_dataType))
            addSilence_right = npinsert(addSilence_right, startSample, npzeros(lengthSample, addSilence_dataType))
        elif addSilence_MODE.lower() == "override":
            addSilence_left = npappend(
                npappend(addSilence_left[:startSample], npzeros(lengthSample, addSilence_dataType)),
                addSilence_left[startSample + lengthSample:])
            addSilence_right = npappend(
                npappend(addSilence_right[:startSample], npzeros(lengthSample, addSilence_dataType)),
                addSilence_right[startSample + lengthSample:])
        else:
            print("[ERROR]Mode error.", file=stderr)
            exit(1)
        return addSilence_left, addSilence_right
    except Exception:
        print("[ERROR]Can not adding silence piece.", file=stderr)
        exit(1)


def isExist(inputName):
    if not os.path.exists(inputName):
        print("[ERROR]File Not Exist.", file=stderr)
        exit(1)


# This API is powered by FFmpeg.
# Copyright (c) 2000~2021 the FFmpeg developers.
# Built with gcc 11.2.0 (Rev5, Built by MSYS2 project)
# 2022/08/28
def transcoding_to_tempDir(inputFile, outputFile):
    print("[STEPS]Transcoding (Using FFmpeg)...")
    device = open(os.devnull, 'w')
    # Using run()! Popen is asynchronous!
    subprocess.run('"' + os.path.abspath(os.path.dirname(os.path.realpath(sys.path[0])) + "\\ffmpeg.exe") + '"' + " -y -i " + '"' + inputFile + '"' + " " + '"' + outputFile + '"', stderr=device)


# Main Function.
if __name__ == '__main__':
    print("[INFOS]" + COPYRIGHT_STATEMENT)
    print("[DESCR]" + NAME_AND_VERSION)
    print("[DESCR]" + BUILDER + " " + BUILD_TIME)
    filterwarnings('ignore')
    CommandLineArguments = argv  # Users should use commandline to start this program.
    parser = ArgumentParser()
    parser.add_argument("--input", "-i", required=False, help="Input file name", nargs=1, dest="inputFileName",
                        metavar="inputFile")
    parser.add_argument("--output", "-o", required=False, help="Output file name", nargs=1, dest="outputFileName",
                        metavar="outputFile")
    parser.add_argument("--output-option", required=False, help="Setup sample rate and bit rate of output file.",
                        nargs=2, dest="output_option", metavar=("sampleRate", "bitDepth"))
    parser.add_argument("--overwrite", required=False, action="store_true", dest="overwrite",
                        help="Rewriting output file if it exists.")
    parser.add_argument("--reverb", required=False,
                        help="Option of reverberation or echo, We suggest that these four parameters are 3, 2, 0.22, 12.0 (int,int,float,float).",
                        nargs=4, dest="reverb",
                        metavar=("numberOfEcho", "maxVolumeDb", "offsetSecond", "soundFieldVol"))
    parser.add_argument("--reverb-default", required=False, action="store_true",
                        help="Using default parameter to execute reverb operation.", dest="reverb_default")
    parser.add_argument("--mix", required=False,
                        help="Mix left and right channels as one. The four values should be a decimal from -1 to 1.",
                        nargs=4, dest="mix",
                        metavar=("left_leftRate", "left_rightRate", "right_leftRate", "right_rightRate"))
    parser.add_argument("--mix-default", action="store_true", required=False,
                        help="Using default parameters to mix channels (Vocal Cut).", dest="mix_default")
    parser.add_argument("--gain", required=False, help="Gain or attenuator.", nargs=3, dest="gain",
                        metavar=("valueLeft", "valueRight", "Unit"))
    parser.add_argument("--gain-default", required=False, help="Default gain operation: Gain for 3 Db.",
                        action="store_true", dest="gain_default")
    parser.add_argument("--pitch", required=False,
                        help="Modify the speed and pitch of the audio. We suggest these two values are decimals from 0.5 to 1.5",
                        nargs=2, dest="pitch", metavar=("pitchFactor", "speedFactor"))
    parser.add_argument("--pitch-default", required=False, help="Using default parameter to execute pitch.",
                        action="store_true", dest="pitch_default")
    parser.add_argument("--trim", required=False, help="Take a clip of the audio. Unit can be 'second' or 'sample'.",
                        nargs=3, dest="trim", metavar=("start", "end", "Unit"))
    parser.add_argument("--joint", required=False,
                        help="Joint(append or mix) two audio as one. Mode can be 'append' or 'mix'", nargs=3,
                        dest="joint", metavar=("inputFile01", "inputFile02", "Mode"))
    parser.add_argument("--addMute", required=False,
                        help="Add a piece of mute into the audio. Unit can be 'second' or 'sample', Mode can be 'insert' or 'override'",
                        nargs=4, dest="addSilence", metavar=("start", "length", "Unit", "Mode"))
    parser.add_argument("--filter", required=False, help="Filter the output wav data. ", action="store_true",
                        dest="filter")
    args = parser.parse_args()
    if len(sys.argv) == 1:
        print("[ERROR]Using command line argument to run this software. Add '-h' or '--help' for more help and usage.",
              file=stderr)
        exit(1)
    if args.joint is None and args.inputFileName is None:
        print("[ERROR]Please add argument of input file name.", file=stderr)
        exit(1)
    print("[STEPS]Loading Basic Modules...")

    inputWavFileName01 = r""
    inputWavFileName02 = r""
    out_SampleRate, out_DataType = 44100, npint16  # Default output parameter.
    SampleRate, dataType = 0, npint16  # Default input parameter.
    out_Format = "WAV"  # Default output format.
    # Initial left and right channels.
    left = nparray([])
    right = nparray([])

    if args.output_option is not None:
        out_SampleRate = int(args.output_option[0])
        out_BitDepth = int(args.output_option[1])
        if out_BitDepth == 16:
            out_DataType = npint16
        elif out_BitDepth == 24:
            out_DataType = npint32
        elif out_BitDepth == 32:
            out_DataType = npfloat32
        else:
            print("[ERROR]Only 32-bit, 24-bit 16-bit output bit depth are supported.", file=stderr)
            exit(1)
        if out_SampleRate > 784000 or out_SampleRate < 800:
            print("[ERROR]output file sample rate should be from 800 to 784000.", file=stderr)
            exit(1)

    if args.joint is not None:
        inputWavFileName01 = args.joint[0]
        inputWavFileName02 = args.joint[1]
        joint_mode = args.joint[2].lower()
        if (not joint_mode == "append") and (not joint_mode == "mix"):
            print("[ERROR]Argument 'MODE' should be 'append' or 'mix'.", file=stderr)
            exit(1)
        SampleRate, left, right, dataType = readPcmWavData(inputWavFileName01)
        SampleRate02, left02, right02, dataType02 = readPcmWavData(inputWavFileName02)
        SampleRate, left, right, dataType = joint(left, left02, right, right02,
                                                  SampleRate, SampleRate02, out_SampleRate, dataType, dataType02,
                                                  out_DataType, joint_MODE=joint_mode)
    else:
        inputWavFileName01 = args.inputFileName[0]
        real_proFile = inputWavFileName01
        if not inputWavFileName01.lower().endswith(".wav"):
            real_proFile = os.environ["temp"]+"\\DuYu-Audio-Processor-Temp.WAV"
            transcoding_to_tempDir(inputWavFileName01, real_proFile)
            format_index = inputWavFileName01.rfind(".")
            if not format_index == -1:
                out_Format = inputWavFileName01[format_index + 1:]
        SampleRate, left, right, dataType = readPcmWavData(real_proFile)

    outputWavFileName = inputWavFileName01 + "_OUT." + out_Format  # Default output file name.
    if args.outputFileName is not None:
        format_index = args.outputFileName[0].rfind(".")
        if not format_index == -1:
            out_Format = args.outputFileName[0][format_index + 1:]
        outputWavFileName = args.outputFileName[0]

    if args.overwrite == False and os.path.exists(outputWavFileName):
        print("[ERROR]File already exists. You can use '--overwrite' option to skip it.", file=stderr)
        exit(1)

    if args.reverb is not None:
        numberOfEcho = int(args.reverb[0])
        maxVolumeDb = int(args.reverb[1])
        offsetSecond = float(args.reverb[2])
        soundFieldVol = float(args.reverb[3])
        if numberOfEcho > 20 or numberOfEcho < 1:
            print("[ERROR]Parameter 'numberOfEcho' should be from 1 to 20.", file=stderr)
            exit(1)
        left, right = reverb(SampleRate, left, right, reverb_numberOfEcho=numberOfEcho, reverb_maxVolumeDb=maxVolumeDb,
                             reverb_offsetSecond=offsetSecond, CONSTANT_REVERB_RANGE=soundFieldVol)

    if args.reverb is None and args.reverb_default == True:
        left, right = reverb(SampleRate, left, right)

    if args.mix is not None:
        left_leftRate, left_rightRate, right_leftRate, right_rightRate = float(args.mix[0]), float(args.mix[1]), float(
            args.mix[2]), float(args.mix[3])
        left, right = mixer(left, right, left_leftRate, left_rightRate, right_leftRate, right_rightRate)

    if args.mix is None and args.mix_default == True:
        left, right = mixer(left, right)

    if args.gain is not None:
        valueLeft, valueRight, Unit = float(args.gain[0]), float(args.gain[1]), args.gain[2]
        if not Unit.lower() == "factor" and not Unit.lower() == "db":
            print("[ERROR]Unit Error.", file=stderr)
            exit(1)
        left, right = gain(left, right, valueLeft, valueRight, mode=Unit)

    if args.gain is None and args.gain_default == True:
        left, right = gain(left, right)

    if args.pitch is not None:
        pitchFactor, speedFactor = float(args.pitch[0]), float(args.pitch[1])
        if pitchFactor < 0.25:
            print("[WARNS]pitchFactor is too small.", file=stderr)
        SampleRate, left, right = pitch(left, right, dataType, SampleRate, out_SampleRate,
                                        pitch_pitchFactor=pitchFactor, pitch_speedFactor=speedFactor)

    if args.pitch is None and args.pitch_default == True:
        SampleRate, left, right = pitch(left, right, dataType, SampleRate, out_SampleRate)

    if args.addSilence is not None:
        start, length, UNIT, MODE = float(args.addSilence[0]), float(args.addSilence[1]), args.addSilence[2], \
                                    args.addSilence[3]
        UNIT = UNIT.lower()
        MODE = MODE.lower()
        if (not UNIT == "second") and (not UNIT == "sample"):
            print("[ERROR]argument 'UNIT' should be 'second' or 'sample'.", file=stderr)
            exit(1)
        if (not MODE == "insert") and (not UNIT == "override"):
            print("[ERROR]argument 'MODE' should be 'insert' or 'override'.", file=stderr)
            exit(1)
        left, right = addSilence(SampleRate, start, length, left, right, dataType, addSilence_UNIT=UNIT,
                                 addSilence_MODE=MODE)

    if args.trim is not None:
        start, end, Unit = float(args.trim[0]), float(args.trim[1]), args.trim[2]
        if (not Unit == "second") and (not Unit == "sample"):
            print("[ERROR]argument 'UNIT' should be 'second' or 'sample'.", file=stderr)
            exit(1)
        left, right = trim(SampleRate, start, end, left, right, trim_UNIT=Unit)

    if not SampleRate == out_SampleRate:
        SampleRate, left, right = resample(SampleRate, out_SampleRate, left, right)
    if not dataType == out_DataType:
        left, right, dataType = changeBitRate(left, right, dataType, out_DataType)
    real_outFile = os.environ["temp"] + "\\DuYu-Audio-Processor-Temp-Output.WAV"
    writePcmWavData(real_outFile, left, right, SampleRate, dataType, args.filter)
    transcoding_to_tempDir(real_outFile, outputWavFileName)

    # Delete temp WAV file.
    print("[STEPS]Clear temporary directory environment...")
    try:
        os.remove(os.environ["temp"] + "\\DuYu-Audio-Processor-Temp.WAV")
    except Exception:
        temp_var=None
    try:
        os.remove(os.environ["temp"] + "\\DuYu-Audio-Processor-Temp-Output.WAV")
    except Exception:
        temp_var = None

    print("[STEPS]Completed. Output File: " + outputWavFileName)

# 2020/05/02 v1.0.0
# 2021/02/17 v1.0.1
# 2021/12/31 v1.0.2
# 2022/04/14 v1.1.0
# 2022/07/05 v2.0.0
# 2022/07/11 v3.0.0
# 2022/07/11 v3.1.0
# 2022/07/12 v3.2.0
# 2022/07/18 v3.3.0
# 2022/07/19 v3.3.1
# 2022/08/28 v3.4.0
