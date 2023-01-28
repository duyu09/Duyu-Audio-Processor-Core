# Copyright (c) 2020~2022 Qilu University of Technology, School of Computer Science & Technology, Duyu (No.202103180009).
# Copyright (c) 2010 Bill Cox. (Sonic Library)
# Version 3.2


import ctypes
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
from pylab import plot
from pylab import figure
from pylab import xlabel
from pylab import ylabel
from pylab import show
from warnings import filterwarnings


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


def reverb(reverb_SampleRate, reverb_left, reverb_right, echoFunction=reverb_funDefault, numberOfEcho=3,
           maxVolumeDb=2, offsetSecond=0.22, CONSTANT_REVERB_RANGE=12.0):
    print("[STEPS]Reverb Processing...")
    offsetSample = int(offsetSecond * reverb_SampleRate)
    reverb_sum = 0.0
    factorArr = nplinspace(maxVolumeDb, maxVolumeDb * CONSTANT_REVERB_RANGE, numberOfEcho)
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


def mixer(mixer_left, mixer_right, left_leftRate=1.0, left_rightRate=-1.0, right_leftRate=-1.0, right_rightRate=1.0):
    print("[STEPS]Mixing...")
    left_out = npadd(mixer_left * left_leftRate, mixer_right * left_rightRate)
    right_out = npadd(mixer_left * right_leftRate, mixer_right * right_rightRate)
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
def pitch(pitch_left, pitch_right, pitch_dataType, input_sampleRate, output_sampleRate, pitchFactor=0.8,
          speedFactor=1.0, lib_file=r"./sonic.dll"):
    print("[STEPS]Pitching...")
    pitch_left, pitch_right, current_datatype = changeBitRate(pitch_left, pitch_right, pitch_dataType, npint16)
    wav_len = len(pitch_left)  # length of input
    out_len = round(wav_len / pitchFactor)  # length of output.
    # Call the compiled C library(sonic.dll).
    lib = ctypes.cdll.LoadLibrary
    sonic_lib = lib(lib_file)
    wav_speech_change = sonic_lib.wavChangeSpeed
    wav_speech_change.argtypes = [ndpointer(ctypes.c_short), ndpointer(
        ctypes.c_short), ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_float]
    wav_speech_change.restypes = None
    # ----left pitch----
    resultLeft = npzeros([out_len], dtype=current_datatype)  # must int16.
    wav_speech_change(pitch_left, resultLeft, 1, input_sampleRate, pitch_left.shape[0], pitchFactor)
    # ----right pitch----
    wav_len = len(pitch_right)  # length of input
    out_len = round(wav_len / pitchFactor)  # length of output.
    # Call the compiled C library(sonic.dll).
    resultRight = npzeros([out_len], dtype=current_datatype)  # must int16.
    wav_speech_change(pitch_right, resultRight, 1, input_sampleRate, pitch_left.shape[0], pitchFactor)
    current_sampleRate, resultLeft, resultRight = resample(int(input_sampleRate / pitchFactor), output_sampleRate, resultLeft, resultRight)
    if not speedFactor == 1.0:
        print("[STEPS]Change speed...")
        resultLeft, resultRight, current_datatype = changeBitRate(resultLeft, resultRight, resultLeft.dtype, npint16)
        # ----left speed----
        wav_len = len(resultLeft)  # length of input
        out_len = round(wav_len / speedFactor)  # length of output.
        out_resultLeft = npzeros([out_len], dtype=current_datatype)  # must int16.
        wav_speech_change(resultLeft, out_resultLeft, 1, current_sampleRate, resultLeft.shape[0], speedFactor)
        # ----right speed----
        wav_len = len(resultRight)  # length of input
        out_len = round(wav_len / speedFactor)  # length of output.
        out_resultRight = npzeros([out_len], dtype=current_datatype)  # must int16.
        wav_speech_change(resultRight, out_resultRight, 1, current_sampleRate, resultRight.shape[0], speedFactor)
    else:
        out_resultLeft, out_resultRight = resultLeft, resultRight
    return current_sampleRate, out_resultLeft, out_resultRight


def showFrequencyDomainWave(sampleRate, waveData, fWave_dataType, fftSize=2048, offset=0, figSize=(8, 4)):
    print("[STEPS]Analysing Frequency Domain Wave...")
    typeFactor = getTypeFactor(fWave_dataType)
    xs = waveData[offset:fftSize + offset] / typeFactor
    xf = nprfft(xs) / fftSize
    freq = nplinspace(0, int(sampleRate / 2), int(fftSize / 2) + 1)
    xfp = 20 * nplog10(npclip(npabs(xf), 1e-20, 1e100))
    figure("Frequency Domain Wave", figsize=figSize)
    plot(freq, xfp)
    xlabel("Frequency(Hz)")
    ylabel("Volume(dB)")
    show()


def showTimeDomainWave(sampleRate, waveData, timeWave_dataType, figSize=(8, 4), xLabelUnit='second'):
    print("[STEPS]Analysing Time Domain Wave...")
    length = len(waveData)
    typeFactor = getTypeFactor(timeWave_dataType)
    figure("Time Domain Wave", figsize=figSize)
    if xLabelUnit.lower() == 'second':
        plot(nparray(range(0, length)) / sampleRate, waveData / typeFactor)
        xlabel("Time(s)")
    else:
        plot(range(0, length), waveData / typeFactor)
        xlabel("Samples")
    ylabel("Amplitude")
    show()


def resample(original_sampleRate, current_sampleRate, resample_left, resample_right):
    print("[STEPS]Resample...")
    factor = original_sampleRate * 1.0/current_sampleRate
    resample_left = npinterp(nparange(0, len(resample_left), factor), nparange(0, len(resample_left)), resample_left)
    resample_right = npinterp(nparange(0, len(resample_right), factor), nparange(0, len(resample_right)), resample_right)
    return current_sampleRate, resample_left, resample_right


# showMode=0:only print copyright information.
# showMode=1:only print usage.
# showMode=2:print both.
def showInformation(showMode):
    if showMode == 0 or showMode == 2:
        print(r"DuyuAudioProcessor-CORE v3.2")
        print(r"Copyright (c) 2020~2022 Qilu University of Technology, School of Computer Science & Technology, "
              r"Duyu (No.202103180009).")
        print(r"Copyright (c) 2010 Bill Cox. (Sonic Library)")
    if showMode == 1 or showMode == 3:
        print(r"Usage of DuyuPCMprocessor-CORE")
        print()
        # Writing usage there.


# UNIT = second or sample.
def trim(sampleRate, start, end, trim_left, trim_right, UNIT="second"):
    print("[STEPS]Trim...")
    startSample = start
    endSample = end
    if UNIT.lower() == "second":
        startSample = start * sampleRate
        endSample = end * sampleRate
    try:
        return trim_left[startSample:endSample + 1], trim_right[startSample:endSample + 1]
    except Exception:
        print("[ERROR]Can not trim.", file=stderr)
        exit(1)


def joint(joint_left01, joint_left02, joint_right01, joint_right02, sampleRate01, sampleRate02, out_sampleRate,
          joint_dataType01, joint_dataType02, out_dataType):
    print("[STEPS]Jointing...")
    currentSampleRate01, joint_left01, joint_right01 = resample(sampleRate01, out_sampleRate, joint_left01, joint_right01)
    joint_left01, joint_right01, currentDataType01 = changeBitRate(joint_left01, joint_right01, joint_dataType01, out_dataType)
    currentSampleRate02, joint_left02, joint_right02 = resample(sampleRate02, out_sampleRate, joint_left02, joint_right02)
    joint_left02, joint_right02, currentDataType02 = changeBitRate(joint_left02, joint_right02, joint_dataType02, out_dataType)
    out_left = npappend(joint_left01, joint_left02)
    out_right = npappend(joint_right01, joint_right02)
    return out_sampleRate, out_left, out_right, out_dataType


def changeBitRate(cbr_left, cbr_right, original_dataType, current_dataType):
    print("[STEPS]Changing bit rate...")
    cbr_left = (cbr_left.astype(npfloat64) / getTypeFactor(original_dataType)) * getTypeFactor(current_dataType)
    cbr_right = (cbr_right.astype(npfloat64) / getTypeFactor(original_dataType)) * getTypeFactor(current_dataType)
    return cbr_left.astype(current_dataType), cbr_right.astype(current_dataType), current_dataType


# UNIT="second" or "sample", MODE="insert" or "override".
def addSilence(sampleRate, start, length, addSilence_left, addSilence_right, addSilence_dataType,
               UNIT="second", MODE="insert"):
    print("[STEPS]Adding silence...")
    startSample = start
    lengthSample = length
    if UNIT.lower() == "second":
        startSample = start * sampleRate
        lengthSample = length * sampleRate
    try:
        if MODE.lower() == "insert":
            addSilence_left = npinsert(addSilence_left, startSample, npzeros(lengthSample, addSilence_dataType))
            addSilence_right = npinsert(addSilence_right, startSample, npzeros(lengthSample, addSilence_dataType))
        elif MODE.lower() == "override":
            addSilence_left = npappend(npappend(addSilence_left[:startSample], npzeros(lengthSample, addSilence_dataType)), addSilence_left[startSample+lengthSample:])
            addSilence_right = npappend(npappend(addSilence_right[:startSample], npzeros(lengthSample, addSilence_dataType)), addSilence_right[startSample+lengthSample:])
        else:
            print("[ERROR]Mode error.", file=stderr)
            exit(1)
        return addSilence_left, addSilence_right
    except Exception:
        print("[ERROR]Can not adding silence piece.", file=stderr)
        exit(1)


# Test & Demo.
if __name__ == '__main__':
    print("[STEPS]Loading Basic Modules...")
    filterwarnings('ignore')
    CommandLineArguments = argv  # We can use commandline to start this program.
    inputWavFileName01 = r"C:\Users\35834\Desktop\CORE.wav"
    # inputWavFileName02 = r"C:\Users\35834\Desktop\CORE02.wav"
    outputWavFileName = r"C:\Users\35834\Desktop\CORE_OUT.wav"
    SampleRate, left, right, dataType = readPcmWavData(inputWavFileName01)
    # SampleRate02, left02, right02, dataType02 = readPcmWavData(inputWavFileName02)
    # SampleRate, left, right, dataType = joint(left, left02, right, right02, SampleRate, SampleRate02, 90000, dataType, dataType02, npint32)
    # left, right = reverb(SampleRate, left, right)
    # left, right = mixer(left, right)
    # left, right = gain(left, right, 2, 2)
    # showFrequencyDomainWave(SampleRate, left, dataType, offset=500000)
    # showTimeDomainWave(SampleRate, left, dataType)
    # SampleRate, left, right = resample(SampleRate, 44100, left, right)
    # showInformation(0)
    # left, right = trim(SampleRate, 20, 30, left, right)
    # left, right, dataType = changeBitRate(left, right, dataType, npfloat32)
    # SampleRate, left, right = pitch(left, right, dataType, SampleRate, 44100)
    # left, right = addSilence(SampleRate, 5, 6, left, right, dataType, MODE="override")
    writePcmWavData(outputWavFileName, left, right, SampleRate, dataType, isFilter=True)
    print("[STEPS]Completed.")

# 2020/05/02 v1.0.0
# 2021/02/17 v1.0.1
# 2021/12/31 v1.0.2
# 2022/04/14 v1.1.0
# 2022/07/05 v2.0.0
# 2022/07/11 v3.0.0
# 2022/07/11 v3.1.0
# 2022/07/12 v3.2.0
