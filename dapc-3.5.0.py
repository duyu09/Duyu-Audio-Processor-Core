"""
_*_ coding:utf-8 _*_
@Version  : 3.5.0
@Time     : 2023/01/27
@Author   : DuYu
@File     : dapc.py
@Describe : Python 3 source code of Duyu Audio Processor Core software system. (abbreviate: DAPC)
@Copyright:
            Copyright (c) 2020~2023 DuYu (No.202103180009), Faculty of Computer Science & Technology, Qilu University of Technology. (DAPC)
            Copyright (c) 2010 Bill Cox. (Sonic Library)
            Copyright (c) 2000~2021 the ffmpeg developers. (FFMPEG)
@Note     :
            1. This software relies on ffmpeg and sonic library.
            2. This source code file can be process with PyInstaller(version>=5.1).
            3. FFMPEG and sonic should be encapsulated in executable file.
            4. ffmpeg, sonic and this file should be in the same directory if you run the source code alone.
            5. We suggest that the version of Python interpreter is >=3.8.0.
"""

# Import Modules and Libraries.
import os
import sys
import ctypes
import tempfile
import platform
import subprocess
import numpy as np
from scipy import __version__ as scipyVersion
from scipy.io import wavfile
from scipy import signal
from scipy.interpolate import interp1d
from argparse import ArgumentParser
from warnings import filterwarnings

# Global const variables declare
NAME_AND_VERSION = "Duyu Audio Processor Core v3.5.0"
BUILDER = "Built and packaged with PyInstaller"
BUILD_TIME = "__January 27 2023."
COPYRIGHT_STATEMENT = "Copyright (c) 2020~2023 DuYu (No.202103180009), Faculty of Computer Science & Technology, " \
                      "Qilu University of Technology. (DAPC)"
ATTACHMENT_COPYRIGHT_STATEMENT_SONIC = "Copyright (c) 2010 Bill Cox (SONIC)"
ATTACHMENT_COPYRIGHT_STATEMENT_FFMPEG = "Copyright (c) 2000~2021 the ffmpeg developers (FFMPEG)"
ATTACHMENT_BUILD_INFORMATION_SONIC = 'Sonic Library: Built with TDM-GCC 4.9.2 64bit.'
ATTACHMENT_BUILD_INFORMATION_FFMPEG = 'FFMPEG: Built with GCC 11.2.0.'
CONSTANT_PI = np.pi  # CONSTANT_PI = 3.14159265358979323846264338327950288419716939937511
TEMPORARY_DIRECTORY = tempfile.gettempdir()
PYTHON_ENVIRONMENT = 'Python ' + platform.python_version() + ', Numpy ' + np.__version__ + ', Scipy ' + scipyVersion
OS_NAME = platform.system()
OS_ARCHITECTURE = platform.architecture()[0]
TEMP_INPUT_WAVFILE_NAME = 'DuYu-Audio-Processor-Temp.WAV'
TEMP_OUTPUT_WAVFILE_NAME = 'DuYu-Audio-Processor-Temp-Output.WAV'


def showInformation():
    """
    @Description:
                  Write information(such as copyright, explains and so on.) to stdout.
    @Parameter  :
                  No parameter.
    @Return     :
                  No return value.
    @Note       :
                  No explanation at present.
    """
    print("[INFOS]" + NAME_AND_VERSION + ' on ' + OS_NAME + ' (' + OS_ARCHITECTURE + ')')
    print("[DESCR]" + COPYRIGHT_STATEMENT)
    print("[DESCR]" + ATTACHMENT_COPYRIGHT_STATEMENT_SONIC)
    print("[DESCR]" + ATTACHMENT_COPYRIGHT_STATEMENT_FFMPEG)
    print("[INFOS]" + BUILDER + " " + BUILD_TIME)
    print("[INFOS]" + ATTACHMENT_BUILD_INFORMATION_SONIC)
    print("[INFOS]" + ATTACHMENT_BUILD_INFORMATION_FFMPEG)
    print("[INFOS]" + 'Based on ' + PYTHON_ENVIRONMENT)


# Some operation to execute when it is executed directly or imported.
showInformation()
filterwarnings('ignore')


# APIs
def getTypeFactor(gtf_dataType):
    """
    @Description:
                  Get the maximum value of the given type of data.
                  for example, the maximum value of numpy.int16 is 2^15.
    @Parameter  :
                  1. gtf_dataType (numpy.int16 or numpy.int32 or numpy.float32)
    @Return     :
                  1. The maximum value (int)
    @Note       :
                  No explanation at present.
    """
    typeFactor = 2 ** 15
    if gtf_dataType == np.int16:
        typeFactor = 2 ** 15
    elif gtf_dataType == np.int32:
        typeFactor = 2 ** 31
    elif gtf_dataType == np.float32:
        typeFactor = 1
    return typeFactor


def reverb_funDefault(musicDataArr, volumeDb, offsetSample):
    """
    @Description:
                  This function is used to be called by reverb function.
                  As we know, an echo is a superposition of multiple identical sounds
                  that are constantly decaying. The purpose of this function is to generate one
                  of the sound waves that can be directly used for stacking based on the offset
                  sample number and attenuation value.We filled the offset part with value 0.0.
    @Parameter  :
                  1. musicDataArr (ndarray): Original wave data.
                  2. volumeDb (float): The relative decibel value of attenuation.
                  3. offsetSample (int): The number of zeros added to the head of the wave data
                     array. The unit is sample.
    @Return     :
                  1. Processed wave data array. (ndarray)
    @Note       :
                  We can modify this function conveniently to change how to attenuate the wave.
    """
    k = 10.0 ** (-abs(volumeDb) / 20.0)
    return np.append(np.zeros(offsetSample, np.float32), musicDataArr * k)[:len(musicDataArr)]


def readPcmWavData(inputWavFileName_Original):
    """
    @Description:
                  This function is used to read WAV format file that encapsulate the PCM data.
                  This software only support 2-channel PCM data, and this function can split
                  these 2 channels to 2 ndarray variable.
                  We use 'scipy.io.wavfile.read()' to read wav(pcm) file.
    @Parameter  :
                  1. inputWavFileName_Original (str): Path and name of WAV file.
    @Return     :
                  1. Sample rate of pcm data. (int).
                  2. Left channel of pcm data. (ndarray).
                  3. Right channel of pcm data. (ndarray).
                  4. Data type(bit depth) of original pcm data. (numpy.int16 or int32 or float32).
    @Note       :
                  1. Main function call FFMPEG to decode other format audio file to wav file
                     firstly (Saved in temp folder as 'DuYu-Audio-Processor-Temp.WAV'),
                     and then call this function to read that file. If original file is wav
                     file, then read it directly.
    """
    print("[STEPS]Reading WAV data...")
    # SampleRate_Original = 0
    # MusicData = []
    try:
        SampleRate_Original, MusicData = wavfile.read(inputWavFileName_Original)
    except Exception:
        print("[ERROR]Read WAV file failed.", file=sys.stderr)
        sys.exit(1)
    readWavData_left = []
    readWavData_right = []
    for item in MusicData:
        readWavData_left.append(item[0])
        readWavData_right.append(item[1])
    Larr, Rarr = np.array(readWavData_left), np.array(readWavData_right)
    if not Larr.dtype == Rarr.dtype:
        raise Exception
    else:
        readWavData_dataType = Larr.dtype
    readWavData_left = Larr * 1.0
    readWavData_right = Rarr * 1.0
    return SampleRate_Original, readWavData_left, readWavData_right, readWavData_dataType


def writePcmWavData(write_outputWavFileName, write_left, write_right, write_SampleRate, write_dataType,
                    sampleRate_Factor=0, isFilter=False, FilterDistance=7):
    """
    @Description:
                  This function is used to write WAV format file that encapsulate the PCM data.
                  This software can only write 2-channel PCM wave!
                  We use 'scipy.io.wavfile.write()' to write wav file.
    @Parameter  :
                  1. write_outputWavFileName (str): Path and name of WAV file to write.
                  2. write_left (ndarray): Left channel of pcm data.
                  3. write_right (ndarray): Right channel of pcm data.
                  4. write_SampleRate (int): Sample rate of pcm data.
                  5. write_dataType (numpy.int16 or int32 or float32):
                     Data type(bit depth) of pcm data.
                  6. sampleRate_Factor (float): The multiple between of output sample rate and
                     input sample rate.It can resample, but it is very 'rough'.
                     So we don't use it after version 3.0.
                  7. isFilter (bool): Whether to use a median filter. We use scipy.signal.medfilt()
                     to filter the pcm wave data. Median filtering can effectively eliminate the
                     abrupt part of the signal.
                  8. FilterDistance (int): Filter radius.
    @Return     :
                  No return value.(Return None.)
    @Note       :
                  1. Main function write pcm data to DuYu-Audio-Processor-Temp-Output.WAV in temporary
                     directory firstly, then call FFMPEG to encode that wav file to other format.
    """
    print("[STEPS]Writing WAV data...")
    if not sampleRate_Factor == 0:
        write_left = write_left[range(0, len(write_left) - 1, sampleRate_Factor)]
        write_right = write_right[range(0, len(write_right) - 1, sampleRate_Factor)]
    else:
        sampleRate_Factor = 1
    if isFilter:
        print("[STEPS]Filtering...")
        signal.medfilt(write_left, FilterDistance)
        signal.medfilt(write_right, FilterDistance)
    MixedData = np.vstack((write_left.astype(write_dataType), write_right.astype(write_dataType))).T
    try:
        wavfile.write(write_outputWavFileName, int(write_SampleRate / sampleRate_Factor), MixedData)
    except Exception:
        print("[ERROR]Read WAV file failed.", file=sys.stderr)
        sys.exit(1)
    return None


def reverb(reverb_SampleRate, reverb_left, reverb_right, echoFunction=reverb_funDefault, reverb_numberOfEcho=3,
           reverb_maxVolumeDb=2.0, reverb_offsetSecond=0.22, CONSTANT_REVERB_RANGE=12.0):
    """
    @Description:
                  This function is used for making reverb and echo effect.
    @Parameter  :
                  1. reverb_SampleRate (int): Sample rate of input pcm data.
                  2. reverb_left (ndarray): Left channel of pcm data.
                  3. reverb_right (ndarray): Right channel of pcm data.
                  4. echoFunction (function):
                       The purpose of this function is to generate one
                     of the sound waves that can be directly used for stacking
                     based on the offset sample number and attenuation value.
                  5. reverb_numberOfEcho (int): The number of times the echo is reflected.
                  6. reverb_maxVolumeDb (float):
                       The volume at which the first reflected echo decays
                     relative to the original sound. The unit is decibel.
                  7. reverb_offsetSecond (float): The time interval between each echo.
                  8. CONSTANT_REVERB_RANGE (float):
                       The multiple of the decibel value of the last echo and
                     the decibel value of the first echo.
    @Return     :
                  1. reverb_left (ndarray): Left channel of pcm data.
                  2. reverb_right (ndarray): Right channel of pcm data.
    @Note       :
                  No note at present.
    """
    print("[STEPS]Reverb Processing...")
    offsetSample = int(reverb_offsetSecond * reverb_SampleRate)
    reverb_sum = 0.0
    factorArr = np.linspace(reverb_maxVolumeDb, reverb_maxVolumeDb * CONSTANT_REVERB_RANGE, reverb_numberOfEcho)
    for a in factorArr:
        reverb_sum += 10.0 ** (-abs(a) / 20.0)
    reverb_left = np.divide(reverb_left, (reverb_sum + 1.0))
    reverb_right = np.divide(reverb_right, (reverb_sum + 1.0))
    c = 1
    for b in factorArr:
        reverb_left += echoFunction(musicDataArr=reverb_left, volumeDb=b, offsetSample=offsetSample * c)
        reverb_right += echoFunction(musicDataArr=reverb_right, volumeDb=b, offsetSample=offsetSample * c)
        c = c + 1
    return reverb_left, reverb_right


def mixer(mixer_left, mixer_right, mixer_left_leftRate=1.0, mixer_left_rightRate=-1.0, mixer_right_leftRate=-1.0,
          mixer_right_rightRate=1.0):
    """
    @Description:
                  This function is used for mixing left and right channel.
    @Parameter  :
                  1. mixer_left (ndarray): Left channel of pcm data.
                  2. mixer_right (ndarray): Right channel of pcm data.
                  3. mixer_left_leftRate (float):
                     The proportion of the original left channel in the new left channel. (-1.0~1.0)
                  4. mixer_left_rightRate (float):
                     The proportion of the original right channel in the new left channel. (-1.0~1.0)
                  5. mixer_right_leftRate (float):
                     The proportion of the original left channel in the new right channel. (-1.0~1.0)
                  6. mixer_right_rightRate (float):
                     The proportion of the original right channel in the new right channel. (-1.0~1.0)
    @Return     :
                  1. left_out (ndarray): Left channel of pcm data.
                  2. right_out (ndarray): Right channel of pcm data.
    @Note       :
                  No note at present.
    """
    print("[STEPS]Mixing...")
    left_out = np.add(mixer_left * mixer_left_leftRate, mixer_right * mixer_left_rightRate)
    right_out = np.add(mixer_left * mixer_right_leftRate, mixer_right * mixer_right_rightRate)
    return left_out, right_out


def gain(gain_left, gain_right, leftFactor=1.0, rightFactor=1.0, leftDB=0.0, rightDB=0.0, mode='factor'):
    """
    @Description:
                  This function is used for gain or attenuate pcm data.
    @Parameter  :
                  1. gain_left (ndarray): Left channel of pcm data.
                  2. gain_right (ndarray): Right channel of pcm data.
                  3. leftFactor (float): The multiple of original left channel sample value.
                  4. rightFactor (float): The multiple of original right channel sample value.
                  5. leftDB (float): The variation of original left channel db value.
                  6. rightDB (float): The variation of original right channel db value.
                  7. mode (str): Mode of gain or attenuation. It can only be 'factor' or 'DB'.
    @Return     :
                  1. left_out (ndarray): Left channel of pcm data.
                  2. right_out (ndarray): Right channel of pcm data.
    @Note       :
                  No note at present.
    """
    print("[STEPS]Gaining...")
    if mode.upper() == 'DB':
        leftFactor = 10.0 ** (-abs(leftDB) / 20.0)
        rightFactor = 10.0 ** (-abs(rightDB) / 20.0)
    return gain_left * leftFactor, gain_right * rightFactor


def pitch(pitch_left, pitch_right, pitch_dataType, input_sampleRate, output_sampleRate, pitch_pitchFactor=0.8,
          pitch_speedFactor=1.0):
    """
    @Description:
                  This function is used for adjusting the speed and pitch.
                  This API is based on SONIC Library.
    @Parameter  :
                  1. pitch_left (ndarray): Left channel of pcm data.
                  2. pitch_right (ndarray): Right channel of pcm data.
                  3. pitch_dataType (np.int16 or np.int32 or np.float32):
                     Input pcm data bit depth.
                  4. input_sampleRate (int): Sample Rate of input pcm data.
                  5. output_sampleRate (int): Sample Rate of output pcm data.
                  6. pitch_pitchFactor (float): Multiple of sound frequency.(Do not change the speed.)
                  7. pitch_speedFactor (float): Multiple of play speed.(Do not change the pitch.)
    @Return     :
                  1. current_sampleRate (int): Sample Rate of output pcm data.
                  2. out_resultLeft (ndarray): Output left channel of pcm data.
                  3. out_resultRight (ndarray): Output right channel of pcm data.
    @Note       :
                  No note at present.
    """
    print("[STEPS]Pitching...")
    if OS_NAME.lower() == 'windows':
        lib_file = r"./sonic.dll"
    elif OS_NAME.lower() == 'darwin':
        lib_file = r"./sonic.dylib"
    else:
        lib_file = r"./sonic.so"
    pitch_left, pitch_right, current_datatype = changeBitRate(pitch_left, pitch_right, pitch_dataType, np.int16)
    wav_len = len(pitch_left)  # length of input
    out_len = round(wav_len / pitch_pitchFactor)  # length of output.
    # Call the compiled C library(sonic.dll).
    lib = ctypes.cdll.LoadLibrary
    sonic_lib = lib(lib_file)
    wav_speech_change = sonic_lib.wavChangeSpeed
    wav_speech_change.argtypes = [np.ctypeslib.ndpointer(ctypes.c_short), np.ctypeslib.ndpointer(
        ctypes.c_short), ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_float]
    wav_speech_change.restypes = None
    # ----left pitch----
    resultLeft = np.zeros([out_len], dtype=current_datatype)  # must int16.
    wav_speech_change(pitch_left, resultLeft, 1, input_sampleRate, pitch_left.shape[0], pitch_pitchFactor)
    # ----right pitch----
    wav_len = len(pitch_right)  # length of input
    out_len = round(wav_len / pitch_pitchFactor)  # length of output.
    # Call the compiled C library(sonic).
    resultRight = np.zeros([out_len], dtype=current_datatype)  # must int16.
    wav_speech_change(pitch_right, resultRight, 1, input_sampleRate, pitch_left.shape[0], pitch_pitchFactor)
    current_sampleRate, resultLeft, resultRight = resample(int(input_sampleRate / pitch_pitchFactor), output_sampleRate,
                                                           resultLeft, resultRight)
    if not pitch_speedFactor == 1.0:
        print("[STEPS]Change speed...")
        resultLeft, resultRight, current_datatype = changeBitRate(resultLeft, resultRight, resultLeft.dtype, np.int16)
        # ----left speed----
        wav_len = len(resultLeft)  # length of input
        out_len = round(wav_len / pitch_speedFactor)  # length of output.
        out_resultLeft = np.zeros([out_len], dtype=current_datatype)  # must int16.
        wav_speech_change(resultLeft, out_resultLeft, 1, current_sampleRate, resultLeft.shape[0], pitch_speedFactor)
        # ----right speed----
        wav_len = len(resultRight)  # length of input
        out_len = round(wav_len / pitch_speedFactor)  # length of output.
        out_resultRight = np.zeros([out_len], dtype=current_datatype)  # must int16.
        wav_speech_change(resultRight, out_resultRight, 1, current_sampleRate, resultRight.shape[0], pitch_speedFactor)
    else:
        out_resultLeft, out_resultRight = resultLeft, resultRight
    out_resultLeft, out_resultRight, current_datatype = changeBitRate(out_resultLeft, out_resultRight,
                                                                      current_datatype, pitch_dataType)
    return current_sampleRate, out_resultLeft, out_resultRight


# resampleAlgorithm can be "zero","nearest","slinear","linear","quadratic","cubic",etc.
def resample(original_sampleRate, current_sampleRate, resample_left, resample_right, resampleAlgorithm="default"):
    """
    @Description:
                  Resample pcm data signal.
                  We used scipy.signal.resample() (Based on FFT algorithm) or
                  several kinds of interpolation algorithms to resample.
    @Parameter  :
                  1. original_sampleRate (int): Sample Rate of input pcm data.
                  2. current_sampleRate (int): Sample Rate of output pcm data.
                  3. resample_left (ndarray): Left channel of pcm data.
                  4. resample_right (ndarray): Right channel of pcm data.
                  5. resampleAlgorithm (str): String describe of resample algorithm.
    @Return     :
                  1. current_sampleRate (int): Sample Rate of output pcm data.
                  2. resample_left (ndarray): Output left channel of pcm data.
                  3. resample_right (ndarray): Output right channel of pcm data.
    @Note       :
                  No note at present.
    """
    print("[STEPS]Resample...")
    factor = original_sampleRate * 1.0 / current_sampleRate
    if resampleAlgorithm.lower() == "default":
        """
        resample_left = np.interp(np.arange(0, len(resample_left), factor), np.arange(0, len(resample_left)),
                                 resample_left)
        resample_right = np.interp(np.arange(0, len(resample_right), factor), np.arange(0, len(resample_right)),
                                  resample_right)
        """
        resample_left = signal.resample(resample_left,
                                        int(len(resample_left) * (current_sampleRate / original_sampleRate)))
        resample_right = signal.resample(resample_right,
                                         int(len(resample_right) * (current_sampleRate / original_sampleRate)))
        return current_sampleRate, resample_left, resample_right
    else:
        resample_x01 = np.linspace(start=1, stop=len(resample_left), num=len(resample_left), dtype=np.float32)
        resample_x02 = np.linspace(start=1, stop=len(resample_left), num=int(len(resample_left) * 1.0 / factor),
                                   dtype=np.float32)
        resample_left_func = interp1d(resample_x01, resample_left, kind=resampleAlgorithm)
        resample_right_func = interp1d(resample_x01, resample_right, kind=resampleAlgorithm)
        resample_left = resample_left_func(resample_x02)
        resample_right = resample_right_func(resample_x02)
        return current_sampleRate, resample_left, resample_right


def trim(sampleRate, trim_start, trim_end, trim_left, trim_right, trim_UNIT="second"):
    """
    This function is used to trim a piece of PCM audio data.
    'trim_UNIT' means the unit of 'sampleRate'. It can be 'second' or 'sample'.
    """
    print("[STEPS]Trim...")
    startSample = trim_start
    endSample = trim_end
    if trim_UNIT.lower() == "second":
        startSample = trim_start * sampleRate
        endSample = trim_end * sampleRate
    try:
        return trim_left[startSample:endSample + 1], trim_right[startSample:endSample + 1]
    except Exception:
        print("[ERROR]Can not trim.", file=sys.stderr)
        sys.exit(1)


def joint(joint_left01, joint_left02, joint_right01, joint_right02, sampleRate01, sampleRate02, out_sampleRate,
          joint_dataType01, joint_dataType02, out_dataType, joint_MODE="append"):
    """
    This function can join two PCM ndarray data as one.
    'joint_MODE' means splicing mode, it can be 'append' or 'mix'.
    """
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
        out_left = np.append(joint_left01, joint_left02)
        out_right = np.append(joint_right01, joint_right02)
    return out_sampleRate, out_left, out_right, out_dataType


def changeBitRate(cbr_left, cbr_right, original_dataType, current_dataType):
    if original_dataType == current_dataType:
        return cbr_left, cbr_right, current_dataType
    else:
        print("[STEPS]Changing bit rate...")
        cbr_left = (cbr_left.astype(np.float32) / getTypeFactor(original_dataType)) * getTypeFactor(current_dataType)
        cbr_right = (cbr_right.astype(np.float32) / getTypeFactor(original_dataType)) * getTypeFactor(current_dataType)
        return cbr_left.astype(current_dataType), cbr_right.astype(current_dataType), current_dataType


def addSilence(sampleRate, addSilence_start, addSilence_length, addSilence_left, addSilence_right, addSilence_dataType,
               addSilence_UNIT="second", addSilence_MODE="insert"):
    """
    This function is used to add a piece of silence to the audio.
    'addSilence_UNIT' means the unit of 'addSilence_start' or 'addSilence_length'. It can be 'second' or 'sample'.
    'addSilence_MODE' can be 'insert' or 'override'.
    """
    print("[STEPS]Adding silence...")
    startSample = addSilence_start
    lengthSample = addSilence_length
    if addSilence_UNIT.lower() == "second":
        startSample = addSilence_start * sampleRate
        lengthSample = addSilence_length * sampleRate
    try:
        if addSilence_MODE.lower() == "insert":
            addSilence_left = np.insert(addSilence_left, startSample, np.zeros(lengthSample, addSilence_dataType))
            addSilence_right = np.insert(addSilence_right, startSample, np.zeros(lengthSample, addSilence_dataType))
        elif addSilence_MODE.lower() == "override":
            addSilence_left = np.append(
                np.append(addSilence_left[:startSample], np.zeros(lengthSample, addSilence_dataType)),
                addSilence_left[startSample + lengthSample:])
            addSilence_right = np.append(
                np.append(addSilence_right[:startSample], np.zeros(lengthSample, addSilence_dataType)),
                addSilence_right[startSample + lengthSample:])
        else:
            print("[ERROR]Mode error.", file=sys.stderr)
            sys.exit(1)
        return addSilence_left, addSilence_right
    except Exception:
        print("[ERROR]Can not adding silence piece.", file=sys.stderr)
        sys.exit(1)


def surround3d(s3d_left, s3d_right, s3d_sampleRate, s3d_inputDataType, s3d_outputDataType,
               s3d_roundLength=12.75, s3d_minimumRangeFactor=0.3, s3d_phase=0.0, s3d_roundLengthUnit="second"):
    """
    This function is used to implement the '3D surround' effect.
    1. 's3d_roundLength' means the length of a whole surround circle.
    2. 's3d_minimumRangeFactor' means the factor of the original volume and
       the minimum volume of the whole circle.
       (0.0<=s3d_minimumRangeFactor<=1.0)
    3. 's3d_phase' means the start position (angle degree) of the surround circle.
       (0.0<=s3d_phase<360.0)
    4. 's3d_roundLengthUnit' means the unit of 's3d_roundLength'.
       (s3d_roundLengthUnit can be "second" or "sample")
    """
    print("[STEPS]Surround 3D processing...")
    s3d_left, s3d_right, s3d_dt = changeBitRate(s3d_left, s3d_right, s3d_inputDataType, np.float32)
    if s3d_roundLengthUnit == "second":
        s3d_roundLength = s3d_roundLength * s3d_sampleRate
    s3d_roundLength = int(s3d_roundLength)
    arrUnitRound = (np.sin(np.arange(0, s3d_roundLength) * (2.0 * CONSTANT_PI / s3d_roundLength)) + 1.0) / 2.0 * \
                   (1.0 - s3d_minimumRangeFactor) + s3d_minimumRangeFactor
    arrRound = np.tile(arrUnitRound, int(len(s3d_left) / s3d_roundLength) + 1)
    arrRoundLeft = np.append(arrUnitRound[:int(len(arrUnitRound) * (s3d_phase / 360.0))], arrRound)[:len(s3d_left)]
    arrRoundRight = 1.0 - (arrRoundLeft - s3d_minimumRangeFactor)
    s3d_left, s3d_right, s3d_dt = changeBitRate(s3d_left * arrRoundLeft, s3d_right * arrRoundRight, s3d_dt,
                                                s3d_outputDataType)
    return s3d_left, s3d_right


def fadeInOut(fio_left, fio_right, fio_sampleRate, fio_inputDataType, fio_outputDataType,
              fio_fadeLength=18.5, fio_fadeMode='sin', fio_fadeLengthUnit='second', fio_mode='out'):
    """
    This function is used to implement the fade in and fade out effect.
    'fio_fadeMode' can be 'sin','cos','liner'.
    'fio_mode' can be 'out','in','both'.
    'fio_fadeLengthUnit' can be 'second','sample'.
    """
    print("[STEPS]Fade in or out processing...")
    fio_left, fio_right, fio_dt = changeBitRate(fio_left, fio_right, fio_inputDataType, np.float32)
    if fio_fadeLengthUnit == "second":
        fio_fadeLength = fio_fadeLength * fio_sampleRate
    fio_fadeLength = int(fio_fadeLength)
    if fio_mode == 'out' or fio_mode == 'both':
        arr_leftHead = fio_left[:-fio_fadeLength]
        arr_rightHead = fio_right[:-fio_fadeLength]
        arr_leftTail = fio_left[-fio_fadeLength:]
        arr_rightTail = fio_right[-fio_fadeLength:]
        arr = np.array([])
        if fio_fadeMode == 'sin':
            arr = np.sin(np.arange(fio_fadeLength) * (0.5 * CONSTANT_PI / fio_fadeLength) + CONSTANT_PI * 0.5)
        elif fio_fadeMode == 'cos':
            arr = np.cos(np.arange(fio_fadeLength) * (0.5 * CONSTANT_PI / fio_fadeLength))
        elif fio_fadeMode == 'liner':
            arr = 1.0 - np.arange(fio_fadeLength) / fio_fadeLength
        fio_left = np.append(arr_leftHead, arr_leftTail * arr)
        fio_right = np.append(arr_rightHead, arr_rightTail * arr)
    elif fio_mode == 'in' or fio_mode == 'both':
        arr_leftHead = fio_left[:fio_fadeLength]
        arr_rightHead = fio_right[:fio_fadeLength]
        arr_leftTail = fio_left[fio_fadeLength:]
        arr_rightTail = fio_right[fio_fadeLength:]
        arr = np.array([])
        if fio_fadeMode == 'sin':
            arr = np.sin(np.arange(fio_fadeLength) * (0.5 * CONSTANT_PI / fio_fadeLength))
        elif fio_fadeMode == 'cos':
            arr = np.cos(np.arange(fio_fadeLength) * (0.5 * CONSTANT_PI / fio_fadeLength) + CONSTANT_PI)
        elif fio_fadeMode == 'liner':
            arr = np.arange(fio_fadeLength) / fio_fadeLength
        fio_left = np.append(arr_leftHead * arr, arr_leftTail)
        fio_right = np.append(arr_rightHead * arr, arr_rightTail)
    fio_left, fio_right, fio_dt = changeBitRate(fio_left, fio_right, fio_dt, fio_outputDataType)
    return fio_left, fio_right


def fftFilter(fft_left, fft_right, fft_sampleRate, fft_inputDataType, fft_outputDataType, fft_frequency=20000,
              fft_frequencyRange=180):
    """
    This function is used to cut off part of the PCM data of some frequency.
    It can cut off the frequency in the range of [fft_frequency-fft_frequencyRange, fft_frequency+fft_frequencyRange]
    It bases on 'np.fft' package.
    """
    print("[STEPS]FFT Filter...")
    fft_left, fft_right, fft_dt = changeBitRate(fft_left, fft_right, fft_inputDataType, np.float32)
    afterFft_left = np.fft.fft(fft_left)
    afterFft_right = np.fft.fft(fft_right)
    fft_index_left_low = int(max(fft_frequency - fft_frequencyRange, 0) / fft_sampleRate * len(fft_left))
    fft_index_left_high = int(min(fft_frequency + fft_frequencyRange, fft_sampleRate) / fft_sampleRate * len(fft_left))
    fft_index_right_low = int(max(fft_frequency - fft_frequencyRange, 0) / fft_sampleRate * len(fft_right))
    fft_index_right_high = int(
        min(fft_frequency + fft_frequencyRange, fft_sampleRate) / fft_sampleRate * len(fft_right))
    for i in range(fft_index_left_low, fft_index_left_high):
        afterFft_left[i] = 0.0 + 0.0j
    for i in range(fft_index_right_low, fft_index_right_high):
        afterFft_right[i] = 0.0 + 0.0j
    afterFft_left, afterFft_right = afterFft_left[:int(len(afterFft_left) / 2)], afterFft_right[
                                                                                 :int(len(afterFft_right) / 2)]
    fft_left = np.fft.irfft(afterFft_left)
    fft_right = np.fft.irfft(afterFft_right)
    fft_left, fft_right, fft_dt = changeBitRate(fft_left, fft_right, fft_dt, fft_outputDataType)
    return fft_left, fft_right, fft_dt


def isExist(inputName):
    if not os.path.exists(inputName):
        print("[ERROR]File Not Exist.", file=sys.stderr)
        sys.exit(1)


def transcoding_to_tempDir(inputFile, outputFile):
    """
    @Description:
                  This function is used to call FFMPEG to decode or encode audio data.
                  This API is powered by FFMPEG.
    @Parameter  :
                  1. write_outputWavFileName (str): Path and name of WAV file to write.
                  2. write_left (ndarray): Left channel of pcm data.

    @Return     :
                  No return value.(Return None.)
    @Note       :
                  Executable file FFMPEG should be in the same directory with this script file,
                  or packaged in root directory in executable file with pyinstaller.
    """
    print("[STEPS]Transcoding (Using FFMPEG)...")
    device = open(os.devnull, 'w')
    # Run the script file separately:
    subprocess.call(os.path.join('.', 'ffmpeg') + " -y -i " + '"' + inputFile + '"' + " " + '"' + outputFile + '"',
                    stderr=device, shell=True)
    # Run binary executable file:

    # Using run()! Popen is asynchronous!
    # subprocess.call(os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(sys.path[0]))), 'ffmpeg') + " -y -i " + '"' + inputFile + '"' + " " + '"' + outputFile + '"', stderr=device)


# Main Function.
if __name__ == '__main__':
    filterwarnings('ignore')
    CommandLineArguments = sys.argv  # Users should use commandline to start this program.
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
    parser.add_argument("--surround3d", required=False,
                        help="Generate 3D surround sound effects. 'roundLengthUnit' can be second or sample.",
                        nargs=4, dest="surround3d",
                        metavar=("roundLength", "minimumRangeFactor", "phase", "roundLengthUnit"))
    parser.add_argument("--surround3d-default", required=False, help="Using default parameter to execute surround 3D.",
                        action="store_true", dest="surround3d_default")
    parser.add_argument("--fade-default", required=False, help="Using default parameter to execute fade in and out.",
                        action="store_true", dest="fade_default")
    parser.add_argument("--fade", required=False, help="Execute fade in and out. "
                                                       "'fadeLength' means the second of a whole surround circle. "
                                                       "'mode' can be 'out' or 'in' or 'both'",
                        dest="fade", nargs=2, metavar=("fadeLength", "mode"))
    parser.add_argument("--fftFilter-default", required=False, help="Using default parameter to execute fft filter.",
                        action="store_true", dest="fftFilter_default")
    parser.add_argument("--fftFilter", required=False, help="Execute FFT Filter, Cut off some frequency. "
                                                            "Cut off the frequency in the range of "
                                                            "[fft_frequency-fft_frequencyRange, fft_frequency+fft_frequencyRange]",
                        dest="fftFilter", nargs=2, metavar=("frequency", "range"))
    args = parser.parse_args()
    if len(sys.argv) == 1:
        print("[ERROR]Using command line argument to run this software. Add '-h' or '--help' for more help and usage.",
              file=sys.stderr)
        sys.exit(1)
    if args.joint is None and args.inputFileName is None:
        print("[ERROR]Please add argument of input file name.", file=sys.stderr)
        sys.exit(1)
    print("[STEPS]Loading Basic Modules...")

    inputWavFileName01 = r""
    inputWavFileName02 = r""
    out_SampleRate, out_DataType = 44100, np.int16  # Default output parameter.
    SampleRate, dataType = 0, np.int16  # Default input parameter.
    out_Format = "WAV"  # Default output format.
    # Initial left and right channels.
    left = np.array([])
    right = np.array([])

    if args.output_option is not None:
        out_SampleRate = int(args.output_option[0])
        out_BitDepth = int(args.output_option[1])
        if out_BitDepth == 16:
            out_DataType = np.int16
        elif out_BitDepth == 24:
            out_DataType = np.int32
        elif out_BitDepth == 32:
            out_DataType = np.float32
        else:
            print("[ERROR]Only 32-bit, 24-bit 16-bit output bit depth are supported.", file=sys.stderr)
            sys.exit(1)
        if out_SampleRate > 784000 or out_SampleRate < 800:
            print("[ERROR]output file sample rate should be from 800 to 784000.", file=sys.stderr)
            sys.exit(1)

    if args.joint is not None:
        inputWavFileName01 = args.joint[0]
        inputWavFileName02 = args.joint[1]
        joint_mode = args.joint[2].lower()
        if (not joint_mode == "append") and (not joint_mode == "mix"):
            print("[ERROR]Argument 'MODE' should be 'append' or 'mix'.", file=sys.stderr)
            sys.exit(1)
        SampleRate, left, right, dataType = readPcmWavData(inputWavFileName01)
        SampleRate02, left02, right02, dataType02 = readPcmWavData(inputWavFileName02)
        SampleRate, left, right, dataType = joint(left, left02, right, right02,
                                                  SampleRate, SampleRate02, out_SampleRate, dataType, dataType02,
                                                  out_DataType, joint_MODE=joint_mode)
    else:
        inputWavFileName01 = args.inputFileName[0]
        real_proFile = inputWavFileName01
        if not inputWavFileName01.lower().endswith(".wav"):
            real_proFile = os.path.join(TEMPORARY_DIRECTORY, TEMP_INPUT_WAVFILE_NAME)
            transcoding_to_tempDir(inputWavFileName01, real_proFile)
            format_index = inputWavFileName01.rfind(".")
            if not format_index == -1:
                out_Format = inputWavFileName01[format_index + 1:]
        SampleRate, left, right, dataType = readPcmWavData(real_proFile)
        if args.output_option is None:
            out_SampleRate, out_DataType = SampleRate, dataType
        left, right, dataType = changeBitRate(left, right, dataType, np.float32)

    outputWavFileName = inputWavFileName01[
                        :inputWavFileName01.rfind(".")] + "_OUT." + out_Format  # Default output file name.
    if args.outputFileName is not None:
        format_index = args.outputFileName[0].rfind(".")
        if not format_index == -1:
            out_Format = args.outputFileName[0][format_index + 1:]
        outputWavFileName = args.outputFileName[0]

    if args.overwrite == False and os.path.exists(outputWavFileName):
        print("[ERROR]File already exists. You can use '--overwrite' option to skip it.", file=sys.stderr)
        sys.exit(1)

    if args.reverb is not None:
        numberOfEcho = int(args.reverb[0])
        maxVolumeDb = int(args.reverb[1])
        offsetSecond = float(args.reverb[2])
        soundFieldVol = float(args.reverb[3])
        if numberOfEcho > 20 or numberOfEcho < 1:
            print("[ERROR]Parameter 'numberOfEcho' should be from 1 to 20.", file=sys.stderr)
            sys.exit(1)
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
            print("[ERROR]Unit Error.", file=sys.stderr)
            sys.exit(1)
        left, right = gain(left, right, valueLeft, valueRight, mode=Unit)

    if args.gain is None and args.gain_default == True:
        left, right = gain(left, right)

    if args.pitch is not None:
        pitchFactor, speedFactor = float(args.pitch[0]), float(args.pitch[1])
        if pitchFactor < 0.25:
            print("[WARNS]pitchFactor is too small.", file=sys.stderr)
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
            print("[ERROR]argument 'UNIT' should be 'second' or 'sample'.", file=sys.stderr)
            sys.exit(1)
        if (not MODE == "insert") and (not UNIT == "override"):
            print("[ERROR]argument 'MODE' should be 'insert' or 'override'.", file=sys.stderr)
            sys.exit(1)
        left, right = addSilence(SampleRate, start, length, left, right, dataType, addSilence_UNIT=UNIT,
                                 addSilence_MODE=MODE)

    if args.trim is not None:
        start, end, Unit = float(args.trim[0]), float(args.trim[1]), args.trim[2]
        if (not Unit == "second") and (not Unit == "sample"):
            print("[ERROR]argument 'UNIT' should be 'second' or 'sample'.", file=sys.stderr)
            sys.exit(1)
        left, right = trim(SampleRate, start, end, left, right, trim_UNIT=Unit)

    if args.surround3d is not None:
        roundLength, minimumRangeFactor, phase, roundLengthUnit = float(args.surround3d[0]), float(args.surround3d[1]), \
                                                                  float(args.surround3d[2]), args.surround3d[3].lower()
        if minimumRangeFactor < 0.0 or minimumRangeFactor > 1.0:
            print("[ERROR]minimumRangeFactor should be in 0.0~1.0", file=sys.stderr)
        if phase < 0.0 or phase > 360.0:
            print("[ERROR]Phase should be in 0.0~360.0", file=sys.stderr)
        left, right = surround3d(left, right, SampleRate, dataType, dataType, s3d_roundLength=roundLength,
                                 s3d_minimumRangeFactor=minimumRangeFactor, s3d_phase=phase,
                                 s3d_roundLengthUnit=roundLengthUnit)
    if args.surround3d is None and args.surround3d_default == True:
        left, right = surround3d(left, right, SampleRate, dataType, dataType)
    if args.fade_default and args.fade is None:
        left, right = fadeInOut(left, right, SampleRate, dataType, dataType)
    if args.fade is not None:
        left, right = fadeInOut(left, right, SampleRate, dataType, dataType,
                                fio_fadeLength=float(args.fade[0]), fio_mode=str(args.fade[1]))
    if args.fftFilter_default and args.fftFilter is None:
        left, right, dataType = fftFilter(left, right, SampleRate, dataType, dataType)
    if args.fftFilter is not None:
        left, right, dataType = fftFilter(left, right, SampleRate, dataType, dataType,
                                          fft_frequency=int(args.fftFilter[0]),
                                          fft_frequencyRange=int(args.fftFilter[1]))

    # The last step.
    if not SampleRate == out_SampleRate:
        SampleRate, left, right = resample(SampleRate, out_SampleRate, left, right)
    if not dataType == out_DataType:
        left, right, dataType = changeBitRate(left, right, dataType, out_DataType)
    real_outFile = os.path.join(TEMPORARY_DIRECTORY, TEMP_OUTPUT_WAVFILE_NAME)
    writePcmWavData(real_outFile, left, right, SampleRate, dataType, isFilter=args.filter)
    transcoding_to_tempDir(real_outFile, outputWavFileName)

    # Delete temporary WAV file.
    print("[STEPS]Clear temporary directory environment...")
    try:
        os.remove(os.path.join(TEMPORARY_DIRECTORY, TEMP_INPUT_WAVFILE_NAME))
    except Exception:
        pass
    try:
        os.remove(os.path.join(TEMPORARY_DIRECTORY, TEMP_OUTPUT_WAVFILE_NAME))
    except Exception:
        pass

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
# 2023/01/27 v3.5.0
