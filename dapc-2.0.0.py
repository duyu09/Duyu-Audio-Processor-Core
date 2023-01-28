# Copyright (c) 2020~2022 Qilu University of Technology, School of Computer Science & Technology, Duyu (No.202103180009).

import numpy as np
from scipy.io import wavfile
from scipy.signal import medfilt


def reverb_funDefault(musicDataArr, volumeDb, offsetSample):
    k = 10.0 ** (-abs(volumeDb) / 20.0)
    return np.append(np.zeros(offsetSample, np.float32), musicDataArr * k)[:len(musicDataArr)]


def readPcmWavData(inputWavFileName_Original):
    SampleRate_Original, MusicData = wavfile.read(inputWavFileName_Original)
    left = []
    right = []
    for item in MusicData:
        left.append(item[0])
        right.append(item[1])
    left = np.array(left) * 1.0
    right = np.array(right) * 1.0
    return SampleRate_Original, left, right


def writePcmWavData(outputWavFileName, left, right, SampleRate, sampleRate_Factor=0, isFilter=False, FilterDistance=7):
    if not sampleRate_Factor == 0:
        left = left[range(0, len(left) - 1, sampleRate_Factor)]
        right = right[range(0, len(right) - 1, sampleRate_Factor)]
    else:
        sampleRate_Factor = 1
    if isFilter:
        medfilt(left, FilterDistance)
        medfilt(right, FilterDistance)
    MixedData = np.vstack((left.astype(np.int16), right.astype(np.int16))).T
    wavfile.write(outputWavFileName, int(SampleRate / sampleRate_Factor), MixedData)


def reverb(SampleRate, left, right, echoFunction=reverb_funDefault, numberOfEcho=3, maxvolumeDb=1, offsetSecond=0.3,
           CONSTANT_REVERB_RANGE=7.0):
    offsetSample = int(offsetSecond * SampleRate)
    sum = 0.0
    factorArr = np.linspace(maxvolumeDb, maxvolumeDb * CONSTANT_REVERB_RANGE, numberOfEcho)
    for a in factorArr:
        sum += 10.0 ** (-abs(a) / 20.0)
    left = np.divide(left, (sum + 1.0))
    right = np.divide(right, (sum + 1.0))
    c = 1
    for b in factorArr:
        left += echoFunction(musicDataArr=left, volumeDb=b, offsetSample=offsetSample * c)
        right += echoFunction(musicDataArr=right, volumeDb=b, offsetSample=offsetSample * c)
        ++c
    return left, right


def mixer(left, right, left_leftRate=1.0, left_rightRate=-1.0, right_leftRate=-1.0, right_rightRate=1.0):
    left_out = np.add(left * left_leftRate, right * left_rightRate)
    right_out = np.add(left * right_leftRate, right * right_rightRate)
    return left_out, right_out


# mode='Factor' or 'DB'
def gain(left, right, leftFactor=1.0, rightFactor=1.0, leftDB=0.0, rightDB=0.0, mode='Factor'):
    if mode.upper() == 'DB':
        leftFactor = 10.0 ** (-abs(leftDB) / 20.0)
        rightFactor = 10.0 ** (-abs(rightDB) / 20.0)
    return left * leftFactor, right * rightFactor


if __name__ == '__main__':
    # You can enter parameters via the command line instead.
    inputWavFileName = r"C:\Users\35834\Desktop\CORE.wav"
    outputWavFileName = r"C:\Users\35834\Desktop\CORE_OUT.wav"
    SampleRate, left, right = readPcmWavData(inputWavFileName)
    left, right = reverb(SampleRate, left, right)
    # left, right = mixer(left, right)
    # left, right = gain(left, right, 2, 2)
    writePcmWavData(outputWavFileName, left, right, SampleRate)

# 2020/05/02 v1.0.0
# 2021/02/17 v1.0.1
# 2021/12/31 v1.0.2
# 2022/04/14 v1.1.0
# 2022/07/05 v2.0.0
