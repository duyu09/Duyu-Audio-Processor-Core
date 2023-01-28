# Copyright (c) 2020~2021 Duyu.
import numpy as np
from scipy.io import wavfile
from scipy import signal

def splitChannel(srcMusicFile):
   # read wav file
   sampleRate, musicData = wavfile.read(srcMusicFile)

   left = []
   right = []
   for item in musicData:
       left.append(item[0])
       right.append(item[1])

   mixed_data = np.array(left) - np.array(right)
   wavfile.write(srcMusicFile + "_output.wav", sampleRate, mixed_data)

splitChannel("1.wav")

# 2020/05/02 v1.0.0
# 2021/02/17 v1.0.1
# 2021/12/31 v1.0.2
