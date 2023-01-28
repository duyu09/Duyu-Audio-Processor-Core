# Copyright (c) 2020~2022 Duyu.
import numpy as np
from scipy.io import wavfile
from scipy import signal

def echo(SrcMusic,Second=0.1,dB=10):
   SampleRate, MusicData = wavfile.read(SrcMusic)
   left = []
   right = []
   for item in MusicData:
       left.append(item[0])
       right.append(item[1])
   left1=np.array(left)
   right1=np.array(right)
   left2=left1.copy()
   right2=right1.copy()
   headSample=int(Second*SampleRate)
   n=10.0**(-abs(dB)/20.0)
   b=0
   for a in range(headSample+1,len(left)):
       left1[a]=(left1[a]*1.0+n*left2[b]).astype(left1.dtype)
       right1[a]=right1[a]*1.0+n*right2[b].astype(right1.dtype)
       b=b+1
   MixedData = np.vstack((left1,right1)).T
   wavfile.write(SrcMusic + "_output.wav", SampleRate, MixedData)

echo(r".\CORE.wav",0.08,3)
# echo(r"C:\Users\86188\Desktop\CORE(HIGH_FRENQUENCE).wav",0.05,1)

# 2020/05/02 v1.0.0
# 2021/02/17 v1.0.1
# 2021/12/31 v1.0.2
# 2022/04/14 v1.1.0
