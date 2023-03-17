# Duyu Audio Processor Core

## Software Introduction  软件简介

软件全称：Duyu Audio Processor Core

软件简称：DAPC (或 dapc)

软件中文全称：DuYu音频处理器核心程序

软件功能：简单处理主流格式的音频文件，诸如：变换采样格式，增益衰减，声道混缩，混响回声效果，调整速度与音调，剪切，混合，添加静音，淡入淡出效果，FFT滤波，3D环绕声音效等。

项目开工时间：2020年05月02日

项目中使用的第三方软件：Sonic Library, FFmpeg.

## Copyright Statement  版权声明

Python source code of Duyu Audio Processor Core (DAPC) software system:

Copyright &copy; 2020~2023 DuYu (No.202103180009), Faculty of Computer Science & Technology, Qilu University of Technology.
            

All rights reserved.

All of the DAPC source code is licensed under the Apache 2.0 license.

<br>

DAPC软件系统Python源代码：

&copy; 2020~2023 齐鲁工业大学 计算机科学与技术学部 杜宇。 

保留所有权利。

DAPC软件系统的所有Python源代码都是在Apache 2.0协议的许可下授权的。

----

C source code of Sonic Library:

Copyright &copy; 2010 Bill Cox.

All of the Sonic library source code is licensed under the Apache 2.0 license.

<br>

Sonic动态链接库C语言源代码：

&copy; 2010 比尔-考克斯

Sonic动态链接库的所有C语言源代码都是在Apache 2.0协议的许可下授权的。

----

C source code of FFmpeg:

Copyright &copy; 2000~2021 the FFmpeg developers.

All of the ffmpeg source code is licensed under the LGPL/GPL.

<br>

FFmpeg音视频转码单元C语言源代码：

&copy; 2000~2021 FFmpeg全体开发者

FFmpeg的所有C语言源代码都是在LGPL/GPL协议的许可下授权的。

----

## Software Update Logs  软件更新日志

### Update on May. 2nd, 2020

基于实现音频消除人声的目的，项目开工。

v1.0.0版本：实现了基于声道“相抵消”方法的音频消除人声和保留伴奏的功能。

### Update on Dec. 31st, 2021

v1.0版本：初步完成了基于Numpy与Scipy模块进行WAV格式的音频声道混缩处理效果。

### Update on July 11th, 2022

v2.0, v3.0与v3.1版本：添加了增益，混响和变换采样等函数，丰富了程序的功能，并且重构了代码，改进了算法，在性能上做了一定的优化。

### Update on July 12th, 2022

v3.2版本：使用了Sonic链接库，可以实现快速的变调与变速处理。使用PyInstaller打包时可以将Sonic添加到根目录下。

### Update on July 18th, 2022

v3.3版本：规范了代码的书写格式。添加了命令行解析的功能，DAPC不仅可以作为Python模块被调用，也可以解析命令行参数，单独运行。

### Update on Sept. 14th, 2022

v3.4版本：添加了调用FFmpeg转码的功能。这样可以使任何主流格式的音频文件都能被解码为封装了PCM数据的WAV波形音频，供DAPC处理。使用PyInstaller打包时可以将FFmpeg添加到根目录下。

### Update on Jan. 28th, 2023

v3.5版本：添加了淡入淡出，3D环绕声和FFT滤波的功能。重构了代码，减少了音频数据在处理时因计算机字长而造成的差错。

### Update on March 14th, 2023

v3.5.1与v3.5.2版本：修复了程序中出现的部分Bug；为后续兼容即时编译(JIT)的Numba库与Cython相关用法做了“铺垫”。

## Statistics of Visiting Numbers  访问次数统计
<div>Number of Total Visits: &nbsp; <img src="https://visitor-badge.glitch.me/badge?page_id=Duyu09_NEW_Audio-Management_Core" /></div> 
