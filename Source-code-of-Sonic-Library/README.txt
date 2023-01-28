gcc编译生成动态库(sonic)的命令行：
gcc -o sonic.dll -shared -fPIC runsonic.c sonic.c sonic.h

不同OS对应的动态链接库拓展名不同：
Windows系统下为sonic.dll
MacOS系统下为sonic.dylib
Linux系统下为sonic.so
