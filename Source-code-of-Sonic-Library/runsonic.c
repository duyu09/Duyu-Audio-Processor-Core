#include "sonic.h"
void wavChangeSpeed(short *wav, short *result, int numChannels, int sampleRate, int wavLength, float speed){
    sonicStream stream;
    stream = sonicCreateStream(sampleRate, numChannels);
    sonicSetSpeed(stream, speed);
    int ret = sonicWriteShortToStream(stream, wav, wavLength);
    int numSamples = wavLength/speed;
    if (ret){
        int new_buffer_size = sonicReadShortFromStream(stream, result, numSamples);
    }
}