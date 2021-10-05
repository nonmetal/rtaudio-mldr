# STFT

STFT(Short Time Fourier Transform), ISTFT(Inverse - Short Time Fourier Transform) for wav,mic input  

provides 25%,50% overlap STFT process.  

## NOTE

```git clone --recursive https://github.com/kooBH/STFT.git```

To build test code you need to clone ```--recursive``` to use submodule  ```koobh/WAV.git```
If you already cloned, then,
```
git submoudle init
git submodule update
```
to use submodule.

## About FFT  
I'm currently using FFT of [Ooura](http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html).  
Since, it is fastest FFT in a single header file.   
But, sometimes (usually  not), there are errors between MATLAB FFT output and Ooura FFT output.     
If you need to perfectly same output as MATLAB, you have to use other FFT library.  

## EXAMPLE

+ Directly from buffer  

```cpp

#include "cpp/STFT.h"

...

// frame/shift must be 4 or 2
STFT process(channels,frame_size,shift_size);
WAV input;
WAV output(channels,sample_rate);
input.OpenFile('input.wav');

short* buf_in; //[channels * shift_size]
double **data; //[channels][frame_size]
short* buf_out; //[channels * shift_size]

while(!input.IsEOF()){
    length = input.ReadUnit(buf_in,shift*channels);
    process.stft(buf_in,length,data);
    process.istft(data,buf_out);
    output.Append(buf_out,shift*channels);
  }

```

+ From 2D array

```cpp

#include "cpp/STFT.h"

...

// frame/shift must be 4 or 2
STFT process(channels,frame_size,shift_size);
WAV input;
WAV output(channels,sample_rate);
input.OpenFile('input.wav');

double **raw; //[channels][shift_size]
double **data; //[channels][frame_size]
short* buf_out; //[channels * shift_size]

while(!input.IsEOF()){
    input.Convert2Array(raw);
    process.stft(raw,data);
    process.istft(data,buf_out);
    output.Append(buf_out,shift*channels);
  }

```

## STATUS

C       : only STFT is implemented, iSTFT is not implemented yet   
C++    : verified (same output as MATLAB routine)  
