#include "RtOutput.h"
#include "WAV.h"
#include "STFT.h"
#include <string>
#include "MLDR.h"
#include "RtInput.h"

/* Set Parameter of Input */
constexpr int ch = 2;
constexpr int out_ch = 1;
constexpr int rate = 16000;
constexpr int frame = 512;
constexpr int shift = 128;

// Read output of AudioProbe() and set manually
#define DEVICE_OUTPUT 2

int main() {
  short* buf_1;
  short* raw = new short[ch * shift];



  AudioProbe();

  WAV input;
  RtInput rtinput(1/*device*/, ch, rate, shift, frame);
  //input.OpenFile("hello_world_ch1.wav");
  int samplerate_in = input.GetSampleRate();
  int samplerate_out = 48000;

  //RtOutput speaker(DEVICE_OUTPUT,1,samplerate_in,samplerate_out,128,512);




  WAV output(ch, rate);
  rtinput.Start();
  output.NewFile("output.wav");
  /*
  input.Rewind();
  buf_1 = new short[input.GetSize()];
  fread(buf_1, sizeof(short), input.GetSize(), input.GetFilePointer());

  speaker.FullBufLoad(buf_1, input.GetSize());
  speaker.Start();
  speaker.Wait();
  */
  int cnt = 0;
  while (cnt<500) {
      if (rtinput.data.stock.load() >= shift * ch) {
          rtinput.GetBuffer(raw);
          output.Append(raw, ch * shift);
          cnt++;
          printf("asd");
      }
      else {
          //printf("else");
          continue;
      }

  }
  output.Finish();
  return 0;
}

// Modification of AudioProbe()
int GetDefaultOutputDevice() {
  RtAudio audio;
  unsigned int devices = audio.getDeviceCount();
  RtAudio::DeviceInfo info;
  for (unsigned int i = 0; i < devices; i++) {
    info = audio.getDeviceInfo(i);
    if (info.probed == true) {
      //if (info.inputChannels != 0) {
      if (true) {
        std::cout << "device = " << i << "\n";
        std::cout << "name = " << info.name << "\n";
        std::cout << "maximum input channels = " << info.inputChannels << "\n";
        std::cout << "maximum output channels = " << info.outputChannels << "\n";
        std::cout << "Samplerates : ";
        for (auto sr : info.sampleRates)
          std::cout << sr << " ";
        std::cout << "\n";
        std::cout << "----------------------------------------------------------" << "\n";
      }
    }
  }
  return 0;
}
/*
int main() {
    _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
    // _crtBreakAlloc = 11917;

    MLDR mldr(frame, ch);

    int length;
    WAV input;
    WAV output(out_ch, rate);
    STFT stft(ch, frame, shift);

    short* buf_in;
    double** data;
    short* buf_out;

    data = new double* [ch];
    for (int i = 0; i < ch; i++) {
        data[i] = new double[frame + 2];
        memset(data[i], 0, sizeof(double) * (frame + 2));
    }

    buf_in = new short[ch * shift];
    memset(buf_in, 0, sizeof(short) * ch * shift);
    buf_out = new short[out_ch * shift];

    input.OpenFile("F01_22GC010A_BUS.wav");
    output.NewFile("output.wav");

    int cnt = 0;
    while (!input.IsEOF()) {
        cnt++;
        length = input.ReadUnit(buf_in, shift * ch);
        stft.stft(buf_in, length, data);
        mldr.Process(data);

        stft.istft(data[0], buf_out);
        output.Append(buf_out, shift * out_ch);
    }
    output.Finish();
    input.Finish();

    for (int i = 0; i < ch; i++)
        delete[] data[i];
    delete[] data;
    delete[] buf_in;
    delete[] buf_out;

    return 0;
}
*/