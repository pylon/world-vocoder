//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/04/01
//
// Summary:
// This example estimates the spectral envelope from an audio file
// and then saves the result to a file.
//
// How to use:
// % spanalysis -h
//
// Related works: f0analysis.cpp, apanalysis.cpp, readandsynthesis.cpp
//-----------------------------------------------------------------------------
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "../../tools/audioio.h"
#include "../../tools/parameterio.h"
#include "world/cheaptrick.h"
#include "world/codec.h"
#include "world/constantnumbers.h"

namespace {

//-----------------------------------------------------------------------------
// Display how to use this program
//-----------------------------------------------------------------------------
void usage(char *argv) {
  printf("\n");
  printf(" %s - spectral envelope estimation by CheapTrick\n", argv);
  printf("\n");
  printf("  usage:\n");
  printf("   %s input.wav input.f0 [options]\n", argv);
  printf("  options:\n");
  printf("   -f f    : FFT size (samples)            [variable]\n");
  printf("           : Default depends on fs (44100 -> 2048, 16000 -> 1024)\n");
  printf("   -q q    : compensation coefficient      [-0.15]\n");
  printf("           : I don't recommend to change this value.\n");
  printf("   -d d    : number of coefficients        [0 (without coding)]\n");
  printf("           : Spectral envelope is decoded by these coefficients.\n");
  printf("           : You must not set this value above the half of\n");
  printf("           : the FFT size.\n");
  printf("   -o name : filename used for output      [output.sp]\n");
  printf("\n");
}

//-----------------------------------------------------------------------------
// Set parameters from command line options
//-----------------------------------------------------------------------------
int SetOption(int argc, char **argv, int *fft_size, float *q1,
    int *number_of_dimensions, char *filename) {
  while (--argc) {
    if (strcmp(argv[argc], "-f") == 0) *fft_size = atoi(argv[argc + 1]);
    if (strcmp(argv[argc], "-q") == 0) *q1 = atof(argv[argc + 1]);
    if (strcmp(argv[argc], "-d") == 0)
      *number_of_dimensions = atof(argv[argc + 1]);
    if (strcmp(argv[argc], "-o") == 0)
      snprintf(filename, 200, argv[argc + 1]);
    if (strcmp(argv[argc], "-h") == 0) {
      usage(argv[0]);
      return 0;
    }
  }
  return 1;
}

void WriteCodedSpectralEnvelope(char *filename,
    const float * const *spectrogram, int fs, int f0_length,
    float frame_period, int number_of_dimensions, int fft_size) {
  float **coded_spectrogram = new float *[f0_length];
  for (int i = 0; i < f0_length; ++i)
    coded_spectrogram[i] = new float[number_of_dimensions];

  CodeSpectralEnvelope(spectrogram, f0_length, fs, fft_size,
      number_of_dimensions, coded_spectrogram);

  WriteSpectralEnvelope(filename, fs, f0_length, frame_period, fft_size,
      number_of_dimensions, coded_spectrogram);
  for (int i = 0; i < f0_length; ++i) delete[] coded_spectrogram[i];
  delete[] coded_spectrogram;
}

}  // namespace

//-----------------------------------------------------------------------------
// This example estimates the spectral envelope from an audio file
// and then saves the result to a file.
//-----------------------------------------------------------------------------
int main(int argc, char **argv) {
  // Command check
  if (argc < 2) return 0;
  if (0 == strcmp(argv[1], "-h")) {
    usage(argv[0]);
    return -1;
  }

  // Read F0 information
  int f0_length = static_cast<int>(GetHeaderInformation(argv[2], "NOF "));
  float frame_period = GetHeaderInformation(argv[2], "FP  ");
  float *f0 = new float[f0_length];
  float *temporal_positions = new float[f0_length];
  ReadF0(argv[2], temporal_positions, f0);

  // Read an audio file
  int x_length = GetAudioLength(argv[1]);
  if (x_length <= 0) {
    if (x_length == 0) {
      printf("error: File not found.\n");
    } else {
      printf("error: File is not .wav format.\n");
    }
    return -1;
  }
  float *x = new float[x_length];
  int fs, nbit;
  wavread(argv[1], &fs, &nbit, x);

  // Default parameters
  CheapTrickOption option = { 0 };
  InitializeCheapTrickOption(fs, &option);
  char filename[200] = "output.sp";
  int number_of_dimensions = 0;

  // Options from command line
  if (SetOption(argc, argv, &option.fft_size, &option.q1,
    &number_of_dimensions, filename) == 0) return -1;

  // Spectral envelope analysis
  float **spectrogram = new float *[f0_length];
  for (int i = 0; i < f0_length; ++i)
    spectrogram[i] = new float[option.fft_size / 2 + 1];
  CheapTrick(x, x_length, fs, temporal_positions, f0, f0_length, &option,
      spectrogram);

  // File output
  if (number_of_dimensions == 0) {
    // If you want to write the raw spectral envelope.
    WriteSpectralEnvelope(filename, fs, f0_length, frame_period,
        option.fft_size, 0, spectrogram);
  } else {
    // If you want to write the coded spectral envelope.
    WriteCodedSpectralEnvelope(filename, spectrogram, fs, f0_length,
        frame_period, number_of_dimensions, option.fft_size);
  }

  // Memory deallocation
  for (int i = 0; i < f0_length; ++i) delete[] spectrogram[i];
  delete[] spectrogram;
  delete[] f0;
  delete[] temporal_positions;
  delete[] x;
  return 0;
}
