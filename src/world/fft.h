//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/02/01
//
// These functions and variables are defined to use FFT as well as FFTW
// Please see fft.cpp to show the detailed information
//-----------------------------------------------------------------------------
#ifndef WORLD_FFT_H_
#define WORLD_FFT_H_

#include "world/macrodefinitions.h"

WORLD_BEGIN_C_DECLS

// Commands for FFT (This is the same as FFTW)
#define FFT_FORWARD 1
#define FFT_BACKWARD 2
#define FFT_ESTIMATE 3

// Complex number for FFT
typedef float fft_complex[2];
// Struct used for FFT
typedef struct {
  int n;
  int sign;
  unsigned int flags;
  fft_complex *c_in;
  float *in;
  fft_complex *c_out;
  float *out;
  float *input;
  int *ip;
  float *w;
} fft_plan;

fft_plan fft_plan_dft_1d(int n, fft_complex *in, fft_complex *out, int sign,
  unsigned int flags);
fft_plan fft_plan_dft_c2r_1d(int n, fft_complex *in, float *out,
  unsigned int flags);
fft_plan fft_plan_dft_r2c_1d(int n, float *in, fft_complex *out,
  unsigned int flags);
void fft_execute(fft_plan p);
void fft_destroy_plan(fft_plan p);

WORLD_END_C_DECLS

#endif  // WORLD_FFT_H_
