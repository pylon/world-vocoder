//-----------------------------------------------------------------------------
// Copyright 2017 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/05/09
//
// Coder/decoder functions for the spectral envelope and aperiodicity.
//-----------------------------------------------------------------------------
#include "world/codec.h"

#include <math.h>

#include "world/constantnumbers.h"
#include "world/fft.h"
#include "world/matlabfunctions.h"

namespace {
//-----------------------------------------------------------------------------
// Aperiodicity is initialized by the value 1.0 - world::kMySafeGuardMinimum.
// This value means the frame/frequency index is aperiodic.
//-----------------------------------------------------------------------------
static void InitializeAperiodicity(int f0_length, int fft_size,
    float **aperiodicity) {
  for (int i = 0; i < f0_length; ++i)
    for (int j = 0; j < fft_size / 2 + 1; ++j)
      aperiodicity[i][j] = 1.0 - world::kMySafeGuardMinimum;
}

//-----------------------------------------------------------------------------
// This function identifies whether this frame is voiced or unvoiced.
//-----------------------------------------------------------------------------
static int CheckVUV(const float *coarse_aperiodicity,
    int number_of_aperiodicities, float *tmp_aperiodicity) {
  float tmp = 0.0;
  for (int i = 0; i < number_of_aperiodicities; ++i) {
    tmp += coarse_aperiodicity[i];
    tmp_aperiodicity[i + 1] = coarse_aperiodicity[i];
  }
  tmp /= number_of_aperiodicities;

  return tmp > -0.5 ? 1 : 0;  // -0.5 is not optimized, but okay.
}

//-----------------------------------------------------------------------------
// Aperiodicity is obtained from the coded aperiodicity.
//-----------------------------------------------------------------------------
static void GetAperiodicity(const float *coarse_frequency_axis,
    const float *coarse_aperiodicity, int number_of_aperiodicities,
    const float *frequency_axis, int fft_size, float *aperiodicity) {
  interp1(coarse_frequency_axis, coarse_aperiodicity,
      number_of_aperiodicities + 2, frequency_axis, fft_size / 2 + 1,
      aperiodicity);
  for (int i = 0; i <= fft_size / 2; ++i)
    aperiodicity[i] = pow(10.0, aperiodicity[i] / 20.0);
}

//-----------------------------------------------------------------------------
// Frequency is converted into its mel representation.
//-----------------------------------------------------------------------------
static inline float FrequencyToMel(float frequency) {
  return world::kM0 * log(frequency / world::kF0 + 1.0);
}

//-----------------------------------------------------------------------------
// Mel is converted into frequency.
//-----------------------------------------------------------------------------
static inline float MelToFrequency(float mel) {
  return world::kF0 * (exp(mel / world::kM0) - 1.0);
}

//-----------------------------------------------------------------------------
// DCT for spectral envelope coding
//-----------------------------------------------------------------------------
static void DCTForCodec(const float *mel_spectrum, int max_dimension,
    const fft_complex *weight, const ForwardRealFFT *forward_real_fft,
    int number_of_dimensions, float *mel_cepstrum) {
  int bias = max_dimension / 2;
  for (int i = 0; i < max_dimension / 2; ++i) {
    forward_real_fft->waveform[i] = mel_spectrum[i * 2];
    forward_real_fft->waveform[i + bias] =
      mel_spectrum[max_dimension - (i * 2) - 1];
  }
  fft_execute(forward_real_fft->forward_fft);

  float normalization = sqrt(forward_real_fft->fft_size);
  for (int i = 0; i < number_of_dimensions; ++i)
    mel_cepstrum[i] = (forward_real_fft->spectrum[i][0] * weight[i][0] -
      forward_real_fft->spectrum[i][1] * weight[i][1]) / normalization;
}

//-----------------------------------------------------------------------------
// IDCT for spectral envelope decoding
//-----------------------------------------------------------------------------
static void IDCTForCodec(const float *mel_cepstrum, int max_dimension,
    const fft_complex *weight, const InverseComplexFFT *inverse_complex_fft,
    int number_of_dimensions, float *mel_spectrum) {
  float normalization = sqrt(inverse_complex_fft->fft_size);
  for (int i = 0; i < number_of_dimensions; ++i) {
    inverse_complex_fft->input[i][0] =
      mel_cepstrum[i] * weight[i][0] * normalization;
    inverse_complex_fft->input[i][1] =
      -mel_cepstrum[i] * weight[i][1] * normalization;
  }
  for (int i = number_of_dimensions; i < max_dimension; ++i) {
    inverse_complex_fft->input[i][0] = 0.0;
    inverse_complex_fft->input[i][1] = 0.0;
  }

  fft_execute(inverse_complex_fft->inverse_fft);

  for (int i = 0; i < max_dimension / 2; ++i) {
    mel_spectrum[i * 2] = inverse_complex_fft->output[i][0];
    mel_spectrum[(i * 2) + 1] =
      inverse_complex_fft->output[max_dimension - i - 1][0];
  }
}

//-----------------------------------------------------------------------------
// Spectral envelope in a frame is coded
//-----------------------------------------------------------------------------
static void CodeOneFrame(const float *log_spectral_envelope,
    const float *frequency_axis, int fft_size, const float *mel_axis,
    const fft_complex *weight, int max_dimension, int number_of_dimensions,
    const ForwardRealFFT *forward_real_fft, float *coded_spectral_envelope) {
  float *mel_spectrum = new float[max_dimension];
  interp1(frequency_axis, log_spectral_envelope, fft_size / 2 + 1,
      mel_axis, max_dimension, mel_spectrum);

  // DCT
  DCTForCodec(mel_spectrum, max_dimension, weight, forward_real_fft,
      number_of_dimensions, coded_spectral_envelope);

  delete[] mel_spectrum;
}

//-----------------------------------------------------------------------------
// Coded spectral envelope in a frame is decoded
//-----------------------------------------------------------------------------
static void DecodeOneFrame(const float *coded_spectral_envelope,
    const float *frequency_axis, int fft_size, const float *mel_axis,
    const fft_complex *weight, int max_dimension, int number_of_dimensions,
    const InverseComplexFFT *inverse_complex_fft, float *spectral_envelope) {
  float *mel_spectrum = new float[max_dimension + 2];

  // IDCT
  IDCTForCodec(coded_spectral_envelope, max_dimension, weight,
      inverse_complex_fft, number_of_dimensions, &mel_spectrum[1]);
  mel_spectrum[0] = mel_spectrum[1];
  mel_spectrum[max_dimension + 1] = mel_spectrum[max_dimension];

  interp1(mel_axis, mel_spectrum, max_dimension + 2, frequency_axis,
      fft_size / 2 + 1, spectral_envelope);

  for (int i = 0; i < fft_size / 2 + 1; ++i)
    spectral_envelope[i] = exp(spectral_envelope[i] / max_dimension);

  delete[] mel_spectrum;
}

//-----------------------------------------------------------------------------
// GetParameters() generates the required parameters.
//-----------------------------------------------------------------------------
static void GetParametersForCoding(float floor_frequency,
    float ceil_frequency, int fs, int fft_size, float *mel_axis,
    float *frequency_axis, fft_complex *weight) {
  int max_dimension = fft_size / 2;
  float floor_mel = FrequencyToMel(floor_frequency);
  float ceil_mel = FrequencyToMel(ceil_frequency);

  // Generate the mel axis and the weighting vector for DCT.
  for (int i = 0; i < max_dimension; ++i) {
    mel_axis[i] = (ceil_mel - floor_mel) * i / max_dimension + floor_mel;
    weight[i][0] = 2.0 * cos(i * world::kPi / fft_size) / sqrt(fft_size);
    weight[i][1] = 2.0 * sin(i * world::kPi / fft_size) / sqrt(fft_size);
  }
  weight[0][0] /= sqrt(2.0);

  // Generate the frequency axis on mel scale
  for (int i = 0; i < max_dimension; ++i)
    frequency_axis[i] = FrequencyToMel(static_cast<float>(i) * fs / fft_size);
}

//-----------------------------------------------------------------------------
// GetParameters() generates the required parameters.
//-----------------------------------------------------------------------------
static void GetParametersForDecoding(float floor_frequency,
    float ceil_frequency, int fs, int fft_size, int number_of_dimensions,
    float *mel_axis, float *frequency_axis, fft_complex *weight) {
  int max_dimension = fft_size / 2;
  float floor_mel = FrequencyToMel(floor_frequency);
  float ceil_mel = FrequencyToMel(ceil_frequency);

  // Generate the weighting vector for IDCT.
  for (int i = 0; i < number_of_dimensions; ++i) {
    weight[i][0] = cos(i * world::kPi / fft_size) * sqrt(fft_size);
    weight[i][1] = sin(i * world::kPi / fft_size) * sqrt(fft_size);
  }
  weight[0][0] /= sqrt(2.0);
  // Generate the mel axis for IDCT.
  for (int i = 0; i < max_dimension; ++i)
    mel_axis[i + 1] =
      MelToFrequency((ceil_mel - floor_mel) * i / max_dimension + floor_mel);
  mel_axis[0] = 0;
  mel_axis[max_dimension + 1] = fs / 2.0;

  // Generate the frequency axis
  for (int i = 0; i < fft_size / 2 + 1; ++i)
    frequency_axis[i] = static_cast<float>(i) * fs / fft_size;
}

}  // namespace

int GetNumberOfAperiodicities(int fs) {
  return  static_cast<int>(MyMinFloat(world::kUpperLimit, fs / 2.0 -
    world::kFrequencyInterval) / world::kFrequencyInterval);
}

void CodeAperiodicity(const float * const *aperiodicity, int f0_length,
    int fs, int fft_size, float **coded_aperiodicity) {
  int number_of_aperiodicities = GetNumberOfAperiodicities(fs);
  float *coarse_frequency_axis = new float[number_of_aperiodicities];
  for (int i = 0; i < number_of_aperiodicities; ++i)
    coarse_frequency_axis[i] = world::kFrequencyInterval * (i + 1.0);

  float *log_aperiodicity = new float[fft_size / 2 + 1];

  for (int i = 0; i < f0_length; ++i) {
    for (int j = 0; j < fft_size / 2 + 1; ++j)
      log_aperiodicity[j] = 20 * log10(aperiodicity[i][j]);
    interp1Q(0, static_cast<float>(fs) / fft_size, log_aperiodicity,
        fft_size / 2 + 1, coarse_frequency_axis, number_of_aperiodicities,
        coded_aperiodicity[i]);
  }

  delete[] coarse_frequency_axis;
  delete[] log_aperiodicity;
}

void DecodeAperiodicity(const float * const *coded_aperiodicity,
    int f0_length, int fs, int fft_size, float **aperiodicity) {
  InitializeAperiodicity(f0_length, fft_size, aperiodicity);
  int number_of_aperiodicities = GetNumberOfAperiodicities(fs);
  float *frequency_axis = new float[fft_size / 2 + 1];
  for (int i = 0; i <= fft_size / 2; ++i)
    frequency_axis[i] = static_cast<float>(fs) / fft_size * i;

  float *coarse_frequency_axis = new float[number_of_aperiodicities + 2];
  for (int i = 0; i <= number_of_aperiodicities; ++i)
    coarse_frequency_axis[i] = i * world::kFrequencyInterval;
  coarse_frequency_axis[number_of_aperiodicities + 1] = fs / 2.0;

  float *coarse_aperiodicity = new float[number_of_aperiodicities + 2];
  coarse_aperiodicity[0] = -60.0;
  coarse_aperiodicity[number_of_aperiodicities + 1] =
    -world::kMySafeGuardMinimum;

  for (int i = 0; i < f0_length; ++i) {
    if (CheckVUV(coded_aperiodicity[i], number_of_aperiodicities,
      coarse_aperiodicity) == 1) continue;
    GetAperiodicity(coarse_frequency_axis, coarse_aperiodicity,
        number_of_aperiodicities, frequency_axis, fft_size, aperiodicity[i]);
  }

  delete[] coarse_aperiodicity;
  delete[] coarse_frequency_axis;
  delete[] frequency_axis;
}

void CodeSpectralEnvelope(const float * const *spectrogram, int f0_length,
    int fs, int fft_size, int number_of_dimensions,
    float **coded_spectral_envelope) {
  float *mel_axis = new float[fft_size / 2];
  float *frequency_axis = new float[fft_size / 2 + 1];
  float *tmp_spectrum = new float[fft_size / 2 + 1];
  fft_complex *weight = new fft_complex[fft_size / 2];

  // Generation of the required parameters
  GetParametersForCoding(world::kFloorFrequency,
      MyMinFloat(fs / 2.0, world::kCeilFrequency), fs, fft_size,
      mel_axis, frequency_axis, weight);

  ForwardRealFFT forward_real_fft = { 0 };
  InitializeForwardRealFFT(fft_size / 2, &forward_real_fft);

  for (int i = 0; i < f0_length; ++i) {
    for (int j = 0; j < fft_size / 2 + 1; ++j)
      tmp_spectrum[j] = log(spectrogram[i][j]);
    CodeOneFrame(tmp_spectrum, frequency_axis, fft_size, mel_axis, weight,
        fft_size / 2, number_of_dimensions, &forward_real_fft,
        coded_spectral_envelope[i]);
  }

  DestroyForwardRealFFT(&forward_real_fft);
  delete[] weight;
  delete[] tmp_spectrum;
  delete[] frequency_axis;
  delete[] mel_axis;
}

void DecodeSpectralEnvelope(const float * const *coded_spectral_envelope,
    int f0_length, int fs, int fft_size, int number_of_dimensions,
    float **spectrogram) {
  float *mel_axis = new float[fft_size / 2 + 2];
  float *frequency_axis = new float[fft_size / 2 + 1];
  float *tmp_spectrum = new float[fft_size / 2 + 1];
  fft_complex *weight = new fft_complex[fft_size / 2];

  // Generation of the required parameters
  GetParametersForDecoding(world::kFloorFrequency,
      MyMinFloat(fs / 2.0, world::kCeilFrequency),
      fs, fft_size, number_of_dimensions, mel_axis, frequency_axis, weight);

  InverseComplexFFT inverse_complex_fft = { 0 };
  InitializeInverseComplexFFT(fft_size / 2, &inverse_complex_fft);

  for (int i = 0; i < f0_length; ++i) {
    DecodeOneFrame(coded_spectral_envelope[i], frequency_axis, fft_size,
        mel_axis, weight, fft_size / 2, number_of_dimensions,
        &inverse_complex_fft, spectrogram[i]);
  }

  DestroyInverseComplexFFT(&inverse_complex_fft);
  delete[] weight;
  delete[] tmp_spectrum;
  delete[] frequency_axis;
  delete[] mel_axis;
}
