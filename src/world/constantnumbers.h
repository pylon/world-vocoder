//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise
// Author: mmorise [at] yamanashi.ac.jp (Masanori Morise)
// Last update: 2017/04/29
//
// This header file only defines constant numbers used for several function.
//-----------------------------------------------------------------------------
#ifndef WORLD_CONSTANT_NUMBERS_H_
#define WORLD_CONSTANT_NUMBERS_H_

namespace world {
  // for StoneMask()
  const float kFloorF0StoneMask = 40.0;

  const float kPi = 3.1415926535897932384;
  const float kMySafeGuardMinimum = 0.000000000001;
  const float kEps = 0.00000000000000022204460492503131;
  const float kFloorF0 = 71.0;
  const float kCeilF0 = 800.0;
  const float kDefaultF0 = 500.0;
  const float kLog2 = 0.69314718055994529;
  // Maximum standard deviation not to be selected as a best f0.
  const float kMaximumValue = 100000.0;
// Note to me (fs: 48000)
// 71 Hz is the limit to maintain the FFT size at 2048.
// If we use 70 Hz as FLOOR_F0, the FFT size of 4096 is required.

  // for D4C()
  const int kHanning = 1;
  const int kBlackman = 2;
  const float kFrequencyInterval = 3000.0;
  const float kUpperLimit = 15000.0;
  const float kThreshold = 0.85;
  const float kFloorF0D4C = 47.0;

  // for Codec (Mel scale)
  // S. Stevens & J. Volkmann,
  // The Relation of Pitch to Frequency: A Revised Scale,
  // American Journal of Psychology, vol. 53, no. 3, pp. 329-353, 1940.
  const float kM0 = 1127.01048;
  const float kF0 = 700.0;
  const float kFloorFrequency = 40.0;
  const float kCeilFrequency = 20000.0;

}  // namespace world

#endif  // WORLD_CONSTANT_NUMBERS_H_
