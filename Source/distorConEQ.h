//
// distorConEQ.h
//
// Code generation for function 'distorConEQ'
//

#ifndef DISTORCONEQ_H
#define DISTORCONEQ_H

// Include files
#include "rtwtypes.h"
#include <cstddef>
#include <cstdlib>

// Type Declarations
struct distorConEQStackData;

// Type Definitions
struct dsp_BiquadFilter_1 {
  int S0_isInitialized;
  double W0_FILT_STATES[4];
  int W1_PreviousNumChannels;
  double P0_ICRTP;
};

struct cell_wrap_2 {
  unsigned int f1[8];
};

namespace coder {
namespace dspcodegen {
class BiquadFilter {
public:
  BiquadFilter();
  ~BiquadFilter();
  boolean_T matlabCodegenIsDeleted;
  int isInitialized;
  boolean_T isSetupComplete;
  dsp_BiquadFilter_1 cSFunObject;
};

} // namespace dspcodegen
class multibandParametricEQ {
public:
  void designFilters();
  multibandParametricEQ();
  ~multibandParametricEQ();
  boolean_T matlabCodegenIsDeleted;
  int isInitialized;
  boolean_T isSetupComplete;
  boolean_T TunablePropsChanged;
  cell_wrap_2 inputVarSize[1];
  double pSampleRateDialog;
  dspcodegen::BiquadFilter SOSFilterObj;
  double NumMatrix[3];
  double DenMatrix[2];
  double NumChannels;
  double privFrequencies[10];
  double privQualityFactors[10];
  double privPeakGains[10];
  boolean_T AreFiltersDesigned;
};

} // namespace coder
class derivedAudioPlugin {
public:
  derivedAudioPlugin();
  ~derivedAudioPlugin();
  boolean_T matlabCodegenIsDeleted;
  double PrivateSampleRate;
  int PrivateLatency;
  double gain;
  double level;
  coder::multibandParametricEQ EQ;
};

// Function Declarations
extern void createPluginInstance(distorConEQStackData *SD,
                                 unsigned long long thisPtr);

extern void distorConEQ_initialize(distorConEQStackData *SD);

extern void distorConEQ_terminate(distorConEQStackData *SD);

extern int getLatencyInSamplesCImpl(distorConEQStackData *SD);

extern void onParamChangeCImpl(distorConEQStackData *SD, int paramIdx,
                               double value);

extern void processEntryPoint(distorConEQStackData *SD, double samplesPerFrame,
                              const double i1_data[], const int i1_size[1],
                              const double i2_data[], const int i2_size[1],
                              double o1_data[], int o1_size[1],
                              double o2_data[], int o2_size[1]);

extern void resetCImpl(distorConEQStackData *SD, double rate);

#endif
// End of code generation (distorConEQ.h)
