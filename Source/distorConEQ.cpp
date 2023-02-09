//
// distorConEQ.cpp
//
// Code generation for function 'distorConEQ'
//

// Include files
#include "distorConEQ.h"
#include "distorConEQ_types.h"
#include "rt_nonfinite.h"
#include "coder_array.h"
#include <cmath>
#include <cstring>

// Function Declarations
namespace coder {
namespace audio {
namespace internal {
static void designHPEQFilter(double G, double GB, double w0, double BW,
                             double B[5], double A[5]);

}
} // namespace audio
} // namespace coder
static derivedAudioPlugin *getPluginInstance(distorConEQStackData *SD);

static void getPluginInstance_free(distorConEQStackData *SD);

static void getPluginInstance_init(distorConEQStackData *SD);

static double rt_powd_snf(double u0, double u1);

// Function Definitions
namespace coder {
void multibandParametricEQ::designFilters()
{
  double minval[10];
  double varargin_2[10];
  double Af[5];
  double Bf[5];
  double GB;
  double unnamed_idx_0;
  double unnamed_idx_0_tmp;
  int k;
  boolean_T x[2];
  GB = this->pSampleRateDialog / 2.0;
  for (k = 0; k < 10; k++) {
    varargin_2[k] = this->privPeakGains[k];
    minval[k] = std::fmin(GB, this->privFrequencies[k]) / GB;
  }
  unnamed_idx_0 = minval[0];
  unnamed_idx_0_tmp = minval[0] / this->privQualityFactors[0];
  if (minval[0] > 1.0) {
    unnamed_idx_0 = 1.0;
  }
  if (unnamed_idx_0_tmp > 1.0) {
    unnamed_idx_0_tmp = 1.0;
  }
  if (unnamed_idx_0_tmp < 0.0) {
    unnamed_idx_0_tmp = 0.0;
  }
  GB = varargin_2[0] / 2.0;
  if (std::isinf(varargin_2[0]) && (varargin_2[0] < 0.0)) {
    GB = -3.0102999566398121;
  }
  this->NumMatrix[0] = 0.0;
  this->NumMatrix[1] = 0.0;
  this->NumMatrix[2] = 0.0;
  this->NumMatrix[0] = 1.0;
  this->DenMatrix[0] = 0.0;
  this->DenMatrix[1] = 0.0;
  if (!(std::abs(varargin_2[0]) <= 2.2204460492503131E-16)) {
    double Gsq;
    Gsq = rt_powd_snf(10.0, varargin_2[0] / 10.0);
    GB = rt_powd_snf(10.0, GB / 10.0);
    if ((!(std::abs(Gsq - GB) <= 2.2204460492503131E-16)) &&
        (!(std::abs(GB - 1.0) <= 2.2204460492503131E-16))) {
      boolean_T exitg1;
      boolean_T guard1{false};
      boolean_T y;
      audio::internal::designHPEQFilter(Gsq, GB, unnamed_idx_0,
                                        unnamed_idx_0_tmp, Bf, Af);
      x[0] = true;
      x[1] = true;
      if (!(Bf[3] == 0.0)) {
        x[0] = false;
      }
      if (!(Bf[4] == 0.0)) {
        x[1] = false;
      }
      y = true;
      k = 0;
      exitg1 = false;
      while ((!exitg1) && (k < 2)) {
        if (!x[k]) {
          y = false;
          exitg1 = true;
        } else {
          k++;
        }
      }
      guard1 = false;
      if (y) {
        x[0] = true;
        x[1] = true;
        if (!(Af[3] == 0.0)) {
          x[0] = false;
        }
        if (!(Af[4] == 0.0)) {
          x[1] = false;
        }
        y = true;
        k = 0;
        exitg1 = false;
        while ((!exitg1) && (k < 2)) {
          if (!x[k]) {
            y = false;
            exitg1 = true;
          } else {
            k++;
          }
        }
        if (y) {
          this->NumMatrix[0] = Bf[0];
          this->NumMatrix[1] = Bf[1];
          this->NumMatrix[2] = Bf[2];
          this->DenMatrix[0] = Af[1];
          this->DenMatrix[1] = Af[2];
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }
      if (guard1) {
        this->NumMatrix[0] = Bf[0];
        this->NumMatrix[1] = Bf[1];
        this->NumMatrix[2] = Bf[2];
        this->DenMatrix[0] = Af[1];
        this->DenMatrix[1] = Af[2];
      }
    }
  }
  this->AreFiltersDesigned = true;
}

namespace audio {
namespace internal {
static void designHPEQFilter(double G, double GB, double w0, double BW,
                             double B[5], double A[5])
{
  double Ahat_idx_1;
  double Ba_idx_0;
  double Bhat_idx_0;
  double Bhat_idx_1;
  double WB;
  double a;
  double c0;
  int Ahat_idx_0;
  signed char ii_size_idx_0;
  signed char ii_size_idx_1;
  WB = std::tan(3.1415926535897931 * BW / 2.0);
  a = std::sqrt((G - GB) / (GB - 1.0));
  Ba_idx_0 = std::sqrt(G) * WB;
  for (Ahat_idx_0 = 0; Ahat_idx_0 < 5; Ahat_idx_0++) {
    B[Ahat_idx_0] = 0.0;
    A[Ahat_idx_0] = 0.0;
  }
  Bhat_idx_0 = 0.0;
  Ahat_idx_0 = 0;
  Bhat_idx_1 = 0.0;
  Ahat_idx_1 = 0.0;
  if (w0 == 0.0) {
    c0 = 1.0;
  } else if (w0 == 1.0) {
    c0 = -1.0;
  } else if (w0 == 0.5) {
    c0 = 0.0;
  } else {
    c0 = std::cos(3.1415926535897931 * w0);
  }
  if (a != 0.0) {
    ii_size_idx_0 = 1;
    ii_size_idx_1 = 1;
    Ahat_idx_1 = WB + a;
    Bhat_idx_0 = (Ba_idx_0 + a) / Ahat_idx_1;
    Bhat_idx_1 = (Ba_idx_0 - a) / Ahat_idx_1;
    Ahat_idx_0 = 1;
    Ahat_idx_1 = (WB - a) / Ahat_idx_1;
  } else {
    ii_size_idx_0 = 0;
    ii_size_idx_1 = 0;
  }
  if ((c0 == 1.0) || (c0 == -1.0)) {
    B[0] = Bhat_idx_0;
    A[0] = Ahat_idx_0;
    B[2] = 0.0;
    A[2] = 0.0;
    B[1] = c0 * Bhat_idx_1;
    A[1] = c0 * Ahat_idx_1;
  } else if ((ii_size_idx_0 != 0) && (ii_size_idx_1 != 0)) {
    B[0] = Bhat_idx_0;
    B[1] = c0 * (Bhat_idx_1 - Bhat_idx_0);
    B[2] = -Bhat_idx_1;
    A[0] = 1.0;
    A[1] = c0 * (Ahat_idx_1 - 1.0);
    A[2] = -Ahat_idx_1;
  }
}

} // namespace internal
} // namespace audio
} // namespace coder
static derivedAudioPlugin *getPluginInstance(distorConEQStackData *SD)
{
  static const short iv[10]{100,  181,  325,  585,   1053,
                            1896, 3415, 6151, 11078, 19953};
  coder::multibandParametricEQ *obj;
  if (!SD->pd->plugin_not_empty) {
    int i;
    //  Pass constructor args to plugin.
    SD->pd->plugin.gain = 1.0;
    SD->pd->plugin.level = 1.0;
    SD->pd->plugin.PrivateLatency = 0;
    if (!SD->pd->thisPtr_not_empty) {
      SD->pd->thisPtr = 0ULL;
      SD->pd->thisPtr_not_empty = true;
    }
    //  Creación del objeto EQ
    obj = &SD->pd->plugin.EQ;
    SD->pd->plugin.EQ.pSampleRateDialog = 44100.0;
    for (i = 0; i < 10; i++) {
      SD->pd->plugin.EQ.privFrequencies[i] = iv[i];
    }
    for (i = 0; i < 10; i++) {
      SD->pd->plugin.EQ.privQualityFactors[i] = 1.6;
    }
    std::memset(&SD->pd->plugin.EQ.privPeakGains[0], 0, 10U * sizeof(double));
    SD->pd->plugin.EQ.isInitialized = 0;
    obj->SOSFilterObj.isInitialized = 0;
    // System object Constructor function: dsp.BiquadFilter
    SD->pd->plugin.EQ.SOSFilterObj.cSFunObject.P0_ICRTP = 0.0;
    obj->SOSFilterObj.matlabCodegenIsDeleted = false;
    SD->pd->plugin.EQ.AreFiltersDesigned = false;
    SD->pd->plugin.EQ.NumChannels = -1.0;
    SD->pd->plugin.EQ.AreFiltersDesigned = false;
    SD->pd->plugin.EQ.matlabCodegenIsDeleted = false;
    SD->pd->plugin.matlabCodegenIsDeleted = false;
    SD->pd->plugin_not_empty = true;
  }
  return &SD->pd->plugin;
}

static void getPluginInstance_free(distorConEQStackData *SD)
{
  coder::dspcodegen::BiquadFilter *b_obj;
  coder::multibandParametricEQ *obj;
  if (!SD->pd->plugin.matlabCodegenIsDeleted) {
    SD->pd->plugin.matlabCodegenIsDeleted = true;
  }
  obj = &SD->pd->plugin.EQ;
  if (!obj->matlabCodegenIsDeleted) {
    obj->matlabCodegenIsDeleted = true;
    if (obj->isInitialized == 1) {
      obj->isInitialized = 2;
      if (obj->isSetupComplete) {
        if (obj->SOSFilterObj.isInitialized == 1) {
          obj->SOSFilterObj.isInitialized = 2;
        }
        obj->NumChannels = -1.0;
      }
    }
  }
  b_obj = &SD->pd->plugin.EQ.SOSFilterObj;
  if (!b_obj->matlabCodegenIsDeleted) {
    b_obj->matlabCodegenIsDeleted = true;
    if (b_obj->isInitialized == 1) {
      b_obj->isInitialized = 2;
    }
  }
}

static void getPluginInstance_init(distorConEQStackData *SD)
{
  SD->pd->plugin_not_empty = false;
  SD->pd->plugin.EQ.SOSFilterObj.matlabCodegenIsDeleted = true;
  SD->pd->plugin.EQ.matlabCodegenIsDeleted = true;
  SD->pd->plugin.matlabCodegenIsDeleted = true;
}

static double rt_powd_snf(double u0, double u1)
{
  double y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = rtNaN;
  } else {
    double d;
    double d1;
    d = std::abs(u0);
    d1 = std::abs(u1);
    if (std::isinf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = std::pow(u0, u1);
    }
  }
  return y;
}

namespace coder {
namespace dspcodegen {
BiquadFilter::BiquadFilter()
{
}

BiquadFilter::~BiquadFilter()
{
}

} // namespace dspcodegen
multibandParametricEQ::multibandParametricEQ()
{
}

multibandParametricEQ::~multibandParametricEQ()
{
}

} // namespace coder
derivedAudioPlugin::derivedAudioPlugin()
{
}

derivedAudioPlugin::~derivedAudioPlugin()
{
}

void createPluginInstance(distorConEQStackData *SD, unsigned long long thisPtr)
{
  if (!SD->pd->thisPtr_not_empty) {
    SD->pd->thisPtr = thisPtr;
    SD->pd->thisPtr_not_empty = true;
  }
  getPluginInstance(SD);
}

void distorConEQ_initialize(distorConEQStackData *SD)
{
  SD->pd->thisPtr_not_empty = false;
  getPluginInstance_init(SD);
}

void distorConEQ_terminate(distorConEQStackData *SD)
{
  getPluginInstance_free(SD);
}

int getLatencyInSamplesCImpl(distorConEQStackData *SD)
{
  derivedAudioPlugin *plugin;
  plugin = getPluginInstance(SD);
  return plugin->PrivateLatency;
}

void onParamChangeCImpl(distorConEQStackData *SD, int paramIdx, double value)
{
  derivedAudioPlugin *plugin;
  double dv[10];
  plugin = getPluginInstance(SD);
  switch (paramIdx) {
  case 0:
    plugin->gain = value;
    break;
  case 1: {
    double x1;
    boolean_T flag;
    //  Al detectar un cambio en
    //  el fader, actualiza la
    flag = (plugin->EQ.isInitialized == 1);
    if (flag) {
      plugin->EQ.TunablePropsChanged = true;
    }
    x1 = plugin->EQ.privFrequencies[0];
    flag = (x1 == value);
    if (!flag) {
      int i;
      for (i = 0; i < 10; i++) {
        dv[i] = plugin->EQ.privFrequencies[i];
      }
      dv[0] = value;
      for (i = 0; i < 10; i++) {
        plugin->EQ.privFrequencies[i] = dv[i];
      }
      plugin->EQ.AreFiltersDesigned = false;
    }
    //  frecuencia en el EQ
  } break;
  case 2: {
    double x1;
    boolean_T flag;
    //  Al detectar un cambio en
    //  el fader, actualiza el
    flag = (plugin->EQ.isInitialized == 1);
    if (flag) {
      plugin->EQ.TunablePropsChanged = true;
    }
    x1 = plugin->EQ.privPeakGains[0];
    flag = (x1 == value);
    if (!flag) {
      int i;
      for (i = 0; i < 10; i++) {
        dv[i] = plugin->EQ.privPeakGains[i];
      }
      dv[0] = value;
      for (i = 0; i < 10; i++) {
        plugin->EQ.privPeakGains[i] = dv[i];
      }
      plugin->EQ.AreFiltersDesigned = false;
    }
    //  pico de amplitud en el EQ
  } break;
  case 3:
    plugin->level = value;
    break;
  }
}

void processEntryPoint(distorConEQStackData *SD, double samplesPerFrame,
                       const double i1_data[], const int i1_size[1],
                       const double i2_data[], const int i2_size[1],
                       double o1_data[], int o1_size[1], double o2_data[],
                       int o2_size[1])
{
  coder::multibandParametricEQ *obj;
  derivedAudioPlugin *plugin;
  coder::array<double, 2U> out;
  coder::array<double, 2U> w;
  cell_wrap_2 varSizes;
  dsp_BiquadFilter_1 *b_obj;
  double A_idx_0;
  double A_idx_1;
  double B_idx_0;
  double B_idx_1;
  double B_idx_2;
  double numAccum;
  double stageOut;
  int k;
  int nx;
  short inSize[8];
  boolean_T exitg1;
  plugin = getPluginInstance(SD);
  B_idx_0 = plugin->gain;
  out.set_size(i1_size[0], 2);
  nx = i1_size[0];
  for (k = 0; k < nx; k++) {
    out[k] = i1_data[k] * B_idx_0;
  }
  nx = i2_size[0];
  for (k = 0; k < nx; k++) {
    out[k + out.size(0)] = i2_data[k] * B_idx_0;
  }
  nx = out.size(0) << 1;
  for (k = 0; k < nx; k++) {
    out[k] = std::atan(out[k]);
  }
  nx = out.size(0) * 2;
  out.set_size(out.size(0), 2);
  for (k = 0; k < nx; k++) {
    out[k] = 0.63661977236758138 * out[k];
  }
  //  1. Soft clipping
  obj = &plugin->EQ;
  if (plugin->EQ.isInitialized != 1) {
    plugin->EQ.isSetupComplete = false;
    plugin->EQ.isInitialized = 1;
    varSizes.f1[0] = static_cast<unsigned int>(out.size(0));
    varSizes.f1[1] = 2U;
    for (k = 0; k < 6; k++) {
      varSizes.f1[k + 2] = 1U;
    }
    plugin->EQ.inputVarSize[0] = varSizes;
    plugin->EQ.NumChannels = 2.0;
    plugin->EQ.designFilters();
    plugin->EQ.isSetupComplete = true;
    plugin->EQ.TunablePropsChanged = false;
    if (obj->SOSFilterObj.isInitialized == 1) {
      // System object Initialization function: dsp.BiquadFilter
      plugin->EQ.SOSFilterObj.cSFunObject.W0_FILT_STATES[0] =
          plugin->EQ.SOSFilterObj.cSFunObject.P0_ICRTP;
      plugin->EQ.SOSFilterObj.cSFunObject.W0_FILT_STATES[1] =
          plugin->EQ.SOSFilterObj.cSFunObject.P0_ICRTP;
      plugin->EQ.SOSFilterObj.cSFunObject.W0_FILT_STATES[2] =
          plugin->EQ.SOSFilterObj.cSFunObject.P0_ICRTP;
      plugin->EQ.SOSFilterObj.cSFunObject.W0_FILT_STATES[3] =
          plugin->EQ.SOSFilterObj.cSFunObject.P0_ICRTP;
      plugin->EQ.SOSFilterObj.cSFunObject.W1_PreviousNumChannels = -1;
    }
  }
  if (plugin->EQ.TunablePropsChanged) {
    plugin->EQ.TunablePropsChanged = false;
    if (!plugin->EQ.AreFiltersDesigned) {
      plugin->EQ.designFilters();
    }
  }
  inSize[0] = static_cast<short>(out.size(0));
  inSize[1] = 2;
  for (k = 0; k < 6; k++) {
    inSize[k + 2] = 1;
  }
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 8)) {
    if (obj->inputVarSize[0].f1[k] != static_cast<unsigned int>(inSize[k])) {
      for (k = 0; k < 8; k++) {
        obj->inputVarSize[0].f1[k] = static_cast<unsigned int>(inSize[k]);
      }
      exitg1 = true;
    } else {
      k++;
    }
  }
  B_idx_0 = plugin->EQ.NumMatrix[0];
  B_idx_1 = plugin->EQ.NumMatrix[1];
  B_idx_2 = plugin->EQ.NumMatrix[2];
  A_idx_0 = plugin->EQ.DenMatrix[0];
  A_idx_1 = plugin->EQ.DenMatrix[1];
  if (obj->SOSFilterObj.isInitialized != 1) {
    obj->SOSFilterObj.isSetupComplete = false;
    obj->SOSFilterObj.isInitialized = 1;
    obj->SOSFilterObj.isSetupComplete = true;
    // System object Initialization function: dsp.BiquadFilter
    plugin->EQ.SOSFilterObj.cSFunObject.W0_FILT_STATES[0] =
        plugin->EQ.SOSFilterObj.cSFunObject.P0_ICRTP;
    plugin->EQ.SOSFilterObj.cSFunObject.W0_FILT_STATES[1] =
        plugin->EQ.SOSFilterObj.cSFunObject.P0_ICRTP;
    plugin->EQ.SOSFilterObj.cSFunObject.W0_FILT_STATES[2] =
        plugin->EQ.SOSFilterObj.cSFunObject.P0_ICRTP;
    plugin->EQ.SOSFilterObj.cSFunObject.W0_FILT_STATES[3] =
        plugin->EQ.SOSFilterObj.cSFunObject.P0_ICRTP;
    plugin->EQ.SOSFilterObj.cSFunObject.W1_PreviousNumChannels = -1;
  }
  b_obj = &obj->SOSFilterObj.cSFunObject;
  // System object Outputs function: dsp.BiquadFilter
  w.set_size(out.size(0), 2);
  nx = 0;
  if (plugin->EQ.SOSFilterObj.cSFunObject.W1_PreviousNumChannels == -1) {
    plugin->EQ.SOSFilterObj.cSFunObject.W1_PreviousNumChannels = 2;
  }
  for (k = 0; k < out.size(0); k++) {
    numAccum = b_obj->W0_FILT_STATES[0];
    stageOut = numAccum + B_idx_0 * out[nx];
    numAccum = b_obj->W0_FILT_STATES[1];
    b_obj->W0_FILT_STATES[0] =
        (numAccum + B_idx_1 * out[nx]) - A_idx_0 * stageOut;
    b_obj->W0_FILT_STATES[1] = B_idx_2 * out[nx] - A_idx_1 * stageOut;
    w[nx] = stageOut;
    nx++;
  }
  for (k = 0; k < out.size(0); k++) {
    numAccum = b_obj->W0_FILT_STATES[2];
    stageOut = numAccum + B_idx_0 * out[nx];
    numAccum = b_obj->W0_FILT_STATES[3];
    b_obj->W0_FILT_STATES[2] =
        (numAccum + B_idx_1 * out[nx]) - A_idx_0 * stageOut;
    b_obj->W0_FILT_STATES[3] = B_idx_2 * out[nx] - A_idx_1 * stageOut;
    w[nx] = stageOut;
    nx++;
  }
  //  2. Aplica EQ
  B_idx_0 = plugin->level;
  nx = w.size(0) * 2;
  for (k = 0; k < nx; k++) {
    w[k] = w[k] * B_idx_0;
  }
  //  3. Nivel de salida
  if (1.0 > samplesPerFrame) {
    nx = 0;
  } else {
    nx = static_cast<int>(samplesPerFrame);
  }
  o1_size[0] = nx;
  for (k = 0; k < nx; k++) {
    o1_data[k] = w[k];
  }
  if (1.0 > samplesPerFrame) {
    nx = 0;
  } else {
    nx = static_cast<int>(samplesPerFrame);
  }
  o2_size[0] = nx;
  for (k = 0; k < nx; k++) {
    o2_data[k] = w[k + w.size(0)];
  }
}

void resetCImpl(distorConEQStackData *SD, double rate)
{
  derivedAudioPlugin *plugin;
  double value;
  boolean_T flag;
  plugin = getPluginInstance(SD);
  plugin->PrivateSampleRate = rate;
  flag = (plugin->EQ.isInitialized == 1);
  if (flag) {
    plugin->EQ.TunablePropsChanged = true;
  }
  value = plugin->PrivateSampleRate;
  plugin->EQ.pSampleRateDialog = value;
  plugin->EQ.AreFiltersDesigned = false;
}

// End of code generation (distorConEQ.cpp)
