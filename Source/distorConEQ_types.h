//
// distorConEQ_types.h
//
// Code generation for function 'onParamChangeCImpl'
//

#ifndef DISTORCONEQ_TYPES_H
#define DISTORCONEQ_TYPES_H

// Include files
#include "distorConEQ.h"
#include "rtwtypes.h"

// Type Definitions
struct distorConEQPersistentData {
  derivedAudioPlugin plugin;
  boolean_T plugin_not_empty;
  unsigned long long thisPtr;
  boolean_T thisPtr_not_empty;
};

struct distorConEQStackData {
  distorConEQPersistentData *pd;
};

#endif
// End of code generation (distorConEQ_types.h)
