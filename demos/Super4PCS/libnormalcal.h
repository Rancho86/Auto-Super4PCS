//
// MATLAB Compiler: 6.4 (R2017a)
// Date: Mon Aug 12 08:44:07 2019
// Arguments:
// "-B""macro_default""-W""cpplib:libnormalcal""-T""link:lib""normalcal.m"
//

#ifndef __libnormalcal_h
#define __libnormalcal_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "D:\Program Files\MATLAB\R2017a\extern\include/mclmcrrt.h"
#include "D:\Program Files\MATLAB\R2017a\extern\include/mclcppclass.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libnormalcal
#define PUBLIC_libnormalcal_C_API __global
#else
#define PUBLIC_libnormalcal_C_API /* No import statement needed. */
#endif

#define LIB_libnormalcal_C_API PUBLIC_libnormalcal_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libnormalcal
#define PUBLIC_libnormalcal_C_API __declspec(dllexport)
#else
#define PUBLIC_libnormalcal_C_API __declspec(dllimport)
#endif

#define LIB_libnormalcal_C_API PUBLIC_libnormalcal_C_API


#else

#define LIB_libnormalcal_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libnormalcal_C_API 
#define LIB_libnormalcal_C_API /* No special import/export declaration */
#endif

extern LIB_libnormalcal_C_API 
bool MW_CALL_CONV libnormalcalInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_libnormalcal_C_API 
bool MW_CALL_CONV libnormalcalInitialize(void);

extern LIB_libnormalcal_C_API 
void MW_CALL_CONV libnormalcalTerminate(void);



extern LIB_libnormalcal_C_API 
void MW_CALL_CONV libnormalcalPrintStackTrace(void);

extern LIB_libnormalcal_C_API 
bool MW_CALL_CONV mlxNormalcal(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__BORLANDC__)

#ifdef EXPORTING_libnormalcal
#define PUBLIC_libnormalcal_CPP_API __declspec(dllexport)
#else
#define PUBLIC_libnormalcal_CPP_API __declspec(dllimport)
#endif

#define LIB_libnormalcal_CPP_API PUBLIC_libnormalcal_CPP_API

#else

#if !defined(LIB_libnormalcal_CPP_API)
#if defined(LIB_libnormalcal_C_API)
#define LIB_libnormalcal_CPP_API LIB_libnormalcal_C_API
#else
#define LIB_libnormalcal_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_libnormalcal_CPP_API void MW_CALL_CONV normalcal(int nargout, mwArray& ptCloudOut, const mwArray& ptCloudIn, const mwArray& R);

#endif
#endif
