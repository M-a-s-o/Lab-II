#ifndef functions_main
#define functions_main

#include "RtypesCore.h"
#include <math.h>

typedef Double_t func(Double_t *, Double_t *);

func func_lin, func_lin_B, func_expo, \
func_RC, func_RL, func_RL_real, \
func_RLC_sott, func_RLC_crit, func_RLC_sovr, \
func_trasf_RC, func_trasf_RL, func_trasf_RL_real, \
func_trasf_RLC_to_VR, func_trasf_RLC_to_VC, func_trasf_RLC_to_VL, \
func_trasf_RLC_to_VR_real, func_trasf_RLC_to_VC_real, func_trasf_RLC_to_VL_real, \
func_trasf_RLC_to_VR_fase, func_trasf_RLC_to_VC_fase, func_trasf_RLC_to_VL_fase, \
func_dist_spec, func_indice_aria, func_righello, \
func_Fraunhofer, func_Cauchy;

#endif