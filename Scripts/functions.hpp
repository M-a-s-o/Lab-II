#ifndef functions_main
#define functions_main

#include "RtypesCore.h"
#include <math.h>

Double_t func_lin (Double_t *x, Double_t *par);
Double_t func_lin_B (Double_t *x, Double_t *par);
Double_t func_expo (Double_t *x, Double_t *par);
Double_t func_RC (Double_t *x, Double_t *par);
Double_t func_RL (Double_t *x, Double_t *par);
Double_t func_RL_real (Double_t *x, Double_t *par);
Double_t func_RLC_sott (Double_t *x, Double_t *par);
Double_t func_RLC_crit (Double_t *x, Double_t *par);
Double_t func_RLC_sovr (Double_t *x, Double_t *par);
Double_t func_trasf_RC (Double_t *x, Double_t *par);
Double_t func_trasf_RL (Double_t *x, Double_t *par);
Double_t func_trasf_RL_real (Double_t *x, Double_t *par);
Double_t func_trasf_RLC_to_VR (Double_t *x, Double_t *par);
Double_t func_trasf_RLC_to_VC (Double_t *x, Double_t *par);
Double_t func_trasf_RLC_to_VL (Double_t *x, Double_t *par);
Double_t func_trasf_RLC_to_VR_real (Double_t *x, Double_t *par);
Double_t func_trasf_RLC_to_VC_real (Double_t *x, Double_t *par);
Double_t func_trasf_RLC_to_VL_real (Double_t *x, Double_t *par);
Double_t func_trasf_RLC_to_VR_fase (Double_t *x, Double_t *par);
Double_t func_trasf_RLC_to_VC_fase (Double_t *x, Double_t *par);
Double_t func_trasf_RLC_to_VL_fase (Double_t *x, Double_t *par);
Double_t func_dist_spec (Double_t *x, Double_t *par);
Double_t func_indice_aria (Double_t *x, Double_t *par);
Double_t func_righello (Double_t *x, Double_t *par);

#endif