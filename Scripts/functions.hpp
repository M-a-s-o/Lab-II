#ifndef functions_main
#define functions_main

#include "RtypesCore.h"
#include <math.h>

// nfit = 0
Double_t func_lin (Double_t *x, Double_t *par) {
    return par[0]+x[0]*par[1];
}

// nfit = 1
Double_t func_lin_B (Double_t *x, Double_t *par) {
    return x[0]*par[0];
}

// nfit = 2
Double_t func_expo (Double_t *x, Double_t *par) {
    return par[0]*(exp(par[1]*x[0])-1)+par[2];
}

// nfit = 3
Double_t func_RC (Double_t *x, Double_t *par) {
    const double T = 1./.8;
    return par[0]*(1-2/(1+exp(-T/(2*par[1])))*exp(-x[0]/par[1]))+par[2]; // V_C
    // par 0 = V_g,     par 1 = tau,    par 2 = A
}

// nfit = 4
Double_t func_RL (Double_t *x, Double_t *par) {
    //double T = 1./100.*1000;
    const double T = 10.; // milli secondi
    return par[0]*(2./(1+exp(-T/(2*par[1])))*exp(-x[0]/par[1]))+par[2]; // V_L
    // par 0 = V_g,     par 1 = tau,    par 2 = A
}

// nfit = 5
Double_t func_RL_real (Double_t *x, Double_t *par) {
    const double T = 1./100., R = 992+50;
    return par[0]*(2*exp(-x[0]/par[2])/(1+exp(-T/(2*par[2])))+par[1]/(par[1]+R)*(1-2*exp(-x[0]/par[2])/(1+exp(-T/(2*par[2])))))+par[3]; // V_L
    // par 0 = V_g,     par 1 = R_L,    par 2 = tau,    par 3 = A
}

// nfit = 6
Double_t func_RLC_sott (Double_t *x, Double_t *par) {
    const double R = 1;// 500;
    return R*par[0]*exp(-par[1]*x[0])*sin(par[2]*x[0]+par[4])+par[3];
    // par 0 = I0,  par 1 = gamma,  par 2 = omega,  par 3 = A,  par 4 = B
}

// nfit = 7
Double_t func_RLC_crit (Double_t *x, Double_t *par) {
    const double R = 1;
    return R*par[1]*x[0]*exp(-par[0]*x[0])+par[2];
    // par 0 = gamma,  par 1 = A,  par 2 = B
}

// nfit = 8
Double_t func_RLC_sovr (Double_t *x, Double_t *par) {
    const double R = 1;
    //return R*par[0]*exp(-par[1]*x[0])*(par[4]*exp(par[2]*x[0])-par[5]*exp(-par[2]*x[0]))+par[3];
    // par 0 = I0,  par 1 = gamma,  par 2 = beta,   par 3 = A,  par 4 = B,  par 5 = C
    return R*par[0]*exp(par[1]*x[0])-R*par[2]*exp(par[3]*x[0])+par[4];
    // par 0 = A,   par 1 = B,  par 2 = C,  par 3 = D,  par 4 = E
}

// nfit = 9
Double_t func_trasf_RC (Double_t *x, Double_t *par) {
    const double R = 2002;
    return 1/sqrt(1+(x[0]*R*par[0])*(x[0]*R*par[0]))+par[1];
    // par 0 = C,   par 1 = A
}

// nfit = 10
Double_t func_trasf_RL (Double_t *x, Double_t *par) {
    const double R = 2002;
    return 1*x[0]*par[0]/sqrt(R*R+(x[0]*par[0])*(x[0]*par[0]))+par[1];
    // par 0 = L,   par 1 = A
}

// nfit = 11
Double_t func_trasf_RL_real (Double_t *x, Double_t *par) {
    const double R = 2002;
    return 1*sqrt((par[1]*par[1]+(x[0]*par[0])*(x[0]*par[0]))/((R+par[1])*(R+par[1])+(x[0]*par[0])*(x[0]*par[0])))+par[2];
    // par 0 = L,   par 1 = R_L,    par 2 = A
}

// nfit = 12
Double_t func_trasf_RLC_to_VR (Double_t *x, Double_t *par) {
    const double R = 2002;
    return R/sqrt((R+50)*(R+50)+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])))+par[2];
    // par 0 = C,   par 1 = L,  par 2 = A
}

// nfit = 13
Double_t func_trasf_RLC_to_VC (Double_t *x, Double_t *par) {
    const double R = 2002;
    return 1/(x[0]*par[0])/sqrt(R*R+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])))+par[2];
    // par 0 = C,   par 1 = L,  par 2 = A
}

// nfit = 14
Double_t func_trasf_RLC_to_VL (Double_t *x, Double_t *par) {
    const double R = 2002;
    return x[0]*par[1]/sqrt(R*R+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])))+par[2];
    // par 0 = C,   par 1 = L,  par 2 = A
}

// nfit = 15
Double_t func_trasf_RLC_to_VR_real (Double_t *x, Double_t *par) {
    const double R = 2002;
    return R/sqrt((R*par[2])*(R+par[2])+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])))+par[3];
    // par 0 = C,   par 1 = L,  par 2 = R_L,    par 3 = A
}

// nfit = 16
Double_t func_trasf_RLC_to_VC_real (Double_t *x, Double_t *par) {
    const double R = 2002;
    return 1/(x[0]*par[0])/sqrt((R+par[2])*(R+par[2])+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])))+par[3];
    // par 0 = C,   par 1 = L,  par 2 = R_L,    par 3 = A
}

// nfit = 17
Double_t func_trasf_RLC_to_VL_real (Double_t *x, Double_t *par) {
    const double R = 2002;
    return sqrt((par[2]*par[2]+(x[0]*par[1])*(x[0]*par[1]))/((R+par[2])*(R+par[2])+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0]))))+par[3];
    // par 0 = C,   par 1 = L,  par 2 = R_L,    par 3 = A
}

// nfit = 18
Double_t func_trasf_RLC_to_VR_fase (Double_t *x, Double_t *par) {
    const double R = 2002;
    return -atan2((x[0]*par[1]-1/(x[0]*par[0])), R)+par[2];
    // par 0 = C,   par 1 = L,  par 2 = A
}

// nfit = 19
Double_t func_trasf_RLC_to_VC_fase (Double_t *x, Double_t *par) {
    const double R = 2002;
    return -M_PI/2.-atan2((x[0]*par[1]-1/(x[0]*par[0])), R)+par[2];
    // par 0 = C,   par 1 = L,  par 2 = A
}

// nfit = 20
Double_t func_trasf_RLC_to_VL_fase (Double_t *x, Double_t *par) {
    const double R = 2002;
    return M_PI/2.-atan2((x[0]*par[1]-1/(x[0]*par[0])), R)+par[2];
    // par 0 = C,   par 1 = L,  par 2 = A
}

#endif