/*
*   g++ main.cpp `root-config --cflags --glibs` -O3 -o main
*   compile with -ffast-math?
*/
#include <iostream>
#include <vector>

#include "TROOT.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCanvas.h"

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
    const double Vg = 1;
    const double R = 2002;
    return Vg/sqrt(1+(x[0]*R*par[0])*(x[0]*R*par[0]));
    // par 0 = C
}

// nfit = 10
Double_t func_trasf_RL (Double_t *x, Double_t *par) {
    const double Vg = 1;
    const double R = 992;
    return Vg*x[0]*par[0]/sqrt(R*R+(x[0]*par[0])*(x[0]*par[0]));
    // par 0 = L
}

// nfit = 11
Double_t func_trasf_RL_real (Double_t *x, Double_t *par) {
    const double Vg = 1;
    const double R = 992;
    return Vg*sqrt((par[1]*par[1]+(x[0]*par[0])*(x[0]*par[0]))/((R+par[1])*(R+par[1])+(x[0]*par[0])*(x[0]*par[0])));
    // par 0 = L,   par 1 = R_L
}

// nfit = 12
Double_t func_trasf_RLC_to_VR (Double_t *x, Double_t *par) {
    const double Vg = 1;
    const double R = 992;
    return Vg*R/sqrt(R*R+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])));
    // par 0 = C,   par 1 = L
}

// nfit = 13
Double_t func_trasf_RLC_to_VC (Double_t *x, Double_t *par) {
    const double Vg = 1;
    const double R = 992;
    return Vg/(x[0]*par[0])/sqrt(R*R+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])));
    // par 0 = C,   par 1 = L
}

// nfit = 14
Double_t func_trasf_RLC_to_VL (Double_t *x, Double_t *par) {
    const double Vg = 1;
    const double R = 992;
    return Vg*x[0]*par[1]/sqrt(R*R+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])));
    // par 0 = C,   par 1 = L
}

// nfit = 15
Double_t func_trasf_RLC_to_VR_real (Double_t *x, Double_t *par) {
    const double Vg = 1;
    const double R = 992;
    return Vg*R/sqrt((R*par[2])*(R+par[2])+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])));
    // par 0 = C,   par 1 = L,  par 2 = R_L
}

// nfit = 16
Double_t func_trasf_RLC_to_VC_real (Double_t *x, Double_t *par) {
    const double Vg = 1;
    const double R = 992;
    return Vg/(x[0]*par[0])/sqrt((R+par[2])*(R+par[2])+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0])));
    // par 0 = C,   par 1 = L,  par 2 = R_L
}

// nfit = 17
Double_t func_trasf_RLC_to_VL_real (Double_t *x, Double_t *par) {
    const double Vg = 1;
    const double R = 992;
    return Vg/sqrt((par[2]*par[2]+(x[0]*par[1])*(x[0]*par[1]))/((R+par[2])*(R+par[2])+(x[0]*par[1]-1/(x[0]*par[0]))*(x[0]*par[1]-1/(x[0]*par[0]))));
    // par 0 = C,   par 1 = L,  par 2 = R_L
}

int main (int argc, char *argv[]) {
    gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");
    const char *file_name;
    int nfit = 0, nskipargs = 3, npar = argc-nskipargs;
    double *init_par;
    if (argc >= nskipargs) {
        file_name = argv[1];// = "./data.txt";
        nfit = atoi(argv[2]);
        if (nfit < 0 || nfit > 17) {
            std::cout << "nfit: 0 - A+Bx, 1 - Bx, 2 - Exp, 3 - RC, 4 - RL, 5 - RL irl, 6 - RLC sott, 7 - RLC crit, 8 - RLC sovr" << std::endl;
            return 1;
        }
        if (npar == 0) {
            std::cout << "Passare parametri iniziali." << std::endl;
            return 1;
        }
        init_par = new double[npar];
        for (int i = 0; i < npar; ++i)
            init_par[i] = atof(argv[i+nskipargs]);
    }
    else {
        std::cout << "Passare: nome file, numero nfit, valori iniziali parametri" << std::endl;
        return 1;
    }

    std::vector<Double_t (*)(Double_t *, Double_t *)> func_names = {func_lin, func_lin_B, func_expo, func_RC, func_RL, func_RL_real, func_RLC_sott, func_RLC_crit, func_RLC_sovr, func_trasf_RC, func_trasf_RL, func_trasf_RL_real, func_trasf_RLC_to_VR, func_trasf_RLC_to_VC, func_trasf_RLC_to_VL, func_trasf_RLC_to_VR_real, func_trasf_RLC_to_VC_real, func_trasf_RLC_to_VL_real};

    double func_min = 0., func_max = 1.; // !
    TF1 *fit = new TF1 ("fit", func_names.at(nfit), func_min, func_max, npar); // min e max non cambiano il range del fit
    fit->SetParameters(init_par);
    delete[] init_par;

    TGraphErrors grapherr (file_name);
    //TFitResultPtr result = grapherr.Fit(fit, "SM", "", func_min, func_max);
    TFitResultPtr result = grapherr.Fit(fit, "SM");

    std::cout << "Fit terminato con successo: " << result->IsValid() << std::endl;
    std::cout << "Chi2: " << std::setprecision(15) << result->Chi2() << std::endl;
    std::cout << "Chi2 red: " << std::setprecision(15) << result->Chi2()/result->Ndf() << std::endl;
    std::cout << "Ndf: " << std::setprecision(15) << result->Ndf() << std::endl;
    std::cout << "Prob: " << std::setprecision(15) << result->Prob() << std::endl;

    for (int i = 0; i < npar; ++i)
        std::cout << fit->GetParName(i) << ": " << std::setprecision(15) << fit->GetParameter(i) << " +/- " << fit->GetParError(i) << std::endl; 
    std::cout << "Matrice di covarianza:" << std::endl;
    for (int i = 0; i < npar; ++i) {
        for (int j = 0; j < i+1; ++j)
            std::cout << std::setprecision(15) << (result->GetCovarianceMatrix())[i][j] << "\t"; 
        std::cout << std::endl;
    }
    std::cout << "Matrice di correlazione:" << std::endl;
    for (int i = 0; i < npar; ++i) {
        for (int j = 0; j < i+1; ++j)
            std::cout << std::setprecision(15) << (result->GetCorrelationMatrix())[i][j] << "\t"; 
        std::cout << std::endl;
    }

    TCanvas canvas ("canvas", "titolo", 2000, 1000);
    grapherr.Draw();
    canvas.Print("fit_cpp.png", "png");
    canvas.Clear();

    delete fit;
    return 0;
}