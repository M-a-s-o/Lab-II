/*
*   g++ main.cpp `root-config --cflags --glibs` -O3 -o main
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
    double T = 1./.8;
    return par[0]*(1-2/(1+exp(-T/(2*par[1])))*exp(-x[0]/par[1]))+par[2]; // V_C
    // par 0 = V_g,     par 1 = tau,    par 2 = A
}

// nfit = 4
Double_t func_RL (Double_t *x, Double_t *par) {
    double T = 1./100.;
    return par[0]*(2./(1+exp(-T/(2*par[1])))*exp(-x[0]/par[1]))+par[2]; // V_L
    // par 0 = V_g,     par 1 = tau,    par 2 = A
}

// nfit = 5
Double_t func_RL_real (Double_t *x, Double_t *par) {
    double T = 1./100., R = 992;
    return par[0]*(2*exp(-x[0]/par[2])/(1+exp(-T/(2*par[2])))+par[1]/(par[1]+R)*(1-2*exp(-x[0]/par[2])/(1+exp(-T/(2*par[2])))))+par[3]; // V_L
    // par 0 = V_g,     par 1 = R_L,    par 2 = tau,    par 3 = A
}

// nfit = 6
Double_t func_RLC_sott (Double_t *x, Double_t *par) {
    double R = 10.;
    return R*par[0]*exp(-par[1]*x[0])*sin(par[2]*x[0]+par[4])+par[3];
    // par 0 = I0,  par 1 = gamma,  par 2 = omega,  par 3 = A,  par 4 = B
    // par init:
    // 1.842207802921942 6.250772331335608 310.30877557486264 19.34482826950244 4.660252861131004
}

// nfit = 7
Double_t func_RLC_crit (Double_t *x, Double_t *par) {
    double R = 54.77;
    return R*par[0]*par[1]*par[1]*(x[0]+par[4])*exp(-par[1]*(x[0]+par[4]))+par[3];
    // par 0 = Q0,  par 1 = gamma,  par 3 = A,  par 4 = B
}

// nfit = 8
Double_t func_RLC_sovr (Double_t *x, Double_t *par) {
    double R = 1.;
    return R*par[0]*(par[1]*par[1]+par[2]*par[2])/(2*par[2])*(exp(-(par[1]-par[2])*(x[0]+par[4]))-exp(-(par[1]+par[2])*(x[0]+par[4])))+par[3];
    // par 0 = Q0,  par 1 = gamma,  par 2 = omega,  par 3 = A,  par 4 = B
}

int main (int argc, char *argv[]) {
    gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");
    const char *file_name = "./data.txt";
    int nfit = 0, npar = argc-3;
    double *init_par;
    if (argc >= 3) {
        //file_name = argv[1];
        nfit = atoi(argv[2]);
        if (nfit < 0 || nfit > 8) {
            std::cout << "nfit: 0 - A+Bx, 1 - Bx, 2 - Exp, 3 - RC, 4 - RL, 5 - RL irl, 6 - RLC sott, 7 - RLC crit, 8 - RLC sovr" << std::endl;
            return 1;
        }
        if (argc > 3) {
            init_par = new double[npar];
            for (int i = 0; i < npar; ++i)
                init_par[i] = atof(argv[i+3]);
        }
    }
    else {
        std::cout << "Passare: nome file, numero nfit, valori iniziali parametri" << std::endl;
        return 1;
    }

    std::vector<Double_t (*)(Double_t *, Double_t *)> func_names = {func_lin, func_lin_B, func_expo, func_RC, func_RL, func_RL_real, func_RLC_sott, func_RLC_crit, func_RLC_sovr};

    double func_min = 0., func_max = 1.;
    TF1 *fit = new TF1 ("fit", func_names.at(nfit), func_min, func_max, npar);
    fit->SetParameters(init_par);
    delete[] init_par;

    TGraphErrors grapherr (file_name);
    TFitResultPtr result = grapherr.Fit(fit, "S");

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
    canvas.Print("fit.png", "png");
    canvas.Clear();

    delete fit;

    return 0;
}