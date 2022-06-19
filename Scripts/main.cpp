/*
*   g++ main.cpp functions.hpp `root-config --cflags --glibs` -O3 -o main
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

#include "functions.hpp"

int main (int argc, char *argv[]) {
    gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");
    const char *file_name;
    int nfit = 0, nskipargs = 3, npar = argc-nskipargs;
    double *init_par;
    std::vector<Double_t (*)(Double_t *, Double_t *)> *func_names = new std::vector<Double_t (*)(Double_t *, Double_t *)> {\
    func_lin, func_lin_B, func_expo, \
    func_RC, func_RL, func_RL_real, \
    func_RLC_sott, func_RLC_crit, func_RLC_sovr, \
    func_trasf_RC, func_trasf_RL, func_trasf_RL_real, \
    func_trasf_RLC_to_VR, func_trasf_RLC_to_VC, func_trasf_RLC_to_VL, \
    func_trasf_RLC_to_VR_real, func_trasf_RLC_to_VC_real, func_trasf_RLC_to_VL_real, \
    func_trasf_RLC_to_VR_fase, func_trasf_RLC_to_VC_fase, func_trasf_RLC_to_VL_fase, \
    func_dist_spec, func_indice_aria, func_righello, \
    func_Fraunhofer, func_Cauchy};

    if (argc >= nskipargs) {
        file_name = argv[1];// = "./data.txt";
        nfit = atoi(argv[2]);
        if (nfit < 0 || nfit > func_names->size()-1) {
            std::cout << "Fit permessi: da 0 a " << func_names->size()-1 << "." << std::endl;
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

    double func_min = 0., func_max = 1.; // !
    TF1 *fit = new TF1 ("fit", func_names->at(nfit), func_min, func_max, npar); // min e max non cambiano il range del fit
    delete func_names;
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
    delete fit;

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

    return 0;
}