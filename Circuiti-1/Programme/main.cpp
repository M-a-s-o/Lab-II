#include <iostream>

#include "TROOT.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCanvas.h"

Double_t func_expo (Double_t *x, Double_t *par) {
    return par[0]*(exp(par[1]*x[0])-1)+par[2];
}

Double_t func_lin (Double_t *x, Double_t *par) {
    return par[0]+x[0]*par[1];
}

Double_t func_lin_B (Double_t *x, Double_t *par) {
    return x[0]*par[0];
}

/*void experiment (int N, std::vector<double> &xs, std::vector<double> &ys, std::vector<double> &xe, std::vector<double> &ye, double *theta) {
    xs.clear(); ys.clear(); xe.clear(); ye.clear();
    for (int i = 0; i < N; ++i) {
        double xerr = RNGGauss_ziggurat(0, .01);
        xs.push_back((i+1)/10.+xerr);
        if (xerr < 0)
            xerr = -xerr;
        xe.push_back(xerr);
        double err = RNGGauss_ziggurat(0, .1);
        if (err < 0)
            err = -err;
        ys.push_back(func_expo(&xs.at(i), theta) + err);
        ye.push_back(err);
    }
}*/

int main (int argc, char *argv[]) {
    gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");
    const char *file_name;
    int nfit = 0;
    double *init_par;
    if (argc >= 3) {
        file_name = argv[1];
        nfit = atoi(argv[2]);
        if (nfit < 0 || nfit > 2) {
            std::cout << "nfit: 0 - A+Bx, 1 - Bx, 2 - Exp" << std::endl;
            return 1;
        }
        if (argc > 3) {
        init_par = new double[argc-3];
        for (int i = 0; i < argc-3; ++i)
                init_par[i] = atof(argv[i+3]);
        }
    }
    else {
        std::cout << "Passare: nome file, numero nfit, valori iniziali parametri" << std::endl;
        return 1;
    }

    TGraphErrors grapherr (file_name);
    TF1 *fit;

    if (nfit == 2) {
        fit = new TF1 ("fit", func_expo, 0., 1., 3);
        fit->SetParameters(init_par);
    }
    else if (nfit == 0) {
        fit = new TF1 ("fit", func_lin, 0., 1., 2);
        fit->SetParameters(init_par);
    }
    else if (nfit == 1) {
        fit = new TF1 ("fit", func_lin_B, 0., 1., 1);
        fit->SetParameters(init_par);
    }
    delete[] init_par;

    TFitResultPtr result = grapherr.Fit(fit, "S");

    std::cout << "Fit terminato con successo: " << result->IsValid() << std::endl;
    std::cout << "Chi2: " << std::setprecision(15) << result->Chi2() << std::endl;
    std::cout << "Chi2 red: " << std::setprecision(15) << result->Chi2()/result->Ndf() << std::endl;
    std::cout << "Ndf: " << std::setprecision(15) << result->Ndf() << std::endl;
    std::cout << "Prob: " << std::setprecision(15) << result->Prob() << std::endl;

    for (int i = 0; i < argc-3; ++i) {
        std::cout << fit->GetParName(i) << ": " << std::setprecision(15) << fit->GetParameter(i) << " +/- " << fit->GetParError(i) << std::endl; 
    }
    std::cout << "Matrice di covarianza:" << std::endl;
    for (int i = 0; i < argc-3; ++i) {
        for (int j = 0; j < i+1; ++j)
            std::cout << std::setprecision(15) << (result->GetCovarianceMatrix())[i][j] << "\t"; 
        std::cout << std::endl;
    }
    std::cout << "Matrice di correlazione:" << std::endl;
    for (int i = 0; i < argc-3; ++i) {
        for (int j = 0; j < i+1; ++j)
            std::cout << std::setprecision(15) << (result->GetCorrelationMatrix())[i][j] << "\t"; 
        std::cout << std::endl;
    }

    TCanvas canvas ("canvas", "titolo", 1280, 720);
    grapherr.Draw();
    canvas.Print("fit.png", "png");
    canvas.Clear();

    delete fit;

    return 0;
}