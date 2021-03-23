#include <vector>
#include "TROOT.h"
#include "TStyle.h"
//#include "TSystem.h"
#include "TString.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TFile.h"
//#include "TTree.h"
#include "TF1.h"
#include "Math/Functor.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/Interpolator.h"
#include "Math/WrappedParamFunction.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"
#include "Math/SpecFuncMathCore.h"

#define MOLAR_MASS 131.293
#define FIELD_MIN 1.

int recomb_test(){

    double energy_er[100];
    double energy_nr[100];
    double r_er[100];
    double deltar_er[100];
    double r_nr[100];
    double deltar_nr[100];

    for(int i = 0 ; i < 100 ; i++){
        energy_er[i] = i + 1.;
        energy_nr[i] = 0.2 * (i + 1);
    }

    //get ER array
    for(int i = 0 ; i < 100 ; i++){
    double gamma = 0.124;
    double omega = 31.;
    double delta = 0.24;
    double q0 = 1.13;
    double q1 = 0.47;
    double energy = energy_er[i];
    double E = 114.2
    double sigma = gamma * exp(-energy / omega) * pow(E, -delta);
    double q2 = 0.04;
    double q3 = 1.7;
    double rmean = (1. - log(1 + Ni * sigma/4. )/(Ni * sigma/4.))/(1 + exp(-(energy - q0)/q1));
    double deltaR = q2 * (1 - exp(-energy/q3));
    }

    TGraph * g0a = new TGraph();
    TGraph * g0b = new TGraph();
    TGraph * g1a = new TGraph();
    TGraph * g1b = new TGraph();
  return 0;
}
