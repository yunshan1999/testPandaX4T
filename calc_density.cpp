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

int calc_density(){
    double Kelvin = 177.;
    double bara = 2.06;
    if (Kelvin < 161.40) {  // solid Xenon
    std::cout << "\nWARNING: SOLID PHASE. IS THAT WHAT YOU WANTED?\n"<<std::endl;
     std::cout<<3.41<<std::endl;  // from Yoo at 157K
    // other sources say 3.100 (Wikipedia, 'maximum') and 3.64g/mL at an unknown
    // temperature
  }

  double VaporP_bar;  // we will calculate using NIST
  if (Kelvin < 289.7)
    VaporP_bar = pow(10., 4.0519 - 667.16 / Kelvin);
  else
    VaporP_bar = DBL_MAX;
  if (bara < VaporP_bar) {
    double density =
        bara * 1e5 / (Kelvin * 8.314);  // ideal gas law approximation, mol/m^3
    density *= MOLAR_MASS * 1e-6;
    std::cerr << "\nWARNING: GAS PHASE. IS THAT WHAT YOU WANTED?\n"<<std::endl;
    return density;  // in g/cm^3
  }

  double density =  2.9970938084691329E+02 * exp(-8.2598864714323525E-02 * Kelvin) -
         1.8801286589442915E+06 * exp(-pow((Kelvin - 4.0820251276172212E+02) /
                                               2.7863170223154846E+01,
                                           2.)) -
         5.4964506351743057E+03 * exp(-pow((Kelvin - 6.3688597345042672E+02) /
                                               1.1225818853661815E+02,
                                           2.)) +
         8.3450538370682614E+02 * exp(-pow((Kelvin + 4.8840568924597342E+01) /
                                               7.3804147172071107E+03,
                                           2.)) -
         8.3086310405942265E+02;  // in grams per cubic centimeter based on
                                  // zunzun fit to NIST data; will add gas later
    std::cout<<density<<std::endl;
    return 0;
}