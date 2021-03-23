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

int calc_driftvelocity(){
    double Kelvin = 177.;
    double bara = 2.06;
    double eField = 114.2;
  double speed =
      0.0;  // returns drift speed in mm/usec. based on Fig. 14 arXiv:1712.08607
  int i, j;
  double vi, vf, slope, Ti, Tf, offset;

  double polyExp[11][7] = {
      {-3.1046, 27.037, -2.1668, 193.27, -4.8024, 646.04, 9.2471},  // 100K
      {-2.7394, 22.760, -1.7775, 222.72, -5.0836, 724.98, 8.7189},  // 120
      {-2.3646, 164.91, -1.6984, 21.473, -4.4752, 1202.2, 7.9744},  // 140
      {-1.8097, 235.65, -1.7621, 36.855, -3.5925, 1356.2, 6.7865},  // 155
      {-1.5000, 37.021, -1.1430, 6.4590, -4.0337, 855.43,
       5.4238},  // 157, merging Miller with Yoo
      {-1.4939, 47.879, 0.12608, 8.9095, -1.3480, 1310.9,
       2.7598},  // 163, merging Miller with Yoo
      {-1.5389, 26.602, -.44589, 196.08, -1.1516, 1810.8, 2.8912},  // 165
      {-1.5000, 28.510, -.21948, 183.49, -1.4320, 1652.9, 2.884},   // 167
      {-1.1781, 49.072, -1.3008, 3438.4, -.14817, 312.12, 2.8049},  // 184
      {1.2466, 85.975, -.88005, 918.57, -3.0085, 27.568, 2.3823},   // 200
      {334.60, 37.556, 0.92211, 345.27, -338.00, 37.346, 1.9834}};  // 230

  double Temperatures[11] = {100., 120., 140., 155., 157., 163.,
                             165., 167., 184., 200., 230.};

  if (Kelvin >= Temperatures[0] && Kelvin < Temperatures[1])
    i = 0;
  else if (Kelvin >= Temperatures[1] && Kelvin < Temperatures[2])
    i = 1;
  else if (Kelvin >= Temperatures[2] && Kelvin < Temperatures[3])
    i = 2;
  else if (Kelvin >= Temperatures[3] && Kelvin < Temperatures[4])
    i = 3;
  else if (Kelvin >= Temperatures[4] && Kelvin < Temperatures[5])
    i = 4;
  else if (Kelvin >= Temperatures[5] && Kelvin < Temperatures[6])
    i = 5;
  else if (Kelvin >= Temperatures[6] && Kelvin < Temperatures[7])
    i = 6;
  else if (Kelvin >= Temperatures[7] && Kelvin < Temperatures[8])
    i = 7;
  else if (Kelvin >= Temperatures[8] && Kelvin < Temperatures[9])
    i = 8;
  else if (Kelvin >= Temperatures[9] && Kelvin <= Temperatures[10])
    i = 9;
  else {
    std::cout << "\nERROR: TEMPERATURE OUT OF RANGE (100-230 K)\n"<<std::endl;
    exit(1);
  }

  j = i + 1;
  Ti = Temperatures[i];
  Tf = Temperatures[j];
  // functional form from http://zunzun.com
  vi = polyExp[i][0] * exp(-eField / polyExp[i][1]) +
       polyExp[i][2] * exp(-eField / polyExp[i][3]) +
       polyExp[i][4] * exp(-eField / polyExp[i][5]) + polyExp[i][6];
  vf = polyExp[j][0] * exp(-eField / polyExp[j][1]) +
       polyExp[j][2] * exp(-eField / polyExp[j][3]) +
       polyExp[j][4] * exp(-eField / polyExp[j][5]) + polyExp[j][6];
  if (Kelvin == Ti) std::cout<<vi<<std::endl;
  if (Kelvin == Tf) std::cout<<vf<<std::endl;
  if (vf < vi) {
    offset = (sqrt((Tf * (vf - vi) - Ti * (vf - vi) - 4.) * (vf - vi)) +
              sqrt(Tf - Ti) * (vf + vi)) /
             (2. * sqrt(Tf - Ti));
    slope = -(sqrt(Tf - Ti) *
                  sqrt((Tf * (vf - vi) - Ti * (vf - vi) - 4.) * (vf - vi)) -
              (Tf + Ti) * (vf - vi)) /
            (2. * (vf - vi));
    speed = 1. / (Kelvin - slope) + offset;
  } else {
    slope = (vf - vi) / (Tf - Ti);
    speed = slope * (Kelvin - Ti) + vi;
  }

  if (speed <= 0.) {
    if (eField < 1e2 && eField >= FIELD_MIN) {
      std::cerr << "\nERROR: DRIFT SPEED NON-POSITIVE -- FIELD TOO LOW\n"<<std::endl;
      exit(1);
    }
    if (eField > 1e4) {
      std::cerr << "\nERROR: DRIFT SPEED NON-POSITIVE -- FIELD TOO HIGH\n"<<std::endl;
      exit(1);
    }
  }
  std::cout<<speed<<std::endl;
  return 0;
}
