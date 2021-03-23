#ifndef CALC_HH
#define CALC_HH 

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
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include <math.h>

using namespace std;

class calc {
 public:
  calc();
  ~calc();
  double rand_exponential(double half_life);
  double GetLinearInte(double x, double x1, double x2, double y1, double y2);
  double Trapezoidal(double * x, double * y);
  long BinomFluct(long N0, double p);
  int SelectRanXeAtom();

};

#endif
