#include "calc.h"

using namespace std;

calc::calc() {}

calc::~calc() {}

double calc::rand_exponential(double half_life){
    double r = gRandom->Uniform(0.,1.);
    return log(1 - r) * (-1) * half_life / log(2.);
}

double calc::GetLinearInte(double x, double x1, double x2, double y1, double y2){
    return (y2 - y1) * (x - x1) / (x2 - x1) + y1;
}

double calc::Trapezoidal(double * x, double * y){
    double n = sizeof(x);
    double Integral = 0;
    for(int i = 0; i < n-1; i++){
        double h = x[i+1] - x[i];
        Integral += 0.5 * h * (y[i+1] + y[i]);
        }
    return Integral;
}

long calc::BinomFluct(long N0, double prob) {
    double mean = N0 * prob;
    double sigma = sqrt(N0 * prob * (1. - prob));
    int N1 = 0;
  
    if (prob <= 0.00) return N1;
    if (prob >= 1.00) return N0;
  
    if (N0 < 10) {
        for (int i = 0; i < N0; i++) {
            if (gRandom->Uniform(0.,1.) < prob) N1++;
        }
    } else {
      N1 = int(floor(gRandom->Gaus(mean,sigma) + 0.5));
    }
  
    if (N1 > N0) N1 = N0;
    if (N1 < 0) N1 = 0;
  
    return N1;
}

int calc::SelectRanXeAtom(){

  int A;
  double isotope = gRandom->Uniform(0.,1.) * 100.;

  if (isotope > 0.000 && isotope <= 0.090)
    A = 124;
  else if (isotope > 0.090 && isotope <= 0.180)
    A = 126;
  else if (isotope > 0.180 && isotope <= 2.100)
    A = 128;
  else if (isotope > 2.100 && isotope <= 28.54)
    A = 129;
  else if (isotope > 28.54 && isotope <= 32.62)
    A = 130;
  else if (isotope > 32.62 && isotope <= 53.80)
    A = 131;
  else if (isotope > 53.80 && isotope <= 80.69)
    A = 132;
  else if (isotope > 80.69 && isotope <= 91.13)
    A = 134;
  else
    A = 136;
  return A;
}