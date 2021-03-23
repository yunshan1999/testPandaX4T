#include "Detector.h"

using namespace std;

Detector::Detector() {Initialization(); }

Detector::~Detector() {}

void Detector::Initialization() {
  // Primary Scintillation (S1) parameters
  g1 = 0.09997;    // phd per S1 phot at dtCntr (not phe). Divide out 2-PE effect
  //g1 = 0.142;
  sPEres = 0.3;  // single phe resolution (Gaussian assumed)
  sPEthr = 0.35;  // POD threshold in phe, usually used IN PLACE of sPEeff
  P_dphe = 0.2;    // chance 1 photon makes 2 phe instead of 1 in Hamamatsu PMT

  coinWind = 100.;  // S1 coincidence window in ns
  coinLevel = 2;   // how many PMTs have to fire for an S1 to count
  numPMTs = 368;    // For coincidence calculation 169t 199b

  // Ionization and Secondary Scintillation (S2) parameters
  SEG = 28. ;
  g2 = 5.09017 * 4. ;
  //g2 = 11.4;
  ExtraEff = g2 / SEG;
  //ExtraEff = 0.9;
  s2Fano = 3.61;  // Fano-like fudge factor for SE width
  deltaG = SEG * 0.25;
  s2_thr = 80.;  // the S2 threshold in phe or PE, *not* phd. Affects NR most
  E_drift = 114.2;
  eLife_us = 600.;  // the drift electron mean lifetime in micro-seconds
  driftvelocity = 1.37824;

  // Thermodynamic Properties
  T_Kelvin = 177.;  // for liquid drift speed calculation
  p_bar = 2.06;     // gas pressure in units of bars, it controls S2 size
  density = 2.8611 ;

  radius = 600.;  // mm

  TopDrift = 1200.;  // mm 
  anode = 1209.;  
  gate = 1201.;  // mm. This is where the E-field changes (higher)
  cathode = 1.00;  // mm. Defines point below which events are gamma-X
}

double Detector::get_g1(){return g1;}
double Detector::get_sPEres(){return sPEres;}
double Detector::get_sPEthr(){return sPEthr;}
double Detector::get_P_dphe(){return P_dphe;}
double Detector::get_coinWind(){return coinWind;}
int Detector::get_coinLevel(){return coinLevel;}
int Detector::get_numPMTs(){return numPMTs;}
double Detector::get_SEG(){return SEG;}
double Detector::get_g2(){return g2;}
double Detector::get_ExtraEff(){return ExtraEff;}
double Detector::get_s2Fano(){return s2Fano;}
double Detector::get_deltaG(){return deltaG;}
double Detector::get_s2_thr(){return s2_thr;}
double Detector::get_E_drift(){return E_drift;}
double Detector::get_eLife_us(){return eLife_us;}
double Detector::get_driftvelocity(){return driftvelocity;}
double Detector::get_T_Kelvin(){return T_Kelvin;}
double Detector::get_p_bar(){return p_bar;}
double Detector::get_density(){return density;}
double Detector::get_radius(){return radius;}
double Detector::get_TopDrift(){return TopDrift;}
double Detector::get_anode(){return anode;}
double Detector::get_gate(){return gate;}
double Detector::get_cathode(){return cathode;}


