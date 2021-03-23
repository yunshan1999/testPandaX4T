#ifndef DETECTOR_H_
#define DETECTOR_H_

#include <vector>
#include "TROOT.h"
class Detector {

 private:
  double g1 ;
  double sPEres ;
  double sPEthr ;
  double P_dphe ;
  double coinWind ;
  int coinLevel ;
  int numPMTs ;
  double SEG ;
  double g2 ;
  double ExtraEff ;
  double s2Fano ;
  double deltaG ;
  double s2_thr ;
  double E_drift ;
  double eLife_us ;
  double driftvelocity ;
  double T_Kelvin ;
  double p_bar ;
  double density ;
  double radius ;
  double TopDrift ;
  double anode ;
  double gate ;
  double cathode ;

  public:
  Detector();
  ~Detector();
  void Initialization();

  // Primary Scintillation (S1) parameters
  double get_g1(); 
  double get_sPEres(); 
  double get_sPEthr(); 
  double get_P_dphe(); 

  double get_coinWind(); 
  int get_coinLevel(); 
  int get_numPMTs(); 

  // Ionization and Secondary Scintillation (S2) parameters
  double get_SEG(); 
  double get_g2(); 
  double get_ExtraEff();
  double get_s2Fano(); 
  double get_deltaG();
  double get_s2_thr(); 
  double get_E_drift(); 
  double get_eLife_us(); 
  double get_driftvelocity();

  // Thermodynamic Properties
  double get_T_Kelvin(); 
  double get_p_bar(); 
  double get_density();

  //Geometry
  double get_radius(); 
  double get_TopDrift(); 
  double get_anode(); 
  double get_gate(); 
  double get_cathode(); 

};

#endif
