#ifndef SIGNAL_H_
#define SIGNAL_H_ 

#include <iostream>
#include <fstream>
#include <cmath>

#include "calc.h"
#include "Detector.h"

using namespace std;

#define W 13.7
#define W_SCINT 8.5e-3
#define AVO 6.0221409e+23
#define ATOM_NUM 54.
#define MOLAR_MASS 131.293
#define ELEC_MASS 9.109e-31
#define NEUC_MASS 939. //neucleon's mass in MeV
#define FIELD_MIN 1.
#define PHE_MIN 1e-6

struct YieldResult{
    double PhotonYield;
    double ElectronYield;
    double ExcitonRatio;
    double Lindhard;
};

struct QuantaResult{
    int photons;
    int electrons;
    int ions;
    int excitons;
};

class Signalcalc{
    private:
        Detector detector;
        calc calculation;
        // prob distrubiton function of B8 solar neutrino and Flux of ER and NR
        double E_B8[160],P_B8[160];
        double E_ER[6],F_ER[6];
        double E_NR[39],F_NR[39];

    public:
        Signalcalc();
        ~Signalcalc();
        void Readin();

        void ReadB8NeutrinoSpectrum(double * x, double * y);
        void ReadER(double * x, double * y);
        void ReadNR(double * x, double * y);
        double GetB8NR(double E, double * param);
        double GetER(double * param);
        double GetNR(double * param);
        YieldResult GetYields_NR(double energy, double * param);
        YieldResult GetYields_ER(double energy);
        QuantaResult GetQuanta(YieldResult yields);
        double GetS1(QuantaResult quanta);
        double GetS2(QuantaResult quanta);
};

#endif