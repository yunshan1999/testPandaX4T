#ifndef SIGNAL_H_
#define SIGNAL_H_ 

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "calc.h"
#include "Detector.h"

using namespace std;

#define W_EV 13.7e-3
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
    double truthz;
    double lindhard;
};

struct Bias{
    double bias;
    double fluctuation;
};

class Signalcalc{
    private:
        Detector detector;
        calc calculation;
        // prob distrubiton function of B8 solar neutrino and Flux of ER and NR
        double E_B8[160],P_B8[160];
        double E_ER[6],F_ER[6];
        double E_NR[89],F_NR[89];
        //double s1_up[69],s2_up[69];
        //double s1_low[66],s2_low[66];
        double accep[17],E_accep[17];
        int Nhit[20];
        double Eff[20];
        double Npe[31],bias[31],fluc[31];

    public:
        Signalcalc();
        ~Signalcalc();

        void Readin();

        void ReadB8NeutrinoSpectrum(double * x, double * y);
        void ReadER(double * x, double * y);
        void ReadNR(double * x, double * y);
        void ReadAcceptance(double * x, double * y);
        void ReadS1Eff(int * x, double * y);
        void ReadBLS(double * x, double * y, double * z);
        //void ReadNRUpperLimit(double * x, double * y);
        //void ReadNRLowerLimit(double * x, double * y);
        double GetB8NR(double E, double * param);
        double GetER(double * param);
        double GetNR(double * param);
        double GetAcceptance(double Enr);
        double GetS1Eff(int nhit);
        Bias GetBLS(double Npe);
        QuantaResult GetQuanta_NR(double energy);
        QuantaResult GetQuanta_ER(double energy);
        vector<double> GetS1(QuantaResult quanta);
        vector<double> GetS2(QuantaResult quanta);
       // void NRCut(double s1, double s2, double * param);
};

#endif