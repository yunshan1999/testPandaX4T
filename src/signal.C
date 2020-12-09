#include "signal.h"

Signalcalc::Signalcalc(){Readin();};
Signalcalc::~Signalcalc(){};

void Signalcalc::Readin(){
    ReadB8NeutrinoSpectrum(E_B8,P_B8);
    ReadER(E_ER,F_ER);
    ReadNR(E_NR,F_NR);
}

void Signalcalc::ReadB8NeutrinoSpectrum(double * x, double * y){
    double temp1, temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code3/code3/B8.dat");
    for(int i = 0 ; i < 160; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();
}

void Signalcalc::ReadER(double * x, double * y){
    double temp1, temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code3/code3/ER.dat");
    for(int i = 0; i < 6; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();
}

void Signalcalc::ReadNR(double * x, double * y){
    double temp1, temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code3/code3/NR.dat");
    for(int i = 0; i < 39; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();
}

double Signalcalc::GetB8NR(double E, double * param){
    double Flux_B8 = 5.95e6;//events per cm2 per second (SNO+2019)
    double dSigma = 1.25e-15;//parameter of diff cross section (A. Drukier and L. Stodolsky in 1984)
    double temp;//y value of spectrum curve
    int n = sizeof(E_B8);
    int i;

    /* find the interval */
    for(i = 0 ; i < n ; i++){
        if(E > E_B8[i]);
        else break;
    }
    if(i==0)temp = P_B8[i];
    else{
        double Emin = E_B8[i-1];
        double Emax = E_B8[i];
        double P_Emin = P_B8[i-1];
        double P_Emax = P_B8[i];
        temp = calculation.GetLinearInte(E, Emin, Emax, P_Emin, P_Emax);
    }
    param[0] *= temp * Flux_B8;//E in MeV

    //let neutrino elastic scatter off a Xe target//
    double Mxe = (double)calculation.SelectRanXeAtom();//in atomic number
    double N = Mxe - ATOM_NUM;
    param[1] = Mxe;
    Mxe *= NEUC_MASS;//in MeV
    double En = E;
    double cth = gRandom->Uniform(-1.0, 1.0);
    double Enf = En * En * (1 - cth) / (Mxe + En * (1 - cth));
    param[0] *= (16. * 2.) * dSigma * N * N * E * E * ( 1 + cth ) / MOLAR_MASS  ;
    return 1000. * Enf;//in keV
}

double Signalcalc::GetER(double * param){
    double Eer = gRandom->Uniform(0.,100.);//ER in keV
    double temp;//y value of spectrum curve
    int n = sizeof(E_ER);
    int i;

    /* find the interval */
    for(i = 0 ; i < n ; i++){
        if(Eer > E_ER[i]);
        else break;
    }
    if(i==0)temp = F_ER[i];
    else{
        double Emin = E_ER[i-1];
        double Emax = E_ER[i];
        double P_Emin = F_ER[i-1];
        double P_Emax = F_ER[i];
        temp = calculation.GetLinearInte(Eer, Emin, Emax, P_Emin, P_Emax);
    }
    param[0] *= (100.-0.) * temp;
    return Eer;
}

double Signalcalc::GetNR(double * param){
    double Enr = gRandom->Uniform(0.1,20.);//NR in keV
    double temp;//y value of spectrum curve
    int n = sizeof(E_NR);
    int i;

    /* find the interval */
    for(i = 0 ; i < n ; i++){
        if(Enr > E_NR[i]);
        else break;
    }
    if(i==0)temp = F_NR[i];
    else{
        double Emin = E_NR[i-1];
        double Emax = E_NR[i];
        double P_Emin = F_NR[i-1];
        double P_Emax = F_NR[i];
        temp = calculation.GetLinearInte(Enr, Emin, Emax, P_Emin, P_Emax);
    }
    param[0] *= (20.-0.1) * temp;
    return Enr;
}

YieldResult Signalcalc::GetYields_NR(double energy, double * param){
    double Param[11] = {11., 1.1, 0.0480, -0.0533, 12.6, 0.3, 2., 0.3, 2., 0.5, 1.};//model parameters
    double Ne = -999;
    double Nph = -999;
    double NexONi = -999;//NexONi = Nex/Ni
    double L = 1.;
    double massNumber = param[1];
    double ScaleFactor[2];

    ScaleFactor[0] = sqrt(MOLAR_MASS / massNumber);
    ScaleFactor[1] = ScaleFactor[0];
    double Nq = Param[0] * pow(energy, Param[1]);
    double ThomasImel = Param[2] * pow(detector.get_E_drift(), Param[3]) * pow(detector.get_density() / 2.8, 0.3);//2.8 g/ml is ref density from nest code
    double Qy = 1. / (ThomasImel * pow(energy + Param[4], Param[9]));
    Qy *= 1. - 1. / pow(1. + pow((energy / Param[5]), Param[6]), Param[10]);
    double Ly = Nq / energy - Qy;
    if (Qy < 0.0) Qy = 0.0;
    if (Ly < 0.0) Ly = 0.0;
    Ne = Qy * energy * ScaleFactor[1];
    Nph = Ly * energy * ScaleFactor[0] * (1. - 1. / (1. + pow((energy / Param[7]), Param[8])));
    Nq = Nph + Ne;
    double Ni = (4. / ThomasImel) * (exp(Ne * ThomasImel / 4.) - 1.);
    double Nex = (-1. / ThomasImel) * (4. * exp(Ne * ThomasImel / 4.) - (Ne + Nph) * ThomasImel - 4.);
    if (fabs(Nex - (Nq - Ni)) > PHE_MIN || fabs(Ni - (Nq - Nex)) > PHE_MIN) {
        std::cout << "\nERROR: Quanta not conserved. Tell Matthew Immediately!\n";
        exit(1);
        }
    NexONi = Nex / Ni;
    L = (Nq / energy) * W * 1e-3;      
          
    assert(Ne != -999 && Nph != -999 && NexONi != -999);
    if (Nph > energy / W_SCINT)
    Nph = energy / W_SCINT;  // yields can never exceed 1 / [ W ~ few eV ]
    if (Ne > energy / W_SCINT) Ne = energy / W_SCINT;
    if (Nph < 0.) Nph = 0.;
    if (Ne < 0.) Ne = 0.;
    // if (NexONi < 0.) NexONi = 0.;
    if (L < 0.) L = 0.;
    if (L > 1.) L = 1.;  // Lindhard Factor
    if (energy < 0.001 * W / L) {
        Nph = 0.;
        Ne = 0.;
    }

    YieldResult result;
    result.PhotonYield = Nph;
    result.ElectronYield = Ne;
    result.ExcitonRatio = NexONi;
    result.Lindhard = L;
    return result;     
}

YieldResult Signalcalc::GetYields_ER(double energy){
    double Ne = -999;
    double Nph = -999;
    double NexONi = -999;
    double alpha = 0.067366 + detector.get_density() * 0.039693;
    double L = 1.;
    double QyLvllowE = 1e3 / W + 6.5 * (1. - 1. / (1. + pow(detector.get_E_drift() / 47.408, 1.9851)));
    double HiFieldQy = 1. + 0.4607 / pow(1. + pow(detector.get_E_drift() / 621.74, -2.2717), 53.502);
    double QyLvlmedE = 32.988 - 32.988 /(1. + pow(detector.get_E_drift() / (0.026715 * exp(detector.get_density() / 0.33926)), 0.6705));
    QyLvlmedE *= HiFieldQy;
    double DokeBirks = 1652.264 + (1.415935e10 - 1652.264) / (1. + pow(detector.get_E_drift() / 0.02673144, 1.564691));
    double Nq = energy * 1e3 / W;  //
    double LET_power = -2.;
    double QyLvlhighE = 28.;
    //      if (density > 3.) QyLvlhighE = 49.; Solid Xe effect from Yoo. But,
    //      beware of enabling this line: enriched liquid Xe for neutrinoless
    //      double beta decay has density higher than 3g/cc;
    double Qy = QyLvlmedE + (QyLvllowE - QyLvlmedE) /pow(1. + 1.304 * pow(energy, 2.1393), 0.35535) + 
                QyLvlhighE / (1. + DokeBirks * pow(energy, LET_power));
    if (Qy > QyLvllowE && energy > 1. && detector.get_E_drift() > 1e4) Qy = QyLvllowE;
    double Ly = Nq / energy - Qy;
    Ne = Qy * energy;
    Nph = Ly * energy;
    NexONi = alpha * erf(0.05 * energy);
    assert(Ne != -999 && Nph != -999 && NexONi != -999);
    if (Nph > energy / W_SCINT)
    Nph = energy / W_SCINT;  // yields can never exceed 1 / [ W ~ few eV ]
    if (Ne > energy / W_SCINT) Ne = energy / W_SCINT;
    if (Nph < 0.) Nph = 0.;
    if (Ne < 0.) Ne = 0.;
    // if (NexONi < 0.) NexONi = 0.;
    if (L < 0.) L = 0.;
    if (L > 1.) L = 1.;  // Lindhard Factor
    if (energy < 0.001 * W / L) {
        Nph = 0.;
        Ne = 0.;
    }
    YieldResult result;
    result.PhotonYield = Nph;
    result.ElectronYield = Ne;
    result.ExcitonRatio = NexONi;
    result.Lindhard = L;
    return result;  
}

QuantaResult Signalcalc::GetQuanta(YieldResult yields){
    double Param[11] = {11., 1.1, 0.0480, -0.0533, 12.6, 0.3, 2., 0.3, 2., 0.5, 1.};//model parameters
    QuantaResult result;
    int Nq_actual, Ne, Nph, Ni, Nex;

    double NexONi = yields.ExcitonRatio, Fano = 1.;
    double Nq_mean = yields.PhotonYield + yields.ElectronYield;

    double elecFrac = yields.ElectronYield / Nq_mean;
    if (elecFrac > 1.) elecFrac = 1.;
    if (elecFrac < 0.) elecFrac = 0.;

    if (NexONi < 0.) NexONi = 0.;
    double alf = 1. / (1. + NexONi);
    double recombProb = 1. - (NexONi + 1.) * elecFrac;
    if (recombProb < 0.) NexONi = 1. / elecFrac - 1.;

    if (yields.Lindhard == 1.) {
        Fano = 0.12707 - 0.029623 * detector.get_density() - 0.0057042 * pow(detector.get_density(), 2.) + 0.0015957 * pow(detector.get_density(), 3.);  
                    // to get it to be ~0.03 for LXe (E Dahl Ph.D. thesis)
                    // Fano factor is  << 1//~0.1 for GXe w/ formula from Bolotnikov et al. 1995
        Nq_actual = int(floor(gRandom->Gaus(Nq_mean, sqrt(Fano * Nq_mean)) + 0.5));
        if (Nq_actual < 0 || Nq_mean == 0.) Nq_actual = 0;
        Ni = calculation.BinomFluct(Nq_actual, alf);
        Nex = Nq_actual - Ni;
    }
    else {
        Fano = 1.;
        Ni = int(floor(gRandom->Gaus(Nq_mean * alf, sqrt(Fano * Nq_mean * alf)) + 0.5));
        if (Ni < 0) Ni = 0;
        Fano = 1.;
        Nex = int(floor(gRandom->Gaus(Nq_mean * NexONi * alf, sqrt(Fano * Nq_mean * NexONi * alf)) + 0.5));
        if (Nex < 0) Nex = 0;
        Nq_actual = Nex + Ni;
    }

    if (Nq_actual == 0) {
      result.ions = 0;
      result.excitons = 0;
      result.photons = 0;
      result.electrons = 0;
      return result;
    }

    if (Nex < 0) Nex = 0;
    if (Ni < 0) Ni = 0;
    if (Nex > Nq_actual) Nex = Nq_actual;
    if (Ni > Nq_actual) Ni = Nq_actual;

    result.ions = Ni;
    result.excitons = Nex;

    if (Nex <= 0 ) recombProb = yields.PhotonYield / double(Ni);
    if (recombProb < 0.) recombProb = 0.;
    if (recombProb > 1.) recombProb = 1.;
    if (std::isnan(recombProb) || std::isnan(elecFrac) || Ni == 0 || recombProb == 0.0) {
        recombProb = 0.0;
        elecFrac = 1.0;
        result.photons = Nex;
        result.electrons = Ni;
        return result;
    }    
      double ef = detector.get_E_drift();
    double cc =
        0.075351 +
        (0.050461 - 0.075351) / pow(1. + pow(ef / 30057., 3.0008), 2.9832e5);
    if (cc < 0.) cc = 0.;
    double bb = 0.54;
    double aa = cc / pow(1. - bb, 2.);
    double omega = -aa * pow(recombProb - bb, 2.) + cc;
    if (omega < 0.0) omega = 0.0;
  
    if (yields.Lindhard < 1.)
      omega = Param[2] * exp(-pow(elecFrac - Param[3], 2.) / Param[4]);
    double Variance =
        recombProb * (1. - recombProb) * Ni + omega * omega * Ni * Ni;
    Ne = int(floor(gRandom->Gaus((1. - recombProb) * Ni, sqrt(Variance)) + 0.5));
    if (Ne < 0) Ne = 0;
    if (Ne > Ni) Ne = Ni;
  
    Nph = Nq_actual - Ne;
    if (Nph > Nq_actual) Nph = Nq_actual;
    if (Nph < Nex) Nph = Nex;
  
    if ((Nph + Ne) != (Nex + Ni)) {
      std::cerr << "\nERROR: Quanta not conserved. Tell Matthew Immediately!\n";
      exit(1);
    }
  
    result.photons = Nph;
    result.electrons = Ne;
  
    return result;  // quanta returned with recomb fluctuations
}

double Signalcalc::GetS1(QuantaResult quanta) {
    int Nph = quanta.photons;
    int nHits = calculation.BinomFluct(Nph, detector.get_g1()), Nphe = 0;
    Nphe = nHits + calculation.BinomFluct(nHits, detector.get_P_dphe());
    return Nphe;
}

double Signalcalc::GetS2(QuantaResult quanta){
    int Ne = quanta.electrons;
    double truthz = gRandom->Uniform(0.,detector.get_TopDrift());
    double dt = (detector.get_TopDrift() - truthz) / detector.get_driftvelocity();
    int Nee = calculation.BinomFluct(Ne, detector.get_ExtraEff() * exp(-dt / detector.get_eLife_us()));
    int nHits = long(floor(gRandom->Gaus(detector.get_SEG() * double(Nee), 
                sqrt(detector.get_s2Fano() * detector.get_SEG() * double(Nee) )) + 0.5));
    int Nphe = nHits + calculation.BinomFluct(nHits, detector.get_P_dphe());
    return Nphe;
}
