#include "signal.h"

Signalcalc::Signalcalc(){Readin();};
Signalcalc::~Signalcalc(){};

void Signalcalc::Readin(){
    ReadB8NeutrinoSpectrum(E_B8,P_B8);
    ReadER(E_ER,F_ER);
    ReadNR(E_NR,F_NR);
    ReadAcceptance(E_accep,accep);
    ReadS1Eff(Nhit,Eff);
    ReadBLS(Npe,bias,fluc);
    //ReadNRUpperLimit(s1_up, s2_up);
    //ReadNRLowerLimit(s1_low,s2_low);
}

void Signalcalc::ReadB8NeutrinoSpectrum(double * x, double * y){
    double temp1, temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code5/B8.dat");
    for(int i = 0 ; i < 160; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();
}

void Signalcalc::ReadER(double * x, double * y){
    double temp1, temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code5/ER.dat");
    for(int i = 0; i < 6; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();
}

void Signalcalc::ReadNR(double * x, double * y){
    double temp1, temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code5/NR.dat");
    for(int i = 0; i < 89; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();
}

/*
void Signalcalc::ReadNRUpperLimit(double * x, double * y){
    double temp1, temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code3/code3/UpperLimitCut.dat");
    for(int i = 0; i < 69; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();
}

void Signalcalc::ReadNRLowerLimit(double * x, double * y){
    double temp1, temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code3/code3/LowerLimitCut.dat");
    for(int i = 0; i < 66; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();
}*/

void Signalcalc::ReadAcceptance(double * x, double *y){
    double temp1, temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code5/Acceptance.dat");
    for(int i = 0; i < 17; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();    
}

void Signalcalc::ReadS1Eff(int * x, double *y){
    int temp1;
    double temp2;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code5/S1_eff.dat");
    for(int i = 0; i < 20; i++){
        infile>>temp1>>temp2;
        x[i] = temp1;
        y[i] = temp2;
    }
    infile.close();    
}

void Signalcalc::ReadBLS(double * x,double * y,double * z){
    double temp1, temp2, temp3;
    ifstream infile("/Users/yunshan/work&study/sdu/neutrino/code5/bls.dat");
    for(int i = 0; i < 31; i++){
        infile>>temp1>>temp2>>temp3;
        x[i] = temp1;
        y[i] = temp2;
        z[i] = temp3;
    }
    infile.close();
}

double Signalcalc::GetB8NR(double E, double * param){

    //calculate flux
    double Flux_B8 = 5.69e6 ;//events per cm2 per day (SNO+2019)
    double temp;//y value of spectrum curve
    int n = sizeof(E_B8)/8;
    int i;

    /* find the interval */
    for(i = 0 ; i < n ; i++){
        if(E > E_B8[i]);
        else break;
    }
    if(i==0)temp = P_B8[i];
    else if(i==n)temp = P_B8[n-1];
    else{
        double Emin = E_B8[i-1];
        double Emax = E_B8[i];
        double P_Emin = P_B8[i-1];
        double P_Emax = P_B8[i];
        temp = calculation.GetLinearInte(E, Emin, Emax, P_Emin, P_Emax);
    }
    param[0] *= temp * Flux_B8;//E in MeV

    //let neutrino elastic scatter off a Xe target//
    double G = 1.1663787e-11;//in MeV-2
    double sin2theta = 0.23867;
    double Mxe = (double)calculation.SelectRanXeAtom();//in atomic number
    double N = Mxe - ATOM_NUM;
    double Q = ATOM_NUM * (4.0 * sin2theta - 1.0) + N;
    param[1] = Mxe;
    Mxe *= NEUC_MASS;//in MeV
    double NRmax = 2 * E * E / (Mxe + 2 * E);
    double Enf = gRandom->Uniform(0., NRmax);
    param[0] *= 16. * (NRmax - 0.) * Q * Q * G * G * Mxe * (1 - Mxe * Enf / (2 * E * E))/ (4 * M_PI)  ;
    param[0] *= 0.197e-10 * 0.197e-10;//convert MeV-2 to cm2
    param[0] *= 1.0e6/131.29 * 6.02e23;//number of xe per ton
    double e = 1000. * Enf;//in keV
    //if(e < 1.1)param[0] *= 0.;
    return e;
}

double Signalcalc::GetER(double * param){
    double Eer = gRandom->Uniform(0.,10.);//ER in keV
    double temp;//y value of spectrum curve
    int n = sizeof(E_ER)/8;
    int i;

    /* find the interval */
    for(i = 0 ; i < n ; i++){
        if(Eer > E_ER[i]);
        else break;
    }
    if(i==0)temp = F_ER[i];
    else if(i==n)temp = F_ER[n-1];
    else{
        double Emin = E_ER[i-1];
        double Emax = E_ER[i];
        double P_Emin = F_ER[i-1];
        double P_Emax = F_ER[i];
        temp = calculation.GetLinearInte(Eer, Emin, Emax, P_Emin, P_Emax);
    }
    param[0] *= (10. - 0.) * 0.05;
    return Eer;
}

double Signalcalc::GetNR(double * param){
    double Enr = gRandom->Uniform(0.,100.);//NR in keV
    double temp;//y value of spectrum curve
    int n = sizeof(E_NR)/8;
    int i;

    /* find the interval */
    for(i = 0 ; i < n ; i++){
        if(Enr > E_NR[i]);
        else break;
    }
    if(i==0)temp = F_NR[i];
    else if(i==n)temp = F_NR[n-1];
    else{
        double Emin = E_NR[i-1];
        double Emax = E_NR[i];
        double P_Emin = F_NR[i-1];
        double P_Emax = F_NR[i];
        temp = calculation.GetLinearInte(Enr, Emin, Emax, P_Emin, P_Emax);
    }
    param[0] *= (100.- 0.) * temp;
    param[1] = (double)calculation.SelectRanXeAtom();
    //if(Enr<1.1)param[0] *= 0.;
    return Enr;
}

double Signalcalc::GetAcceptance(double Enr){
    if (Enr < 0.5)return 0.;
    else{
        int n = sizeof(accep)/8;
        double temp;
        int i;
        for(i = 0 ; i < n ; i++){
            if(Enr > E_accep[i]);
            else break;
        }
        if(i==0)temp = accep[i];
        else if(i==n)temp = accep[n-1];
        else{
            double Emin = E_accep[i-1];
            double Emax = E_accep[i];
            double accep_min = accep[i-1];
            double accep_max = accep[i];
            temp = calculation.GetLinearInte(Enr, Emin, Emax, accep_min, accep_max);
        }
        
        return temp;
    }
}

double Signalcalc::GetS1Eff(int nhit){
    if (nhit < 3)return 0.;
    else{
        int n = sizeof(Nhit)/8;
        double temp;
        int i;
        for(i = 0 ; i < n ; i++){
            if(nhit == Nhit[i]){
                temp = Eff[i];
                break;
            }
        }
        if(n == i)temp = Eff[n-1];
       // std::cout<<temp<<std::endl;
        return temp;
    }
}

Bias Signalcalc::GetBLS(double npe){
    int n = sizeof(Npe)/8;
    double temp1,temp2;
    int i;
    for(i = 0 ; i < n ; i++){
        if(npe > Npe[i]);
        else break;
    }
    if(i==0){
        temp1 = bias[i];
        temp2 = fluc[i];
    }
    else if(i==n){
        temp1 = bias[n-1];
        temp2 = fluc[n-1];
    }
    else{
        double nmin = Npe[i-1];
        double nmax = Npe[i];
        double bias_nmin = bias[i-1];
        double bias_nmax = bias[i];
        double fluc_nmin = fluc[i-1];
        double fluc_nmax = fluc[i];

        temp1 = calculation.GetLinearInte(npe, nmin, nmax, bias_nmin, bias_nmax);
        temp2 = calculation.GetLinearInte(npe, nmin, nmax, fluc_nmin, fluc_nmax);
    }

    Bias result;
    result.bias = temp1;
    result.fluctuation = temp2;
    
    return result;
}

QuantaResult Signalcalc::GetQuanta_NR(double energy){
    //calculate lindhard factor//
    double epsilon = 11.5 * energy * pow(54, -7./3.);
    double g = 3 * pow(epsilon, 0.15) + 0.7 * pow(epsilon, 0.6) + epsilon;
    double kappa = 0.1394;
    double L = kappa * g / (1 + kappa * g);

    //fano fluctuations and lindhard fluctuations
    double Fano = 1.; 
    int Nq_mean = energy / W_EV;
    int Nq_actual = int(floor(gRandom->Gaus(Nq_mean, sqrt(Fano * Nq_mean)) + 0.5));
    if(Nq_actual < 0.)Nq_actual = 0;
    int Nq = calculation.BinomFluct(Nq_actual, L);

    //get exciton ratio and do fluctuations
    double alpha = 1.24;
    double zeta = 0.0472;
    double beta = 239.;
    double NexONi = alpha * pow(detector.get_E_drift(),-zeta) * (1 - exp(-beta * epsilon));//NexONi = Nex/Ni
    int Ni = calculation.BinomFluct(Nq,1/(1 + NexONi));
    int Nex = Nq - Ni;
    
    //get recombination fraction fluctuations
    double gamma = 0.01385;
    double delta = 0.0620;
    double sigma = gamma * pow(detector.get_E_drift(), -delta);
    double q2 = 0.04;
    double q3 = 1.7;
    double rmean = 1. - log(1 + Ni * sigma )/(Ni * sigma );
    double deltaR = q2 * (1 - exp(-energy/q3));
    double r = gRandom->Gaus(rmean,deltaR);
    if(r < 0.)r = 0.;
    
    //get photon and electron numbers
    int Ne = calculation.BinomFluct(Ni, 1 - r);
    int Nph = Ni + Nex - Ne;
        
    if (Nph > energy / W_SCINT)
    Nph = energy / W_SCINT;  // yields can never exceed 1 / [ W ~ few eV ]
    if (Ne > energy / W_SCINT) Ne = energy / W_SCINT;
    if (Nph < 0.) Nph = 0.;
    if (Ne < 0.) Ne = 0.;
    if (L < 0.) L = 0.;
    if (L > 1.) L = 1.;  // Lindhard Factor
    if (energy <  W_EV / L) {
        Nph = 0.;
        Ne = 0.;
    }

    QuantaResult result;
    result.photons = Nph;
    result.electrons = Ne;
    result.ions = Ni;
    result.excitons = Nex;
    result.truthz = gRandom->Uniform(0.,detector.get_TopDrift());
    result.lindhard = L;
    return result;     
}

QuantaResult Signalcalc::GetQuanta_ER(double energy){
    double L = 1.;

    int Nq_mean = energy / W_EV;
    double Fano = 0.12707 - 0.029623 * detector.get_density() - 0.0057042 * pow(detector.get_density(), 2.)
     + 0.0015957 * pow(detector.get_density(), 3.); 
    Fano += 0.0015 * sqrt(Nq_mean) * pow(detector.get_E_drift(), 0.5);
    int Nq_actual = int(floor(gRandom->Gaus(Nq_mean, sqrt(Fano * Nq_mean)) + 0.5));
    if(Nq_actual < 0.)Nq_actual = 0;
    int Nq = calculation.BinomFluct(Nq_actual, L);

    //get exciton ratio and do fluctuations
    double alpha = 0.067366 + detector.get_density() * 0.039693;
    double NexONi = alpha * erf(0.05 * energy);//NexONi = Nex/Ni
    int Ni = calculation.BinomFluct(Nq,1/(1 + NexONi));
    int Nex = Nq - Ni;
    
    //get recombination fraction fluctuations
    double gamma = 0.124;
    double omega = 31.;
    double delta = 0.24;
    double q0 = 1.13;
    double q1 = 0.47;
    double sigma = gamma * exp(-energy / omega) * pow(detector.get_E_drift(), -delta);
    double q2 = 0.04;
    double q3 = 1.7;
    double rmean = (1. - log(1 + Ni * sigma/4. )/(Ni * sigma/4.))/(1 + exp(-(energy - q0)/q1));
    double deltaR = q2 * (1 - exp(-energy/q3));
    double r = gRandom->Gaus(rmean,deltaR);
    if(r < 0.)r = 0.;
    
    //get photon and electron numbers
    int Ne = calculation.BinomFluct(Ni, 1 - r);
    int Nph = Ni + Nex - Ne;
        
    if (Nph > energy / W_SCINT)
    Nph = energy / W_SCINT;  // yields can never exceed 1 / [ W ~ few eV ]
    if (Ne > energy / W_SCINT) Ne = energy / W_SCINT;
    if (Nph < 0.) Nph = 0.;
    if (Ne < 0.) Ne = 0.;
    if (L < 0.) L = 0.;
    if (L > 1.) L = 1.;  // Lindhard Factor
    if (energy <  W_EV / L) {
        Nph = 0.;
        Ne = 0.;
    }

    QuantaResult result;
    result.photons = Nph;
    result.electrons = Ne;
    result.ions = Ni;
    result.excitons = Nex;
    result.truthz = gRandom->Uniform(0.,detector.get_TopDrift());
    result.lindhard = L;
    return result;  
}


vector<double> Signalcalc::GetS1(QuantaResult quanta) {
    vector<double> s1(9);
    int Nph = quanta.photons;

    double A = 5.64008e2;
    double B = 3.92408e-1;
    double C = -1.65968e-4;
    double D = 6.69171e-7;
    double dt = (detector.get_TopDrift() - quanta.truthz) / detector.get_driftvelocity();
    double qS1maxmean = 803.319;
    double qS1max = A + B * dt + C * dt * dt + D * dt * dt * dt;

    double g1mean = detector.get_g1();
    double g1_true = g1mean * qS1max/qS1maxmean;
    int nHits = calculation.BinomFluct(Nph, g1_true);

    double eff = GetS1Eff(nHits);
    double dice = gRandom->Rndm();
    int nHits_eff = nHits;
    if(dice > eff) nHits_eff *= 0.;
    double Nphe = nHits_eff + calculation.BinomFluct(nHits_eff, detector.get_P_dphe());
    double pulseArea = gRandom->Gaus(Nphe,detector.get_sPEres() * sqrt(Nphe));
    if(pulseArea < 0.)pulseArea = 0.;
       
    Bias BLS = GetBLS(pulseArea);
    double bias = gRandom->Gaus(BLS.bias,BLS.fluctuation);
    if(bias < 0.) bias = 0.;
    double pulseArea_bias = pulseArea * bias;

    double pulseArea_corr = pulseArea_bias * qS1maxmean/qS1max;

    s1[0] = Nph;//photon number
    s1[1] = nHits;//phd number
    s1[2] = nHits_eff;//effective phd number
    s1[3] = Nphe;//pe number
    s1[4] = pulseArea;//pe number after PMT resolution
    s1[5] = pulseArea_bias;//pe number after reconstruction bias
    s1[6] = pulseArea_corr;//pe number after position correction-z
    s1[7] = quanta.truthz;//truthz
    s1[8] = dt;//drift time
    
    return s1;
}

vector<double> Signalcalc::GetS2(QuantaResult quanta){
    vector<double> s2(9);

    int Ne = quanta.electrons;
    double dt = (detector.get_TopDrift() - quanta.truthz) / detector.get_driftvelocity();
    int Nee = calculation.BinomFluct(Ne, detector.get_ExtraEff() * exp(-dt / detector.get_eLife_us()));
    int nHits = long(floor(gRandom->Gaus(detector.get_SEG() * double(Nee), 
                sqrt(double(Nee)) * detector.get_deltaG()) + 0.5));
    if(nHits < 0.)nHits = 0;
    double Nphe = nHits + calculation.BinomFluct(nHits, detector.get_P_dphe());
    double pulseArea = gRandom->Gaus(Nphe,detector.get_sPEres() * sqrt(Nphe));
    if(pulseArea < 0.)pulseArea = 0.;
    Bias BLS = GetBLS(pulseArea);
    double bias = gRandom->Gaus(BLS.bias,BLS.fluctuation);
    if(bias < 0.)bias = 0.;
    double pulseArea_bias = pulseArea * bias;
    //if(pulseArea_bias < detector.get_s2_thr()) pulseArea_bias = 0. ;
    double pulseArea_corr = pulseArea_bias * exp(dt / detector.get_eLife_us());

    s2[0] = Ne;//electron number
    s2[1] = Nee;//extracted electron number
    s2[2] = nHits;//phd number
    s2[3] = Nphe;//pe number
    s2[4] = pulseArea;//pe number after PMT resolution
    s2[5] = pulseArea_bias;//pe number after reconstruction bias and analysis threshold
    s2[6] = pulseArea_corr;//pe number after position correction
    s2[7] = quanta.truthz;
    s2[8] = dt;
    
    return s2;
}


/*void Signalcalc::NRCut(double s1, double s2, double * param){
    double temp_up, temp_low;
    int n_up = sizeof(s1_up);
    int n_low = sizeof(s1_low);
    int i;

    for(i = 0 ; i < n_up ; i++){
        if(s1 > s1_up[i]);
        else break;
    }
    if(i==0)temp_up = s2_up[i];
    else{
        double s1min = s1_up[i-1];
        double s1max = s1_up[i];
        double s2min = s2_up[i-1];
        double s2max = s2_up[i];
        temp_up = calculation.GetLinearInte(s1, s1min, s1max, s2min, s2max);
    }

    for(i = 0 ; i < n_low ; i++){
        if(s1 > s1_low[i]);
        else break;
    }
    if(i==0)temp_low = s2_low[i];
    else{
        double s1min = s1_low[i-1];
        double s1max = s1_low[i];
        double s2min = s2_low[i-1];
        double s2max = s2_low[i];
        temp_low = calculation.GetLinearInte(s1, s1min, s1max, s2min, s2max);
    }
    if(s2>temp_low&&s2<temp_up) param[0] *= 1.;
    else param[0] *= 0.; 
}*/
