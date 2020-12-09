#include <iostream>
#include "signal.h"

int main(){
    Long64_t Nsim = 10000000;
    TString filename = "testPandaX4T.root";
    TFile * fs = new TFile(filename.Data(), "RECREATE");

    TLorentzVector ki[2],kf[2];
    double params[2];//[0] is weight, [1] is Xe mass number  
	double year = 3600. * 24. * 365.;

	TH1D * h_test = new TH1D("8B_Event_Rate", "",40, 0., 4.); //draw 8B spectrum to test whether weight is correct
    TH2D * h0a = new TH2D("8B_log10(S2/S1)_vs_S1","",50,0.,50.,50,0.,5.); 
	TH2D * h0b = new TH2D("8B_log10S2_vs_S1","",50,0.,50.,50,0.,5.); 
	TH1D * h0c = new TH1D("8B_without_signal_energy","",40, 0., 4.);
	TH1D * h0d = new TH1D("8B_with_signal_energy","",40, 0., 4.);
	TH2D * h1a = new TH2D("ER_log10(S2/S1)_vs_S1","",50,0.,50.,50,0.,5.);
	TH2D * h1b = new TH2D("ER_log10S2_vs_S1","",50,0.,50.,50,0.,5.);
	TH1D * h1c = new TH1D("ER_without_signal_energy","",100, 0., 100.);
	TH1D * h1d = new TH1D("ER_with_signal_energy","",100, 0., 100.);
	TH2D * h2a = new TH2D("NR_log10(S2/S1)_vs_S1","",50,0.,50.,50,0.,5.);
	TH2D * h2b = new TH2D("NR_log10S2_vs_S1","",50,0.,50.,50,0.,5.);
	TH1D * h2c = new TH1D("NR_without_signal_energy","",20, 0., 20.);
	TH1D * h2d = new TH1D("NR_with_signal_energy","",20, 0., 20.);

	TH1D * h1 = new TH1D("h1","",100,0,100);
	TH1D * h0 = new TH1D("h0","",100,0,20);

	h_test->SetDirectory(fs);
	h0a->SetDirectory(fs);
	h0b->SetDirectory(fs);
	h0c->SetDirectory(fs);
	h0d->SetDirectory(fs);
	h1a->SetDirectory(fs);
	h1b->SetDirectory(fs);
	h1c->SetDirectory(fs);
	h1d->SetDirectory(fs);
	h2a->SetDirectory(fs);
	h2b->SetDirectory(fs);
	h2c->SetDirectory(fs);
	h2d->SetDirectory(fs);
	h1 ->SetDirectory(fs);
	h0 ->SetDirectory(fs);
	
	Signalcalc signal;

    for(Long64_t i = 0; i < Nsim ; i++){//8B event by event
        if(i%1000000==0)std::cout<<i<<std::endl;
        params[0] = 1.0;
		double En = gRandom->Uniform(0., 16.);//in MeV
		double Enr = signal.GetB8NR(En, params);// in keV
        h_test->Fill(Enr,params[0]);
		YieldResult yields;
        QuantaResult quanta;
		yields = signal.GetYields_NR(Enr,params);
		quanta = signal.GetQuanta(yields);
		int S1Npe = signal.GetS1(quanta);
		int S2Npe = signal.GetS2(quanta);
		if(S2Npe!=0){
			h0b->Fill(S1Npe,log10(S2Npe),params[0]);
			if(S1Npe!=0){
				h0a->Fill(S1Npe,log10(S2Npe/S1Npe),params[0]);
				h0d->Fill(Enr,params[0]);
			}
		}else if(S1Npe==0)
		h0c->Fill(Enr,params[0]);
    }

	h_test->Scale( year/Nsim/h_test->GetBinWidth(1));
	h0a->Scale(1. /Nsim/h0a->GetXaxis()->GetBinWidth(1)/h0a->GetYaxis()->GetBinWidth(1));
	h0b->Scale(1. /Nsim/h0b->GetXaxis()->GetBinWidth(1)/h0b->GetYaxis()->GetBinWidth(1));
	h0c->Scale(1. /Nsim/h0c->GetBinWidth(1));
	h0d->Scale(1. /Nsim/h0d->GetBinWidth(1));

	for(Long64_t i = 0; i < Nsim ; i++){//ER event by event
        if(i%1000000==0)std::cout<<i<<std::endl;
        params[0] = 1.0;
		double Eer = signal.GetER(params);// in keV
		YieldResult yields;
        QuantaResult quanta;
		yields = signal.GetYields_ER(Eer);
		quanta = signal.GetQuanta(yields);
		int S1Npe = signal.GetS1(quanta);
		int S2Npe = signal.GetS2(quanta);
		if(S2Npe!=0){
			h1b->Fill(S1Npe,log10(S2Npe),params[0]);
			if(S1Npe!=0){
				h1a->Fill(S1Npe,log10(S2Npe/S1Npe),params[0]);
				h1d->Fill(Eer,params[0]);
			}
		}else if(S1Npe==0)
		h1c->Fill(Eer,params[0]);
		h1->Fill(Eer,params[0]);
    }

	h1a->Scale(1. /Nsim/h1a->GetXaxis()->GetBinWidth(1)/h1a->GetYaxis()->GetBinWidth(1));
	h1b->Scale(1. /Nsim/h1b->GetXaxis()->GetBinWidth(1)/h1b->GetYaxis()->GetBinWidth(1));
	h1c->Scale(1. /Nsim/h1c->GetBinWidth(1));
	h1d->Scale(1. /Nsim/h1d->GetBinWidth(1));

	for(Long64_t i = 0; i < Nsim ; i++){//NR event by event
        if(i%1000000==0)std::cout<<i<<std::endl;
        params[0] = 1.0;
		double Enr = signal.GetNR(params);// in keV
		YieldResult yields;
        QuantaResult quanta;
		yields = signal.GetYields_ER(Enr);
		quanta = signal.GetQuanta(yields);
		int S1Npe = signal.GetS1(quanta);
		int S2Npe = signal.GetS2(quanta);
		if(S2Npe!=0){
			h2b->Fill(S1Npe,log10(S2Npe),params[0]);
			if(S1Npe!=0){
				h2a->Fill(S1Npe,log10(S2Npe/S1Npe),params[0]);
				h2d->Fill(Enr,params[0]);
			}
		}else if(S1Npe==0)
		h2c->Fill(Enr,params[0]);
		h0->Fill(Enr,params[0]);
    }

	h2a->Scale(1. /Nsim/h2a->GetXaxis()->GetBinWidth(1)/h2a->GetYaxis()->GetBinWidth(1));
	h2b->Scale(1. /Nsim/h2b->GetXaxis()->GetBinWidth(1)/h2b->GetYaxis()->GetBinWidth(1));
	h2c->Scale(1. /Nsim/h2c->GetBinWidth(1));
	h2d->Scale(1. /Nsim/h2d->GetBinWidth(1));

    fs->Write();

	return 0;
}
