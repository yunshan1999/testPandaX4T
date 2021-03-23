#include <iostream>
#include "signal.h"
#include <TGraphErrors.h>
#include <fstream>
#include "TLine.h"

int main(){
	std::cout<<"code5 running..."<<std::endl;
    Long64_t Nsim = 1000000;
    TString filename = "testPandaX4T.root";
    TFile * fs = new TFile(filename.Data(), "RECREATE");

    TLorentzVector ki[2],kf[2];
    double params[2];//[0] is weight, [1] is Xe mass number  
	double year = 3600. * 24. * 365.;

	TH1D * h0 = new TH1D("8B_Event_Rate", "",4000, 0., 4.5); //draw 8B rate_vs_recoil energy
    TH2D * h0a = new TH2D("8B_Band","",800,0.,80.,500,0.,5.); 
	TH1D * h0b = new TH1D("8B_Event_keVee","",100, 0., 5.);
	TH1D * h0c = new TH1D("8B_S2","",300,0.,3000.);
	TH1D * h1 = new TH1D("ER_Event_Rate","",4000, 0., 10.);
	TH2D * h1a = new TH2D("ER_Band","",800,0.,80.,500,0.,5.);
	TH1D * h1b = new TH1D("ER_Event_keVee","",100, 0., 5.);
	TH1D * h1c = new TH1D("ER_S2","",300,0.,3000.);
	TH1D * h2 = new TH1D("NR_Event_Rate","",4000,0.,50.);
	TH2D * h2a = new TH2D("NR_Band","",800,0.,80.,500,0.,5.);
	TH1D * h2b = new TH1D("NR_Event_keVee","",100, 0., 5.);
	TH1D * h2c = new TH1D("NR_S2","",300,0.,3000.);

	h0->SetDirectory(fs);
	h0a->SetDirectory(fs);
	h0b->SetDirectory(fs);
	h0c->SetDirectory(fs);
	h1->SetDirectory(fs);
	h1a->SetDirectory(fs);
	h1b->SetDirectory(fs);
	h1c->SetDirectory(fs);
	h2->SetDirectory(fs);
	h2a->SetDirectory(fs);
	h2b->SetDirectory(fs);
	h2c->SetDirectory(fs);

	Signalcalc signal;

	//generate 8B rate 
	for(Long64_t i = 0; i < Nsim ; i++){
	    params[0] = 1.0;
		double En = gRandom->Uniform(0., 16.);//in MeV
		double Enr = signal.GetB8NR(En, params);// in keV
	    h0->Fill(Enr,params[0]);
	}
	h0->Scale(year/ Nsim);

	for(Long64_t i = 0; i < Nsim ; i++){
        params[0] = 1.0;
		double Eer = signal.GetER(params);// in keV
		h1->Fill(Eer,params[0]);
    }
	h1->Scale(365./ Nsim);

	for(Long64_t i = 0; i < Nsim ; i++){
	    params[0] = 1.0;
		double Enr = signal.GetNR(params);// in keV
		h2->Fill(Enr,params[0]);
	}		
	h2->Scale(365./ Nsim);

	//do event-by-event simulation
	double Nsim_8B = double(h0->Integral());
	double Nsim_er = double(h1->Integral());
	double Nsim_nr = double(h2->Integral());
	std::cout<<Nsim_8B<<" "<<Nsim_er<<" "<<Nsim_nr<<std::endl;

	for(int i = 0; i < Nsim  ; i++){
		if(i%Nsim==0)std::cout<<i<<std::endl;
		double Enr = h0->GetRandom();
	    QuantaResult quanta;
		quanta = signal.GetQuanta_NR(Enr);
		vector<double> S1 = signal.GetS1(quanta);
		vector<double> S2 = signal.GetS2(quanta);
		if(S1[6]!=0.&&S2[6]!=0.){
			h0a->Fill(S1[6],log10(S2[6]/S1[6]));
		}
		if(S2[5]!=0.&&S1[6]==0.){
			h0b->Fill(Enr * quanta.lindhard);
		}
			//h0c->Fill(S2[5]);
	}

	std::ofstream outfile;
    outfile.open("./cpu_outputs.dat");
	for(int i = 0; i < Nsim  ; i++){
		double Eer = h1->GetRandom();
	    QuantaResult quanta;
		quanta = signal.GetQuanta_ER(Eer);
		vector<double> S1 = signal.GetS1(quanta);
		vector<double> S2 = signal.GetS2(quanta);
		if (outfile.is_open()) 
        {
            outfile <<i<<" "<< S1[6] << " " << S2[6] << std::endl;
	    }
	    else{
			std::cout<<"error in file opening"<<std::endl;
		}

		if(S1[6]!=0.&&S2[6]!=0.){
			h1a->Fill(S1[6],log10(S2[6]/S1[6]));
		}
		//double Er = W_EV * (S1[6]/ 1.2/ 0.09997 + S2[6]/ 1.2/ 28./ 0.727);
		//h1b->Fill(Er);
		if(S2[5]!=0.&&S1[6]==0.){
			//h1b->Fill(S2[5]);
			h1b->Fill(Eer * quanta.lindhard);
		}

	}
	outfile.close();

	for(int i = 0; i < Nsim  ; i++){
		double Enr = h2->GetRandom();
	    QuantaResult quanta;
		quanta = signal.GetQuanta_NR(Enr);
		vector<double> S1 = signal.GetS1(quanta);
		vector<double> S2 = signal.GetS2(quanta);
		if(S1[6]!=0.&&S2[6]!=0.){
			h2a->Fill(S1[6],log10(S2[6]/S1[6]));
		}
		if(S2[5]!=0.&&S1[6]==0.){
			h2b->Fill(Enr * quanta.lindhard);
		}
		//double Er = W_EV * (S1[6]/ 1.2/ 0.09997 + S2[6]/ 1.2/ 28.);
		//h2e->Fill(Er);
		//std::cout<<Enr<<" "<<Er<<std::endl;
	}	

	TH1D * hfit = h1a->ProjectionY("test",101,110);
	TF1 *f1 = new TF1("f1","gaus",1.5,4.);
	f1->SetParameters(hfit->GetMaximum(), hfit->GetMean(), hfit->GetRMS() ); 
	hfit->Fit("f1");

	/*TCanvas * c1 = new TCanvas("c1");
    //gStyle->SetOptStat(0);

    double contours[8];
	contours[0] = 2;
	contours[1] = 3;
	contours[2] = 5;
	contours[3] = 7;
	contours[4] = 10;
    contours[5] = 15;
	contours[6] = 20;
	contours[7] = 25;

    h0a->SetContour(8, contours);
	h0a->Draw("CONT Z LIST");
	c1->Print( "/Users/yunshan/work&study/sdu/neutrino/code5/plots/cont.pdf)","pdf");
    //c1->Update();

   TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
   TList* contLevel = NULL;
   TGraph* curv     = NULL;
   TGraph* gc       = NULL;
 
   Int_t nGraphs    = 0;
   Int_t TotalConts = 0;
 
   if (conts == NULL){
      printf("*** No Contours Were Extracted!\n");
      TotalConts = 0;
      return 0;
   } else {
      TotalConts = conts->GetSize();
   }
 
   printf("TotalConts = %d\n", TotalConts);
 
   for(int i = 0; i < TotalConts; i++){
      contLevel = (TList*)conts->At(i);
      printf("Contour %d has %d Graphs\n", i, contLevel->GetSize());
      nGraphs += contLevel->GetSize();
   }*/
	
		//get median of NR
	double Ex[770], Ey[770];
	for(int i = 0; i < 770; i++){
        int firstbin = i + 31;
        int lastbin = i + 31;
        double aprob[5] = {0.001,0.1585,0.5,0.8415,0.999};
        double median[5];
        TH1D * h_temp = h2a->ProjectionY("proj",firstbin,lastbin);
		if(h_temp->Integral()==0.)continue;
        h_temp->GetQuantiles(5,median,aprob);
        Ex[i] = 0.1 * (firstbin - 0.5) ;
        Ey[i] = median[2];
        }


	//do ER rejection
	/*double perc[770];
	std::vector<double> xs1;
	std::vector<double> LR;
	std::vector<double> LR_err;
	std::vector<double> LR_errx;
	double p = 0.;
	for(int i = 0 ; i < 770 ; i++){
		int firstbin = i + 31;
		int lastbin = i + 31;
		TH1D * h_temp = h1a->ProjectionY("proj",firstbin,lastbin);
		if(h_temp->Integral()==0.)std::cout<<"error"<<std::endl;
		p +=  h_temp->Integral(1,h_temp->FindBin(Ey[i]));
		perc[i] = h_temp->Integral(1,h_temp->FindBin(Ey[i]))/ h_temp->Integral() * 100.;
		if(i<50)std::cout<<h_temp->Integral(1,h_temp->FindBin(Ey[i]))<<std::endl;
	}
	std::cout<<p/h1a->Integral(31,800,1,500) * 100.<<" "<<sqrt(p)/h1a->Integral(31,800,1,500) * 100.<<std::endl;
	std::cout<<"denominator:"<<p<<std::endl;
	std::cout<<"numerator:"<<h1a->Integral(31,800,1,500)<<std::endl;
	
	for(int i = 0; i < 77; i++){
			double temp = 0.;
			xs1.push_back(i+3);
			for(int j = 0; j < 10; j++){
				int firstbin = 10 * i + j + 31;
				int lastbin = 10 * i + j + 31;
				TH1D * h_temp = h1a->ProjectionY("proj",firstbin,lastbin);
				temp += h_temp->Integral() * perc[10 * i + j];
			}
			LR.push_back(temp/h1a->Integral(10 * i + 31, 10 * i + 31 + 9, 1, 500));
			LR_err.push_back(sqrt(temp)/h1a->Integral(10 * i + 31, 10 * i + 31 + 9, 1, 500));
			LR_errx.push_back(0.);
	}
	double * xs11 = xs1.data();
	double * LRR = LR.data();
	double * LRR_err = LR_err.data();
	double * LRR_errx = LR_errx.data();*/

	/*ofstream myfile ("../plots/example.txt");
  	if (myfile.is_open())
  	{
  	  for(int count = 0; count < xs1.size(); count ++){
  	      myfile << xs1[count] << " "<<LRR[count]<<" "<<LRR_err[count]<<"\n" ;
  	  }
  	  myfile.close();
  	}
  	else cout << "Unable to open file";*/

	
	//get 2-d quantiles of B8
	/*std::vector<double> s1;
	std::vector<double> ymin;
	std::vector<double> ymax;
	for(int i = 30; i < 50; i++){
        int firstbin = i + 1;
        int lastbin = i + 1;
        double aprob[5] = {0.05,0.1585,0.5,0.8415,0.95};
        double median[5];
        TH1D * h_temp = h0a->ProjectionY("proj",firstbin,lastbin);
		if(h_temp->Integral()<10.)continue;
        else{
			h_temp->GetQuantiles(5,median,aprob);
        	s1.push_back(0.1 * (firstbin - 0.5));
        	ymin.push_back(median[0]);
        	ymax.push_back(median[4]);
		}
    }

	double ss1[s1.size()],yymin[s1.size()],yymax[s1.size()];
	for(int i = 0; i < s1.size(); i++){
		ss1[i] = s1[i];
		yymin[i] = ymin[i];
		yymax[i] = ymax[i];
	}

	//calculate bkg event within B8 quantiles
	double Nb8 = 0., Ner = 0., Nnr = 0.;
	for(int i = 0; i < s1.size(); i++){
		int firstbin = i + 1;
		int lastbin = i + 1;
		TH1D * h_temp0 = h0a->ProjectionY("proj0",firstbin, lastbin );
		TH1D * h_temp1 = h1a->ProjectionY("proj1",firstbin, lastbin );
		TH1D * h_temp2 = h2a->ProjectionY("proj2",firstbin, lastbin );
		Nb8 += h_temp0->Integral(h_temp0->FindBin(yymin[i]),h_temp0->FindBin(yymax[i]));
		Ner += h_temp1->Integral(h_temp1->FindBin(yymin[i]),h_temp1->FindBin(yymax[i]));
		Nnr += h_temp2->Integral(h_temp2->FindBin(yymin[i]),h_temp2->FindBin(yymax[i]));
	}

	std::cout<<"NB8:"<<Nb8 * Nsim_8B / Nsim<<" "<<sqrt(Nb8)/Nsim * Nsim_8B<<endl;
	std::cout<<"Ner:"<<Ner * Nsim_er / Nsim<<" "<<sqrt(Ner)/Nsim * Nsim_er<<endl;
	std::cout<<"Nnr:"<<Nnr * Nsim_nr / Nsim<<" "<<sqrt(Nnr)/Nsim * Nsim_nr<<endl;*/


	//TGraph * g0 = new TGraph(s1.size(), ss1, yymin);
	//TGraph * g1 = new TGraph(s1.size(), ss1, yymax);
	//TGraph * g1 = new TGraph(770, Ex, Ey);
	//TGraph * g2 = new TGraph(500,Ex,perc);

	//TCanvas * c1 = new TCanvas();
	//c1->SetGrid();
	//c1->SetTitle("leakage ratio");
	//g1->GetXaxis()->SetTitle("S1[PE]");
	//g2->GetYaxis()->SetTitle("leakage ratio[%]");
	//TLine * l = new TLine(3.,p/h1a->Integral(31,500,1,500) * 100.,80.,p/h1a->Integral(31,500,1,500) * 100.);
	//l->SetLineColorAlpha(kBlue,0.7);
	//l->SetLineStyle(2);
	//l->SetLineWidth(2);
	/*TPad* p1 = new TPad("p1", "p1", 0.0, 0.0, 1., 0.3, 0); p1 ->Draw();
   	TPad* p2 = new TPad("p2", "p2", 0.0, 0.3, 1., 1., 0); p2 ->Draw();
	p2->cd();*/
	//h1a->SetStats(0);
	//h0a->GetXaxis()->SetRangeUser(0.,20.);
	//h0a->Draw("colz");
	//h1a->SetXTitle("cS1[PE]");
	//h1a->SetYTitle("log(cS2/cS1)");
	//h1a->Draw("colz");
	//g0->SetLineColor(2);
	//g0->SetLineStyle(1);
	//g0->Draw("SAME");
	//g1->SetLineColor(2);
	//g1->SetLineStyle(1);
	//g1->Draw("SAME");
	//g1->SetLineColor(1);
	//g1->SetLineStyle(2);
	//g1->Draw();
	//g1->Write();

	//TLine *line = new TLine(ss1[s1.size()-1],yymin[s1.size()-1],ss1[s1.size()-1],yymax[s1.size()-1]);
	//line->SetLineStyle(2);
	//line->SetLineColor(2);
	//line->Draw("SAME");
	//p1->cd();
	//g1->GetXaxis()->SetRangeUser(0,50);
	//g1->Draw();
	//g2->Draw();
	//l->Draw("-");
	//double size = sizeof(xs11)/sizeof(xs11[0]);

	//TGraphErrors * gr = new TGraphErrors(xs1.size(),xs11,LRR,LRR_errx,LRR_err); 
	//gr->SetMarkerColor(4); 
	//gr->SetMarkerStyle(20); 
	//gr->Draw("ACP");

	//l->Draw("same");
	//c1->Print( "/Users/yunshan/work&study/sdu/neutrino/code5/plots/leakage ratio.pdf)","pdf");

    fs->Write();
	return 0;
}
