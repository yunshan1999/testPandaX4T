double interp1d(double x, int n, double xx[], double yy[]);

void s2_plotter()
{   

    int Nsim = 10000000;
    double pandaXScaleFactor = (200. - 20.)/(4.245e-1 - 2.647e-2);
    double x[13],r_ER[13],r_CEvNS[13],r_cathode[13];
    double s2[5] = {90.,200.,500.,1000.,3000.};
    double keVee[5] = {0.1452,0.2306,0.4791,0.8509,3.6999};

    TFile * f1 = new TFile("../build/testPandaX4T.root","read");
    TFile * f2 = new TFile("dataSpectrum.root","read");

    TH1F * h0 = (TH1F*)f1->Get("8B_Event_keVee");
    TH1F * h1 = (TH1F*)f1->Get("ER_Event_keVee");
    TH1F * h2 = (TH1F*)f1->Get("NR_Event_keVee");
    TH1F * h3 = (TH1F*)f2->Get("Data");
    TH1F * h4 = new TH1F("data","",h3->GetNbinsX(),h3->GetXaxis()->GetXmin()/pandaXScaleFactor,h3->GetXaxis()->GetXmax()/pandaXScaleFactor);

    TH1F * ha0 = new TH1F("xenon_ER","",100,0.,5.);
    TH1F * ha1 = new TH1F("xenon_CEvNS","",100,0.,5.);
    TH1F * ha2 = new TH1F("xenon_cathode","",100,0.,5.);

    //fill pandax's data hist
    for(int i = 1 ; i <= h4->GetNbinsX() ; i++){
        h4->SetBinContent(i, h3->GetBinContent(i) * pandaXScaleFactor);
    }

    //get xenon1t's results in to a list
    double temp0, temp1, temp2, temp3;
    ifstream infile("./XENON1TS2only.dat");
    for(int i = 0 ; i < 14; i++){
        infile>>temp0>>temp1>>temp2>>temp3;
        x[i] = temp0;
        r_ER[i] = temp1;
        r_CEvNS[i] = temp2;
        r_cathode[i] = temp3;
    }
    infile.close();

    //get hists

    for(int i = 1 ;i <= 100 ; i++){
        double energy = i * 5./100;
        double PE = interp1d(energy,5, keVee, s2);
        ha0->SetBinContent(i,interp1d(PE,13,x,r_ER));
        ha1->SetBinContent(i,interp1d(PE,13,x,r_CEvNS));
        ha2->SetBinContent(i,interp1d(PE,13,x,r_cathode));
        std::cout<<interp1d(PE,13,x,r_cathode)<<std::endl;
    }

    TCanvas * c1 = new TCanvas("c1");
    c1->SetLogy();
    c1->SetLogx();

    h0->SetLineColorAlpha(kBlack,0.9);
    h1->SetLineColorAlpha(kRed,0.9);
    h2->SetLineColorAlpha(kBlue,0.9);
    h4->SetLineColorAlpha(kMagenta,0.9);

    ha0->SetLineColorAlpha(kOrange-2,0.9);
    ha1->SetLineColorAlpha(kGreen,0.9);
    ha2->SetLineColorAlpha(kCyan,0.9);

    //std::cout<<h0->Integral()<<" "<<h1->Integral()<<" "<<h2->Integral()<<std::endl;

    h0->Scale(1023.93 / Nsim /h0->GetBinWidth(1)/365. );
    h1->Scale(2045.76 / Nsim /h1->GetBinWidth(1)/365. );
    h2->Scale(0.789236 / Nsim/h2->GetBinWidth(1)/365. );
    h4->Scale(1./h3->GetBinWidth(1));

    vector<TH1F*> vth1f = {h0, h1, h2, h4,ha0,ha1,ha2};
    for (const auto h: vth1f) h->SetLineWidth(2);
    
    h0->SetStats(0);
    h0->GetXaxis()->SetRangeUser(0.05,5.);
    h0->GetYaxis()->SetRangeUser(1e-5,2e4);
    h0->SetXTitle("Recoil Energy [keV]");
    h0->SetYTitle("Signal Rates(counts/day/ton/keV)");
    h0->SetTitle("S2-only");

    h0->Draw();
    h1->Draw("same");
    h2->Draw("same"); 
    h4->Draw("same");

    ha0->Draw("same");
    ha1->Draw("same");
    ha2->Draw("same");
    
    /*double topmax = 3000. / converter;

   //draw an axis on the top side
    TGaxis * axis = new TGaxis(10.,16.494819, 3000.,16.494819,0,topmax,510,"-"); 
    axis->SetLineColor(kBlack);
    axis->SetLabelColor(kBlack);
    axis->Draw();*/

    TLegend * leg = new TLegend(0.7,0.55,0.89,0.89);
    leg->AddEntry(h0,"Boron-8");
    leg->AddEntry(h1,"ER");
    leg->AddEntry(h2,"NR");
    leg->AddEntry(h4,"PandaX-II Data");

    leg->AddEntry(ha0,"XENON1T ER");
    leg->AddEntry(ha1,"XENON1T Boron-8");
    leg->AddEntry(ha2,"XENON1T cathode");
    leg->SetLineColor(kWhite);
    leg->Draw("same");
    
    //c1->SaveAs("S2-only rates.pdf");
    std::cout<<h0->Integral(1,6, "width")<<std::endl;
    std::cout<<h1->Integral(1,6,"width")<<std::endl;
   
}

double interp1d(double x,int n, double xx[], double yy[]){
    double temp;
    int i;
    for(i = 0 ; i < n  ; i++){
        if(x > xx[i]);
        else break;
    }
    if(i == 0) return yy[0];
    if(i == n) return yy[n-1];
    else{
        double x1 = xx[i-1];
        double y1 = yy[i-1];
        double x2 = xx[i];
        double y2 = yy[i];
        return (y2 - y1) * (x - x1) / (x2 - x1) + y1;
    }
}