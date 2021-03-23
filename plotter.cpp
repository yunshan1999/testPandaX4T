#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"

int plotter(){
    TFile * fs = new TFile("/Users/yunshan/work&study/sdu/neutrino/code3/code3/build/testPandaX4T.root", "r");
	TH1D * h0 = (TH1D *) fs->Get("8B_Event_Rate");
    TH1D * h1 = (TH1D *) fs->Get("8B_without_signal_energy");
    TH1D * h2 = (TH1D *) fs->Get("8B_with_signal_energy");
    h0->SetLineColor(1);
	h0->SetLineWidth(2);
	h1->SetLineColor(2);
	h2->SetLineColor(3);

    TCanvas * c0 = new TCanvas();
    c0->cd();
	c0->SetLogy();

    h0->Draw("HIST");
    h1->DrawClone("HIST SAME");
    h2->DrawClone("HIST SAME");

    auto legend = new TLegend(0.6,0.7,0.9,0.9);
    legend->AddEntry(h1,"without signal","l");
    legend->AddEntry(h2,"with signal","l");
    legend->Draw();

    return 0;

}
	