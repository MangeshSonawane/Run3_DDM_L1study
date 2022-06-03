#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TPaveText.h"
#include <string>

using namespace ROOT;

void DR_eta_script(){

        int k = 1;
        int contours[] = {0,0};
        TString Dxy = "";
        TString upt = "UPT=(6,4)";

        TFile *file_0 = TFile::Open("DisplacedMuons_NuGun_rate_added_rate_ER2_1M_DQ_6_4.root");
        TH2F *h1;
        if (k==0){ 
          h1 = (TH2F*)file_0->Get("h_rate_dR_vs_eta_dxy00");
          Dxy = "Dxy=(0,0)";
        }
        else if (k==1) {
          h1 = (TH2F*)file_0->Get("h_rate_dR_vs_eta_dxy10");
          Dxy = "Dxy=(1,0)";
        }

        h1->SetTitle("Added Rate wrt ER2, "+Dxy+", DQ, "+upt);
        // h1->GetXaxis()->SetTitle("Lead Muon Pt Threshold [GeV]");
        // h1->GetXaxis()->SetTitle("Sublead Muon Pt Threshold [GeV]");

        Double_t contour[] = {0, 0.3, 0.5, 0.8, 1, 2, 3, 4};
        h1->SetMaximum(10);
        h1->SetMinimum(0);
        h1->SetContour(8, contour);

        TCanvas *c1 = new TCanvas("c1", "c1", 1000, 800);
        c1->SetGridx();
        c1->SetGridy();
        
        TPaveText* pt00 = new TPaveText(1, 18.2, 17.8, 25);
        pt00->AddText("Neutrino Gun : 1M events");

        h1->Draw("colz");
        pt00->Draw();

        TPaveText* pt[10];
        for (int i = contours[0]; i < contours[1]; ++i){
        	pt[i] = new TPaveText(10, 10, 12, 11);
        	pt[i]->Draw();
        } 

        TGraph *gr1;
        TGraph *gr0p8;
        TGraph *gr0p5;

        if (k==0){

          // Dxy (0,0)

          Double_t x00_1[] = {0., 0.8, 0.8, 1.0, 1.0};
          Double_t y00_1[] = {1.2, 1.2, 0.8, 0.8, 0.};
          gr1 = new TGraph((sizeof(x00_1)/sizeof(x00_1[0])), x00_1, y00_1);

          Double_t x00_0p8[] = {0., 0.8, 0.8, 1.0, 1.0};
          Double_t y00_0p8[] = {1.0, 1.0, 0.6, 0.6, 0.0};
          gr0p8 = new TGraph((sizeof(x00_0p8)/sizeof(x00_0p8[0])), x00_0p8, y00_0p8);

          Double_t x00_0p5[] = {0.0, 0.8, 0.8};
          Double_t y00_0p5[] = {0.4, 0.4, 0.};
          gr0p5 = new TGraph((sizeof(x00_0p5)/sizeof(x00_0p5[0])), x00_0p5, y00_0p5);
        }

        else if (k==1){

        // Dxy (1,0)

          Double_t x00_1[] = {0.0,0.8,0.8,1.0,1.0,1.6,1.6,1.8,1.8};
          Double_t y00_1[] = {2.8,2.8,1.4,1.4,1.0,1.0,0.4,0.4,0.0};
          gr1 = new TGraph(sizeof(x00_1)/sizeof(x00_1[0]), x00_1, y00_1);

          Double_t x00_0p8[] = {0.0,0.8,0.8,1.0,1.0,1.6,1.6,1.8,1.8};        
          Double_t y00_0p8[] = {2.4,2.4,1.2,1.2,0.8,0.8,0.2,0.2,0.0};        
          gr0p8 = new TGraph(sizeof(x00_0p8)/sizeof(x00_0p8[0]), x00_0p8, y00_0p8);

          Double_t x00_0p5[] = {0.0,0.8,0.8,1.0,1.0,1.4,1.4,1.6,1.6};
          Double_t y00_0p5[] = {1.0,1.0,0.8,0.8,0.6,0.6,0.4,0.4,0.0};
          gr0p5 = new TGraph(sizeof(x00_0p5)/sizeof(x00_0p5[0]), x00_0p5, y00_0p5);
        }

        
        gr1->SetLineWidth(2);
        // gr1->Draw("L same");

        gr0p8->SetLineWidth(2);
        gr0p8->SetLineColor(kRed);

        // gr0p8->Draw("L same");

        gr0p5->SetLineWidth(2);
        gr0p5->SetLineColor(kGreen);

        // gr0p5->Draw("L same");

        TLegend*leg = new TLegend(0.13, 0.43, 0.33, 0.63);
        leg->SetHeader("Rate contours");
        leg->AddEntry(gr1, "1 kHz");
        leg->AddEntry(gr0p8, "0.8 kHz");
        leg->AddEntry(gr0p5, "0.5 kHz");
        leg->SetLineWidth(0);
        leg->SetFillStyle(0);
        // leg->Draw("same");

}