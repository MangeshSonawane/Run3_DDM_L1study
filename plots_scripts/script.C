#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TPaveText.h"
#include <string>

using namespace ROOT;

void script(){

        int dxy = 1;
        int contours = 9;
        int scenario = 1;
        int relative = 1;
        int sample   = 3;

        double scalematrix[3][3] = {  0.528, 0.411, 0.361,
                                      0.767, 0.595, 0.360,
                                      0.744, 0.497, 0.168};

        double scale = 1.;

        if (relative) {scale = scalematrix[scenario-1][sample-1];}

        TString str[] = {"125_12_900", "125_25_1500", "125_50_3000"};

        TString sc_str = to_string(scenario);
        TFile *file_0 = TFile::Open("signal_"+str[sample-1]+"/DisplacedMuons_signal_"+str[sample-1]+"_sc"+sc_str+"_OR_BR2_40k.root");
        TH2F *h1;
        if (dxy==0){ h1 = (TH2F*)file_0->Get("h_eff_l1pt1_vs_l1pt2_dxy00");}
        else if (dxy==1) {h1 = (TH2F*)file_0->Get("h_eff_l1pt1_vs_l1pt2_dxy10");}
        h1->GetXaxis()->SetTitle("Lead Muon Pt Threshold [GeV]");
        h1->GetXaxis()->SetTitle("Sublead Muon Pt Threshold [GeV]");

        if (relative) {

          h1->Scale(1./scale);
          h1->SetMaximum(4);
          h1->SetMinimum(1);
          Double_t effcontours[] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0};
          h1->SetContour(15, effcontours);
        }

        else {
          h1->SetMaximum(1);
          h1->SetMinimum(0);
          h1->SetContour(10);
        }

        TCanvas *c1 = new TCanvas("c1", "c1", 1000, 800);
        c1->SetGridx();
        c1->SetGridy();
        
        TPaveText* pt00 = new TPaveText(1, 18.2, 17.8, 25);
        if (relative) pt00->AddText("Eff UPT OR BR2 relative to BR2: Scenario "+sc_str);
        else pt00->AddText("Eff UPT OR BR2 : Scenario "+sc_str);
        pt00->AddText("Sample : HTo2XTo4mu_"+str[sample-1]);
        pt00->AddText("Den : Gen Dimuons in acceptance");
        if (scenario == 2) pt00->AddText("Gen Pt >= 23 GeV for both muons");
        if (scenario == 3) pt00->AddText("Gen Pt >= 23 GeV, Gen Lxy >= 60 cm for both muons");
        pt00->AddText("Num : Den, L1 Pt >= Pt Threshold");
        if (dxy==1) {pt00->AddText("L1 Dxy (lead, sublead) >= (1,0)");}
        pt00->AddText("OR BR2");

        h1->Draw("colz");
        pt00->Draw();

        TPaveText* pt[10];
        for (int i = 0; i < contours; ++i){
        	pt[i] = new TPaveText(10, 10, 12, 11);
          TString tag = to_string(i+1);
          if (relative) pt[i]->AddText("1."+tag);
          else pt[i]->AddText("0."+tag);
        	pt[i]->Draw();
        } 

        TGraph *gr1;
        TGraph *gr0p8;
        TGraph *gr0p5;

        if (dxy==0){

          // Dxy (0,0)

          int x00_1[] = {8,8,10,10,13,13,15,15,18,18,23,23,25};
          int y00_1[] = {9,8,8,7,7,6,6,5,5,4,4,3,3};
          gr1 = new TGraph((sizeof(x00_1)/sizeof(x00_1[0])), x00_1, y00_1);

          int x00_0p8[] = {9,9,12,12,14,14,17,17,19,19,25};
          int y00_0p8[] = {10,8,8,7,7,6,6,5,5,4,4};
          gr0p8 = new TGraph((sizeof(x00_0p8)/sizeof(x00_0p8[0])), x00_0p8, y00_0p8);

          int x00_0p5[] = {10,10,12,12,13,13,15,15,18,18,21,21,25};
          int y00_0p5[] = {11,10,10,9,9,8,8,7,7,6,6,5,5};
          gr0p5 = new TGraph((sizeof(x00_0p5)/sizeof(x00_0p5[0])), x00_0p5, y00_0p5);
        }

        else if (dxy==1){

        // // Dxy (1,0)

          int x00_1[] = {5,6,6,10,10,13,13,15,15,25};
          int y00_1[] = {5,5,4,4,3,3,1,1,0,0};
          gr1 = new TGraph((sizeof(x00_1)/sizeof(x00_1[0])), x00_1, y00_1);

          int x00_0p8[] = {5,5,7,7,11,11,17,17,19,19,25};
          int y00_0p8[] = {6,5,5,4,4,3,3,1,1,0,0};
          gr0p8 = new TGraph((sizeof(x00_0p8)/sizeof(x00_0p8[0])), x00_0p8, y00_0p8);

          int x00_0p5[] = {6,7,7,11,11,15,15,25};
          int y00_0p5[] = {6,6,5,5,4,4,3,3};
          gr0p5 = new TGraph((sizeof(x00_0p5)/sizeof(x00_0p5[0])), x00_0p5, y00_0p5);
        }

        // // Dxy (1,1)

        // int x11[] = {5, 6};
        // int y11[] = {3, 0};
        // TGraph *gr = new TGraph(2, x10, y10);

        // // Dxy (2)

        // int x20[] = {4, 7, 8};
        // int y20[] = {4, 3, 0};
        // TGraph *gr = new TGraph(3, x10, y10);

        
        gr1->SetLineWidth(2);
        gr1->Draw("L same");

        gr0p8->SetLineWidth(2);
        gr0p8->SetLineColor(kRed);

        gr0p8->Draw("L same");

        gr0p5->SetLineWidth(2);
        gr0p5->SetLineColor(kGreen);

        gr0p5->Draw("L same");

        TLegend*leg = new TLegend(0.13, 0.43, 0.33, 0.63);
        leg->SetHeader("Rate contours");
        leg->AddEntry(gr1, "1 kHz");
        leg->AddEntry(gr0p8, "0.8 kHz");
        leg->AddEntry(gr0p5, "0.5 kHz");
        leg->SetLineWidth(0);
        leg->SetFillStyle(0);
        leg->Draw("same");

}