#include "TFile.h"
#include "TString.h"
#include "TH2.h"
#include "TPaveText.h"
#include <string>

using namespace ROOT;

void script(){

        int k = 1;
        int contours = 0;

        TFile *file_0 = TFile::Open("DisplacedMuons_private_NuGun_rate_added_rate_ER2_1M.root");
        TH2F *h1;
        if (k==0){ h1 = (TH2F*)file_0->Get("h_rate_l1pt1_vs_l1pt2_dxy0");}
        else if (k==1) {h1 = (TH2F*)file_0->Get("h_rate_l1pt1_vs_l1pt2_dxy10");}
        h1->GetXaxis()->SetRangeUser(0, 25);
        h1->GetYaxis()->SetRangeUser(0, 25);
        h1->GetXaxis()->SetTitle("Lead Muon UPT threshold [GeV]");
        h1->GetYaxis()->SetTitle("Sublead Muon UPT threshold [GeV]");

        h1->SetMaximum(10);
        Double_t contour[] = {0, 0.3, 0.5, 0.8, 1, 2, 3, 4, 10};
        h1->SetContour(9, contour);

        TCanvas *c1 = new TCanvas("c1", "c1", 1000, 800);
        c1->SetGridx();
        c1->SetGridy();
        
        TPaveText* pt00 = new TPaveText(1, 22, 15, 25);
        pt00->AddText("Neutrino Gun : 1M events");
        pt00->AddText("Requiring |#eta_{#mu}| < 2.0");
//        pt00->AddText("With PU in [48,58] : 174346 events");

        h1->Draw("colz");
        pt00->Draw();

        TPaveText* pt[10];
        for (int i = 0; i < contours; ++i){
        	pt[i] = new TPaveText(10, 10, 11, 11);
          // pt[i]->AddText("0")
        	pt[i]->Draw();
        } 

        TGraph *gr1;
        TGraph *gr0p8;
        TGraph *gr0p5;

        if (k==0){

          // Dxy (0,0)

          //Without eta restriction
          // int x00_1[] = {8,8,10,10,13,13,15,15,18,18,23,23,25};
          // int y00_1[] = {9,8,8,7,7,6,6,5,5,4,4,3,3};

          // With eta restriction
          int x00_1[] = {6,8,8,11,11,12,12,16,16,22,22,25};
          int y00_1[] = {7,7,6,6,5,5,4,4,3,3,1,1};
          gr1 = new TGraph((sizeof(x00_1)/sizeof(x00_1[0])), x00_1, y00_1);

          //Without eta restriction
          // int x00_0p8[] = {9,9,12,12,14,14,17,17,19,19,25}; 
          // int y00_0p8[] = {10,8,8,7,7,6,6,5,5,4,4};

          //With eta restriction
          int x00_0p8[] = {7,7,9,9,12,12,14,14,18,18,25}; 
          int y00_0p8[] = {8,7,7,6,6,5,5,4,4,3,3};
          gr0p8 = new TGraph((sizeof(x00_0p8)/sizeof(x00_0p8[0])), x00_0p8, y00_0p8);

          //Without eta restriction
          //int x00_0p5[] = {10,10,12,12,13,13,15,15,18,18,21,21,25};
          //int y00_0p5[] = {11,10,10,9,9,8,8,7,7,6,6,5,5};

          //With eta restriction
          int x00_0p5[] = {8,8,12,12,14,14,17,17,23,23,25};
          int y00_0p5[] = {9,7,7,6,6,5,5,4,4,3,3};
          gr0p5 = new TGraph((sizeof(x00_0p5)/sizeof(x00_0p5[0])), x00_0p5, y00_0p5);
        }

        else if (k==1){

        // // Dxy (1,0)

          //int x00_1[] = {5,6,6,10,10,13,13,15,15,25};
          //int y00_1[] = {5,5,4,4,3,3,1,1,0,0};
          int x00_1[] = {4,6,6,7,7,8,8,25};
          int y00_1[] = {4,4,3,3,1,1,0,0};
          gr1 = new TGraph((sizeof(x00_1)/sizeof(x00_1[0])), x00_1, y00_1);

          //int x00_0p8[] = {5,5,7,7,11,11,17,17,19,19,25};
          //int y00_0p8[] = {6,5,5,4,4,3,3,1,1,0,0};
          int x00_0p8[] = {4,4,7,7,8,8,9,9,25};
          int y00_0p8[] = {5,4,4,3,3,1,1,0,0};
          gr0p8 = new TGraph((sizeof(x00_0p8)/sizeof(x00_0p8[0])), x00_0p8, y00_0p8);

          //int x00_0p5[] = {6,7,7,11,11,15,15,25};
          //int y00_0p5[] = {6,6,5,5,4,4,3,3};
          int x00_0p5[] = {5,5,7,7,9,9,10,10,12,12,25};
          int y00_0p5[] = {6,5,5,4,4,3,3,1,1,0,0};
          gr0p5 = new TGraph((sizeof(x00_0p5)/sizeof(x00_0p5[0])), x00_0p5, y00_0p5);

        }

        
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