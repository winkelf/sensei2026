/*

This is a simple code that plots SENSEI images. It takes as inputs rootfiles named something like:

  hits_skp_sensei_2022-09-06_135K_run5_commissioning_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR10800_1_16.root

*/

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <random>
#include <cmath>
#include "TCanvas.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TBox.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "Math/ProbFunc.h"
#include <fitsio.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <TH2F.h>
#include <TCanvas.h>
using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////// Main function /////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int hduPlotter() {

  vector<int> maskList = {
                          //1,     // Neighbor
                          //4,     // Bleed
                          //8,     // Halo
                          //16,    // Crosstalk
                          //64,    // Edge
                          //128,   // Serial register
                          //512,   // Hot pixel
                          //1024,  // Hot column
                          //4096,  // Extended bleed
                          //8192,  // LEC
                          //8194,  // LEC + event?
                          //16384, // Full well mask
                          //32768, // Cluster shape
                         };

  gROOT->SetBatch(kTRUE);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetCanvasColor(kBlack);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetPadColor(kBlack);
  gStyle->SetStatColor(kBlack);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(kBlack);
  gStyle->SetTitleTextColor(kWhite);
  gStyle->SetLabelColor(kWhite, "XYZ");
  gStyle->SetTitleColor(kWhite, "XYZ");
  gStyle->SetOptStat(0);
  gStyle->SetHistFillColor(kBlack);
  gStyle->SetHistLineColor(kWhite);
  gStyle->SetFuncColor(kWhite);
  gStyle->SetFrameLineColor(kWhite);

  ////////////////////////////////////////////////////////////////////
  ////////////////////// Build grid //////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  TH2D* h  = new TH2D();
  TH2D* mh = new TH2D();

  vector<vector<int>> matrix_data(600, vector<int>(3500, 0));

  cout<<"Uploading SENSEI rootfile" <<endl;

  const char* fileName = "../inputs/hits_corr_proc_hadded.root";
  //const char* fileName = "../inputs/skim_good_exp216000.root";
  const char* treeName = "hitSumm";

  // Open the ROOT file
  TFile* file = TFile::Open(fileName, "READ");

  // Get the TTree
  TTree* tree = static_cast<TTree*>(file->Get(treeName));

  // Set up the branches
  int nSavedPix = 0;
  const int maxPix = 100000;  // Assume a reasonable upper limit for the array size

  float nElectrons;
  float xBary;
  float yBary;
  int   xPix[maxPix];
  int   yPix[maxPix];
  float ePix[maxPix];
  int   ohdu[maxPix];
  int   lta[maxPix];
  int   flag[maxPix];

  tree->SetBranchAddress("nSavedPix", &nSavedPix  );
  tree->SetBranchAddress("e",         &nElectrons );
  tree->SetBranchAddress("xBary",     &xBary      );
  tree->SetBranchAddress("yBary",     &yBary      );
  tree->SetBranchAddress("xPix",      xPix        );
  tree->SetBranchAddress("yPix",      yPix        );
  tree->SetBranchAddress("ePix",      ePix        );
  tree->SetBranchAddress("ohdu",      &ohdu       );
  tree->SetBranchAddress("flag",      &flag       );
  tree->SetBranchAddress("LTANAME",   lta         );

  Long64_t nEntries = tree->GetEntries();
  
  std::set<std::pair<int,int>> pairs;
  
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i);
    for (int j = 0; j < nSavedPix; j++)  pairs.insert({lta[j], ohdu[j]});
  }

  for (auto const& p : pairs) {
  
    int LTA  = 11;
    //int LTA  = p.first;
    int OHDU = p.second;
  
    std::cout << "Processing LTA=" << LTA  << " OHDU=" << OHDU << std::endl;
  
    TH2D h (Form("h_LTA%d_OHDU%d",LTA,OHDU), Form("LTA=%d, OHDU=%d",LTA,OHDU), 3200,0,3200,520,0,520 );
    TH2D mh(Form("mh_LTA%d_OHDU%d",LTA,OHDU),  "",                             3200,0,3200,520,0,520 );
  
    // refill histograms
    for (Long64_t i = 0; i < nEntries; i++) {
  
      tree->GetEntry(i);
  
      // keep only this pair
      if (lta[0]  != LTA)  continue;
      if (ohdu[0] != OHDU) continue;
  
      int mask = flag[0];
  
      for (int j = 0; j < nSavedPix; j++) {
  
        int x = xPix[j];
        int y = yPix[j];
        double e = ePix[j];
  
        if (x < 0 || x >= 3200 || y < 0 || y >= 520)
            continue;
  
        for (int m=0; m<maskList.size(); ++m) {
          if ((mask & maskList[m]) == maskList[m]) mh.SetBinContent(x+1,y+1,1);
        }
  
        h.SetBinContent(x+1,y+1,(e<0.68)?0:e);
      }
    }
  
    TCanvas c("c","c",2000,600);
  
    h.Draw("COLZ");
    h.GetYaxis()->SetAxisColor(kWhite);
    h.GetXaxis()->SetAxisColor(kWhite);

    const Int_t nLevels = 6;
    Double_t levels[nLevels] = {0.68, 1.5, 2.5, 3.5, 4.5, 5}; // Set levels at desired z-values
    
    // Apply the levels to the histogram
    h.SetContour(nLevels, levels);
    h.SetStats(0);
    h.GetYaxis()->SetAxisColor(kWhite);
    h.GetXaxis()->SetAxisColor(kWhite);
    h.GetZaxis()->SetAxisColor(kWhite);
    h.SetMinimum(0); // Set the minimum z value
    h.SetMaximum(5); // Set the maximum z value
    // Set the number of divisions and ticks for the z-axis
    h.GetZaxis()->SetNdivisions(14);   // Set number of divisions for the ticks
    h.GetZaxis()->SetTickLength(0.03); // Set the length of the ticks

    // Add custom labels near the color bar
    TLatex *text = new TLatex();
    text->SetTextSize(0.03); // Adjust text size
    text->SetTextAlign(12); // Align left, center vertically
    text->SetTextColor(kWhite);
  
    // Set label positions manually near the color bar (adjust for your canvas layout)
    text->DrawLatex(3250, 30,  " 0 e^{-}");
    text->DrawLatex(3250, 120, " 1 e^{-}"); 
    text->DrawLatex(3250, 230, " 2 e^{-}"); 
    text->DrawLatex(3250, 340, " 3 e^{-}"); 
    text->DrawLatex(3250, 440, " 4 e^{-}"); 
    text->DrawLatex(3250, 500, "+5 e^{-}");

    mh.Draw("BOX SAME");
  
    c.SaveAs(Form("pngs/LTA_%02d_OHDU_%02d.png",LTA,OHDU));
  }  


//  int LTA  = 18; 
//  int OHDU = 2;
//
//  set<int> ltaValues;
//
//  Long64_t nEntries = tree->GetEntries();
//  
//  for (Long64_t i = 0; i < nEntries; i++) {
//    tree->GetEntry(i);
//  
//    for (int j = 0; j < nSavedPix; j++) ltaValues.insert(lta[j]);
//  }
//  
//  std::cout << "Possible LTA values:\n";
//  for (auto v : ltaValues)  std::cout << v << std::endl;
//
//  exit(1);
//
//  // Create the histogram
//  *h = TH2D("", Form("skim_good_exp0.root    -     LTA: %1.i, OHDU: %1.i", LTA, OHDU), 3200, 0, 3200, 520, 0, 520);
//  h->SetStats(0);
//  h->GetYaxis()->SetAxisColor(kWhite);
//  h->GetXaxis()->SetAxisColor(kWhite);
//  
//  *mh = TH2D("", "skim_good_exp0.root", 3200, 0, 3200, 520, 0, 520);
//  mh->SetStats(0);
//  mh->SetFillColor(kRed);
//
//  int counter = 0;
//  float countE    = 0;
//  float countTot  = 0;
//
//  cout<< nEntries <<endl;
//
//  // Main loop (filling TH2)
//  for (Long64_t i = 0; i < nEntries; i++) {
//    tree->GetEntry(i);
//    counter+=1;
//
//    if ( (*lta)  != LTA ) continue; 
//    if ( (*ohdu) != OHDU) continue;
//
//    int mask = *flag;
//
//    for (int j = 0; j < nSavedPix; j++) {
//      int x    = xPix[j];
//      int y    = yPix[j];
//      double e = ePix[j];
//
//      // Ensure we don't go out of bounds
//      if (x >= 0 && x < 3500 && y >= 0 && y < 600) {
//
//        for (int m=0; m<maskList.size(); m++) {
//          if ( (mask&maskList.at(m))==maskList.at(m) ){
//            mh->SetBinContent(x+1, y+1, 1);
//	  } 
//	}
//
//        // 0.68 is the minimum charge to be considered a hit      
//        if (e<0.68){
//          h->SetBinContent(x+1, y+1, 0);
//          countE++;
//        }
//        else{ h->SetBinContent(x+1, y+1, e); }
//
//        countTot++;
//      } 
//    }
//  }// Main loop
//
//  ////////////////////////////////////////////////////////////////////
//  //////////////////// Plotting grid /////////////////////////////////
//  ////////////////////////////////////////////////////////////////////
//
//  TCanvas* c0;
//  TH2D* h2    = new TH2D();
//
//  c0  = new TCanvas("c0", "Grid", 2000, 600);	  
//  *h2 = *h;	  
//  c0->SetFillColor(kBlack);
//  c0->SetBatch(kTRUE);
//
//  const Int_t nLevels = 6;
//  Double_t levels[nLevels] = {0.68, 1.5, 2.5, 3.5, 4.5, 5}; // Set levels at desired z-values
//  
//  // Apply the levels to the histogram
//  h2->SetContour(nLevels, levels);
//  h2->SetStats(0);
//  h2->GetYaxis()->SetAxisColor(kWhite);
//  h2->GetXaxis()->SetAxisColor(kWhite);
//  h2->GetZaxis()->SetAxisColor(kWhite);
//  h2->SetMinimum(0); // Set the minimum z value
//  h2->SetMaximum(5); // Set the maximum z value
//
//  // Set the number of divisions and ticks for the z-axis
//  h2->GetZaxis()->SetNdivisions(14);   // Set number of divisions for the ticks
//  h2->GetZaxis()->SetTickLength(0.03); // Set the length of the ticks
//
//  //TH2D* hx = (TH2D*)h2->Clone();
//  h2->Draw("COLZ");
//  mh->SetLineColor(kRed);
//  mh->SetFillColor(kRed);
//  mh->Draw("BOX SAME");
//
//  // Add custom labels near the color bar
//  TLatex *text = new TLatex();
//  text->SetTextSize(0.03); // Adjust text size
//  text->SetTextAlign(12); // Align left, center vertically
//  text->SetTextColor(kWhite);
//
//  // Set label positions manually near the color bar (adjust for your canvas layout)
//  text->DrawLatex(3250, 30,  " 0 e^{-}");
//  text->DrawLatex(3250, 120, " 1 e^{-}"); 
//  text->DrawLatex(3250, 230, " 2 e^{-}"); 
//  text->DrawLatex(3250, 340, " 3 e^{-}"); 
//  text->DrawLatex(3250, 440, " 4 e^{-}"); 
//  text->DrawLatex(3250, 500, "+5 e^{-}"); 
//
//  c0->SaveAs("hdu.png");

  return 0;
}

int main(){
  hduPlotter();      	
  return 0;
}
