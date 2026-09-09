#include <iostream>
#include <vector>
#include <set>
#include <map>
#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLatex.h"

using namespace std;

struct PixelInfo {
  short x;
  short y;
};

int hduPlotter() {

  // Select masks to display
  vector<int> maskList = {
                          //1,     // Neighbor
                          4,     // Bleed
                          8,     // Halo
                          16,    // Crosstalk
                          64,    // Edge
                          128,   // Serial register
                          512,   // Hot pixel
                          1024,  // Hot column
                          4096,  // Extended bleed
                          8192,  // LEC
                          16384, // Full well mask
                          //32768, // Cluster shape
		          //65536  // binomialMask
                         };

  gROOT ->SetBatch(kTRUE);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetCanvasColor(kBlack);
  gStyle->SetFrameFillColor(kBlack);
  gStyle->SetPadColor(kBlack);
  gStyle->SetStatColor(kBlack);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(kBlack);
  gStyle->SetTitleTextColor(kWhite);
  gStyle->SetLabelColor(kWhite,"XYZ");
  gStyle->SetTitleColor(kWhite,"XYZ");
  gStyle->SetHistFillColor(kBlack);
  gStyle->SetHistLineColor(kWhite);

  ////////////////////////////////////////////////////////////
  // Open file and get trees
  ////////////////////////////////////////////////////////////

  //const char* fileName = "../inputs/hits_corr_proc_hadded_EXPOSURE0_11.root"; //--------------------------------------------------------------------------------------------------------------------------
  //const char* fileName = "../inputs/hits_corr_proc_hadded_EXPOSURE216000_11.root";
  //const char* fileName = "/home/fwinkel/Desktop/SENSEI/SENSEI2023/commissioning/hits/hits_corr_proc_skp_sensei_2023-02-14_135K_run3_commissioning_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR10800_6_64.root";
  const char* fileName = "/home/fwinkel/Desktop/SENSEI/SENSEI2023/commissioning/hits/hits_corr_proc_skp_sensei_2023-02-14_135K_run1_commissioning_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_CLEAR10800_6_47.root";
  TFile* file          = TFile::Open(fileName,"READ");
  TTree* hitTree       = (TTree*)file->Get("hitSumm");
  TTree* calPixTree    = (TTree*)file->Get("calPixTree");

  std::string exposure(fileName);
  std::string expo;
  size_t pos = exposure.find("EXPOSURE");
  if (pos != std::string::npos) {
    size_t end = exposure.find('_', pos);
    expo = exposure.substr(pos, end - pos);
    std::cout << expo << std::endl;
  }

  ////////////////////////////////////////////////////////////
  // hitSumm branches
  ////////////////////////////////////////////////////////////

  const int maxPix = 100000;

  int nSavedPix;

  int runID;
  int lta;
  int ohdu;
  int flag;

  float e;
  float xBary;
  float yBary;

  int   xPix[maxPix];
  int   yPix[maxPix];
  float ePix[maxPix];

  hitTree->SetBranchAddress("runID",    &runID    );
  hitTree->SetBranchAddress("LTANAME",  &lta      );
  hitTree->SetBranchAddress("ohdu",     &ohdu     );
  hitTree->SetBranchAddress("flag",     &flag     );
  hitTree->SetBranchAddress("e",        &e        );
  hitTree->SetBranchAddress("xBary",    &xBary    );
  hitTree->SetBranchAddress("yBary",    &yBary    );
  hitTree->SetBranchAddress("nSavedPix",&nSavedPix);
  hitTree->SetBranchAddress("xPix",     xPix      );
  hitTree->SetBranchAddress("yPix",     yPix      );
  hitTree->SetBranchAddress("ePix",     ePix      );

  ////////////////////////////////////////////////////////////
  // calPixTree branches
  ////////////////////////////////////////////////////////////

  int cal_x;
  int cal_y;
  int cal_mask;
  int cal_lta;
  int cal_ohdu;

  calPixTree->SetBranchAddress("x",      &cal_x   );
  calPixTree->SetBranchAddress("y",      &cal_y   );
  calPixTree->SetBranchAddress("mask",   &cal_mask);
  calPixTree->SetBranchAddress("LTANAME",&cal_lta );
  calPixTree->SetBranchAddress("ohdu",   &cal_ohdu);

  // Discover available LTA/OHDU pairs
  set<pair<int,int>> pairs;

  Long64_t nEvents = hitTree->GetEntries();

  for (Long64_t i=0;i<nEvents;i++) {
    hitTree->GetEntry(i);
    pairs.insert({lta,ohdu});
  }

  // Build mask cache
  map< pair<int,int>,  vector<PixelInfo> > maskCache;

  Long64_t nPixels = calPixTree->GetEntries();

  cout << "Building mask cache from " << nPixels << " pixels" << endl;

  for (Long64_t i=0;i<nPixels;i++) {

    if (i%1000000==0) cout  << i << " / " << nPixels << endl;

    calPixTree->GetEntry(i);

    bool keepPixel = false;

    for (size_t m=0;m<maskList.size();m++) {

      if ((cal_mask & maskList[m]) == maskList[m]) {
        keepPixel = true;
        break;
      }
    }

    if (!keepPixel) continue;
    if (cal_x < 0 || cal_x >= 3200)  continue;
    if (cal_y < 0 || cal_y >= 520)   continue;

    maskCache[ {cal_lta,cal_ohdu} ].push_back( {(short)cal_x,(short)cal_y} );
  }

  cout << "Mask cache complete" << endl;

  // Produce plots
  for (auto const& p : pairs) {

    int LTA  = p.first;
    int OHDU = p.second;

    cout << "Processing LTA=" << LTA << " OHDU=" << OHDU << endl;

    TH2D h(  Form("h_LTA%d_OHDU%d %s",  LTA, OHDU, expo.c_str()), Form("LTA=%d OHDU=%d %s", LTA, OHDU, expo.c_str()), 3200,0,3200, 520,0,520 );
    TH2D mh( Form("mh_LTA%d_OHDU%d %s", LTA, OHDU, expo.c_str()), "",                                         3200,0,3200, 520,0,520 );

    // Fill charge image
    for (Long64_t i=0;i<nEvents;i++) {

      hitTree->GetEntry(i);
      if (lta  != LTA)  continue;
      if (ohdu != OHDU) continue;

      for (int j=0;j<nSavedPix;j++) {

        int x = xPix[j];
        int y = yPix[j];

        if (x < 0 || x >= 3200) continue;
        if (y < 0 || y >= 520)  continue;

        double charge = ePix[j];

        h.SetBinContent( x+1, y+1, (charge < 0.68) ? 0 : charge );
      }
    }

    // Fill mask overlay
    auto it = maskCache.find({LTA,OHDU});

    if (it != maskCache.end()) {
      for (const auto& pix : it->second) {
        mh.SetBinContent( pix.x+1, pix.y+1, 1 );
      }
    }

    // Draw
    TCanvas c( Form("c_%d_%d",LTA,OHDU), "", 2000, 600 );

    const Int_t nLevels      = 6;
    Double_t levels[nLevels] = { 0.68, 1.5, 2.5, 3.5, 4.5, 5.0 };

    h.SetContour( nLevels, levels );
    h.SetMinimum(0);
    h.SetMaximum(5);
    h.GetYaxis()->SetAxisColor(kWhite);
    h.GetXaxis()->SetAxisColor(kWhite);
    h.GetZaxis()->SetAxisColor(kWhite);
    h.SetStats(0);
    h.Draw("COLZ");

    mh.SetFillColor(kRed);
    mh.SetLineColor(kRed);
    mh.Draw("BOX SAME");

    TLatex text;
    text.SetTextColor(kWhite);
    text.SetTextSize(0.03);
    text.DrawLatex(3250,30 ," 0 e^{-}");
    text.DrawLatex(3250,120," 1 e^{-}");
    text.DrawLatex(3250,230," 2 e^{-}");
    text.DrawLatex(3250,340," 3 e^{-}");
    text.DrawLatex(3250,440," 4 e^{-}");
    text.DrawLatex(3250,500,"+5 e^{-}");

    c.SaveAs( Form( "pngs_2023/LTA_%02d_OHDU_%02d_%s.png", LTA, OHDU, expo.c_str() ) );
  }

  return 0;
}

int main() {
  hduPlotter();
  return 0;
}
