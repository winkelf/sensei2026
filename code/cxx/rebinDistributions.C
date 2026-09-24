// Redo the multichannel.C threshold=0.1 diagnostic distribution plots
// (k1, k2, k3, N, q_1e..q_4e) for each scanned radius, but at coarser
// binning: half of the original number of bins, and half of that again
// (a quarter of the original). Rebin(2)/Rebin(4) merges adjacent bins of
// the histograms already filled and saved by multichannel.C, which is
// exactly equivalent to having filled with half/quarter as many
// same-range bins to begin with -- no need to re-scan the FITS images.
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <iostream>

using namespace std;

void styleForDark(){
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
  gStyle->SetHistFillColor(kBlack);
  gStyle->SetHistLineColor(kWhite);
  gStyle->SetFuncColor(kWhite);
  gStyle->SetFrameLineColor(kWhite);
  gStyle->SetTitleFontSize(0.025);
}

void canvasMakerRebin(TH1D* hist, const char* title, const char* pdfName){
  TCanvas* canv = new TCanvas(title);
  canv->SetFillColor(kBlack);
  canv->SetBatch(kTRUE);
  hist->SetStats(0);
  hist->GetXaxis()->SetTitle(title);
  hist->GetYaxis()->SetTitle("Number of entries");
  hist->GetYaxis()->SetAxisColor(kWhite);
  hist->GetXaxis()->SetAxisColor(kWhite);
  hist->Draw("HIST");
  canv->SetLogy();
  canv->SaveAs(pdfName);
  canv->Close();
}

void rebinDistributions(){

  styleForDark();

  vector<int> radii = {15, 20, 25, 30};

  // {histName in file, x-axis title, output base name}
  vector<std::array<string,3>> hists = {
    {"hN_01",    "Distribution of N",     "hN_thr0p1"},
    {"h_k1_01",  "Distribution of k_{1}", "k1_thr0p1"},
    {"h_k2_01",  "Distribution of k_{2}", "k2_thr0p1"},
    {"h_k3_01",  "Distribution of k_{3}", "k3_thr0p1"},
    {"hq_1e_01", "q = #frac{k_{1}}{N}",   "q_statistic_1e_thr0p1"},
    {"hq_2e_01", "q = #frac{k_{2}}{N}",   "q_statistic_2e_thr0p1"},
    {"hq_3e_01", "q = #frac{k_{3}}{N}",   "q_statistic_3e_thr0p1"},
    {"hq_4e_01", "q = #frac{k_{4}}{N}",   "q_statistic_4e_thr0p1"},
  };

  // {rebin group size, suffix, human label} -- applied to the ORIGINAL
  // (already saved) histogram each time, so both derived versions share
  // the same original bin edges/range.
  vector<std::pair<int,string>> rebins = {
    {2, "_half"},
    {4, "_quarter"},
  };

  for (int radius : radii) {
    string inPath  = "../root/multichannel/r" + to_string(radius) + "/multichannel_thr0p1.root";
    string outDir  = "../results/multichannel/r" + to_string(radius);
    gSystem->mkdir(outDir.c_str(), true);

    TFile* f = TFile::Open(inPath.c_str(), "READ");
    if (!f || f->IsZombie()) {
      cout << "Could not open " << inPath << ", skipping r" << radius << endl;
      continue;
    }

    cout << "== r" << radius << " ==" << endl;

    for (const auto& h : hists) {
      const string& histName = h[0];
      const string& axisTitle = h[1];
      const string& baseName = h[2];

      TH1D* orig = (TH1D*)f->Get(histName.c_str());
      if (!orig) {
        cout << "  missing histogram " << histName << " in " << inPath << endl;
        continue;
      }

      for (const auto& rb : rebins) {
        int groupSize = rb.first;
        const string& suffix = rb.second;

        TH1D* clone = (TH1D*)orig->Clone((histName + suffix).c_str());
        clone->Rebin(groupSize);

        string pdfName = outDir + "/" + baseName + suffix + ".pdf";
        canvasMakerRebin(clone, axisTitle.c_str(), pdfName.c_str());
        cout << "  wrote " << pdfName
             << " (" << clone->GetNbinsX() << " bins)" << endl;

        delete clone;
      }
    }

    f->Close();
    delete f;
  }
}
