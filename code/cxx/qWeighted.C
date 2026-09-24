#include <yaml-cpp/yaml.h>
#include "helperFunctions.C"
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cmath>
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include "TSystem.h"
#include <chrono>
#include <sstream>

using namespace chrono;
using namespace std;

// Extract the LTA number from a filename of the form
// "..._EXPOSURE216000_<LTA>_<runId>.root"
int extractLTA(const string& fileName) {
  size_t dot = fileName.find_last_of('.');
  size_t lastUnderscore = fileName.find_last_of('_', dot - 1);
  size_t ltaUnderscore = fileName.find_last_of('_', lastUnderscore - 1);
  string ltaStr = fileName.substr(ltaUnderscore + 1, lastUnderscore - ltaUnderscore - 1);
  return stoi(ltaStr);
}

// Extract the runId from the same filename pattern.
int extractRunId(const string& fileName) {
  size_t dot = fileName.find_last_of('.');
  size_t lastUnderscore = fileName.find_last_of('_', dot - 1);
  return stoi(fileName.substr(lastUnderscore + 1, dot - lastUnderscore - 1));
}

// Only these LTAs/ohdus (the "good quads") are processed; any other LTA is skipped.
static const map<int, vector<int>> GOOD_QUADS = {
  {3,  {2}},
  {9,  {1, 2}},
  {11, {1, 2, 3}},
  {13, {2, 3, 4}},
  {14, {1, 3}},
  {15, {2}},
  {16, {3}},
  {17, {1, 2, 3, 4}},
  {18, {1, 2, 3}},
};

// Dynamic mask bit used by this algorithm.
// It is OR-ed with the pre-existing mask so previous mask bits are preserved.
static constexpr int BINOMIAL_MASK_BIT     = 65536;
static constexpr int PREEXISTING_MASK_BITS = 0x56FC;

// Builds the distribution of the single combined statistic
//   q = (k1 + w2*k2 + w3*(k3 + k4 + k5)) / N
// over the sliding-window grid of one quadrant, and returns the events/
// pixels that survive the single combined-q mask built along the way.
//
// Unlike multichannel.C's analyze() (one mask per channel, so each
// channel's q can be thresholded independently), this statistic mixes
// all channels together, so there is only one mask: a window that trips
// the q > threshold cut is masked out of every later window's N,
// regardless of which channel drove q up. Same sequential bookkeeping
// as qAppletExport.C / q1MaskGif.C otherwise.
//
// Return vector layout:
//   [0] surviving 1e events
//   [1] surviving 2e events
//   [2] surviving (3+4+5)e events
//   [3] unmasked pixels (single mask, shared by every channel)
vector<int> analyzeWeighted(float threshold,
                             float w2,
                             float w3,
                             string fileName,
                             int ohdu,
                             TH1D* hq,
                             TH1D* hN,
                             TH1D* h_k1,
                             TH1D* h_k2,
                             TH1D* h_k345,
                             TH1D* hK,
                             TH1D* h_k1w,
                             TH1D* h_k2w,
                             TH1D* h_k345w,
                             int radius) {

  int n   = 2*radius;
  float r = radius;

  auto [matrix_data, baseMask] = loadImage3(fileName, ohdu);

  vector<vector<int>> mask = baseMask;

  int Nx = 3200;
  int Ny = 520;

  for (int xcoord = 0; xcoord <= Nx - n; xcoord = xcoord + n/2) {
    for (int ycoord = 0; ycoord <= Ny - n; ycoord = ycoord + n/2) {

      // Coordinates of the circle center
      int x_center = xcoord + r;
      int y_center = ycoord + r;

      int k1   = 0;
      int k2   = 0;
      int k345 = 0;
      int trials = 0;

      for (int i = ycoord; i < ycoord + n; ++i) {
        for (int j = xcoord; j < xcoord + n; ++j) {

          bool isInside = isInsideCircle(r, x_center, y_center, j, i);
          if (!isInside) continue;

          int maskInPixel = mask.at(i).at(j);

          // Reject pixels with any standard SENSEI mask.
          if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;

          // Reject pixels already masked by a previous window.
          if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;

          trials += 1;

          int eventInPixel = matrix_data.at(i).at(j);
          if (eventInPixel == 1) {
            k1 += 1;
          } else if (eventInPixel == 2) {
            k2 += 1;
          } else if (eventInPixel == 3 || eventInPixel == 4 || eventInPixel == 5) {
            k345 += 1;
          }
        }
      }

      // Fully masked window: nothing left to compute a statistic from.
      if (trials == 0) continue;

      long double k1w   = static_cast<long double>(k1);
      long double k2w   = static_cast<long double>(w2) * k2;
      long double k345w = static_cast<long double>(w3) * k345;
      long double K     = k1w + k2w + k345w;
      long double q     = K / static_cast<long double>(trials);

      hN    ->Fill(trials);
      h_k1  ->Fill(k1);
      h_k2  ->Fill(k2);
      h_k345->Fill(k345);
      hq    ->Fill(static_cast<double>(q));

      // Combined statistic (before dividing by N) and its per-channel
      // weighted contributions, so the three can be stacked in one plot.
      hK     ->Fill(static_cast<double>(K));
      h_k1w  ->Fill(static_cast<double>(k1w));
      h_k2w  ->Fill(static_cast<double>(k2w));
      h_k345w->Fill(static_cast<double>(k345w));

      if (q > threshold) {
        for (int i = ycoord; i < ycoord + n; ++i) {
          for (int j = xcoord; j < xcoord + n; ++j) {

            bool insideMask = isInsideCircle(r, x_center, y_center, j, i);
            if (!insideMask) continue;

            // Preserve all existing mask bits and add the binomial bit.
            mask.at(i).at(j) |= BINOMIAL_MASK_BIT;
          }
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////
  ///////// Counting n-electron events that survived the mask //////////
  /////////////////////////////////////////////////////////////////////

  int nUnmaskedPixels = countUnmaskedPixels2(mask, true);

  int n1EventsUnmasked   = 0;
  int n2EventsUnmasked   = 0;
  int n345EventsUnmasked = 0;

  for (int x = 0; x < Ny; ++x) {
    for (int y = 0; y < Nx; ++y) {

      int event = matrix_data.at(x).at(y);
      if (event < 1 || event > 5) continue;

      int maskInPixel = mask.at(x).at(y);

      // Standard SENSEI masks.
      if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;

      // Binomial mask.
      if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;

      if (event == 1) {
        n1EventsUnmasked += 1;
      } else if (event == 2) {
        n2EventsUnmasked += 1;
      } else {
        n345EventsUnmasked += 1;
      }
    }
  }

  return {n1EventsUnmasked, n2EventsUnmasked, n345EventsUnmasked, nUnmaskedPixels};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////                Main           /////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Config: ../configs/qweighted_cfg.yaml (w2, w3 -- used only when w2Arg/w3Arg
//         are left at their default, negative sentinel)
//         ../configs/config.yaml        (input_files)
//
// w2Arg/w3Arg let several qWeighted(radius, w2, w3) jobs run at once (e.g.
// a weight scan) without each overwriting the same shared config file.
int qWeighted(int radius = 30, float w2Arg = -1.0f, float w3Arg = -1.0f){

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
  gStyle->SetTitleFontSize(0.025);

  // Make list of thresholds/cut values to loop over: finer spacing near
  // the start (where the survival curves move fastest), coarser further
  // out, up to 0.02 -- past the old 0.005 ceiling, since at w2/w3=10 a
  // single 2e/3+e hit alone contributes ~10/N (~0.02) to q.
  vector<float> thresholds = {};
  for (float i=0.0012; i<0.002;  i+=0.0001) thresholds.push_back(static_cast<float>(i));
  for (float i=0.002;  i<0.005;  i+=0.0005) thresholds.push_back(static_cast<float>(i));
  for (float i=0.005;  i<0.0201; i+=0.005)  thresholds.push_back(static_cast<float>(i));

  YAML::Node configFile = YAML::LoadFile("../configs/config.yaml");
  vector<string> inputFiles = configFile["input_files"].as<vector<string>>();

  YAML::Node wConfig = YAML::LoadFile("../configs/qweighted_cfg.yaml");
  float w2 = (w2Arg >= 0.0f) ? w2Arg : wConfig["w2"].as<float>();
  float w3 = (w3Arg >= 0.0f) ? w3Arg : wConfig["w3"].as<float>();

  cout<<"             " <<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"             " <<endl;
  cout<<"Input files: " <<endl;
  for (const auto& file : inputFiles) cout <<"  " <<file << endl;
  cout<<"             " <<endl;
  cout<<"Radius: "   <<radius   <<endl;
  cout<<"w2: "<<w2<<"   w3: "<<w3<<endl;
  cout<<"q = (k1 + w2*k2 + w3*(k3+k4+k5)) / N"<<endl;
  cout<<"             " <<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"             " <<endl;
  cout<<"             " <<endl;

  ostringstream w2Str, w3Str;
  w2Str << w2;
  w3Str << w3;
  string outTag = "qweighted/r" + to_string(radius)
                  + "_w2-" + w2Str.str()
                  + "_w3-" + w3Str.str();
  string txtDir  = "../txts/"    + outTag;
  string pdfDir  = "../results/" + outTag;
  string rootDir = "../root/"    + outTag;
  gSystem->mkdir(txtDir.c_str(), true);
  gSystem->mkdir(pdfDir.c_str(), true);
  gSystem->mkdir(rootDir.c_str(), true);

  ////////////////////////////////////////////////////////////////////////////
  // Threshold scan: masking uses only the single combined q, so every
  // channel shares the same mask/N at a given threshold (unlike
  // multichannel.C, which masks each channel independently).
  ////////////////////////////////////////////////////////////////////////////

  // Histograms below are only needed to satisfy analyzeWeighted()'s
  // signature during the scan; their content mixes every threshold/file
  // together and is not written out. The threshold=0.1 diagnostic pass
  // further down makes fresh histograms of its own for the actual plots.
  TH1D* hN_scan     = new TH1D("","Distribution of N per window",     50, 0, 3000);
  TH1D* hK1_scan    = new TH1D("","Distribution of k_{1} per window", 20, 0, 20);
  TH1D* hK2_scan    = new TH1D("","Distribution of k_{2} per window", 20, 0, 20);
  TH1D* hK345_scan  = new TH1D("","Distribution of k_{3+4+5} per window", 20, 0, 20);
  TH1D* hQ_scan     = new TH1D("","Distribution of q", 200, 0, 0.10);
  TH1D* hK_scan     = new TH1D("","Distribution of K", 40, 0, 20);
  TH1D* hK1w_scan   = new TH1D("","Contribution of k_{1}", 40, 0, 20);
  TH1D* hK2w_scan   = new TH1D("","Contribution of w_{2}k_{2}", 40, 0, 20);
  TH1D* hK345w_scan = new TH1D("","Contribution of w_{3}(k_{3}+k_{4}+k_{5})", 40, 0, 20);

  // Clock
  auto aScan = high_resolution_clock::now();

  vector<float> n1EventsVec   = {};
  vector<float> n2EventsVec   = {};
  vector<float> n345EventsVec = {};
  vector<float> unmaskedPixelsVec = {};
  vector<float> rate1Vec   = {};
  vector<float> rate2Vec   = {};
  vector<float> rate345Vec = {};

  // Per-image survivors at every threshold, one row per (threshold, file, ohdu).
  std::ofstream outPerImage((txtDir + "/qweighted_perImage.txt").c_str());
  if (!outPerImage) { cerr << "Error opening file for writing!" << endl; }
  outPerImage << "#threshold\tlta\tohdu\trunId\tn1\tn2\tn345\tnPix" << endl;

  for (auto wp : thresholds) {
    cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Threshold = "<< wp <<endl;

    float n1Events   = 0;
    float n2Events   = 0;
    float n345Events = 0;
    float unmaskedPixels = 0;

    for (const auto& rf : inputFiles) {

      int lta = extractLTA(rf);
      auto ltaIt = GOOD_QUADS.find(lta);
      if (ltaIt == GOOD_QUADS.end()) continue;
      const vector<int>& fileOhdus = ltaIt->second;

      for (int oh=0; oh<fileOhdus.size(); oh++){
        int ohdu = fileOhdus.at(oh);

        vector<int> qs = analyzeWeighted(wp, w2, w3, rf, ohdu,
                                          hQ_scan, hN_scan, hK1_scan, hK2_scan, hK345_scan,
                                          hK_scan, hK1w_scan, hK2w_scan, hK345w_scan, radius);

        outPerImage << wp << "\t" << lta << "\t" << ohdu << "\t" << extractRunId(rf)
                    << "\t" << qs.at(0) << "\t" << qs.at(1) << "\t" << qs.at(2)
                    << "\t" << qs.at(3) << endl;

        n1Events       += qs.at(0);
        n2Events       += qs.at(1);
        n345Events     += qs.at(2);
        unmaskedPixels += qs.at(3);
      }
    }

    n1EventsVec.push_back(n1Events);
    n2EventsVec.push_back(n2Events);
    n345EventsVec.push_back(n345Events);
    unmaskedPixelsVec.push_back(unmaskedPixels);

    rate1Vec.push_back(  unmaskedPixels > 0 ? n1Events/unmaskedPixels   : 0.0f);
    rate2Vec.push_back(  unmaskedPixels > 0 ? n2Events/unmaskedPixels   : 0.0f);
    rate345Vec.push_back(unmaskedPixels > 0 ? n345Events/unmaskedPixels : 0.0f);
  }

  outPerImage.close();

  // chrono
  auto bScan = high_resolution_clock::now();
  cout << " "<<  endl;
  cout << "================== "<<  endl;
  cout << "Scan took " << duration_cast<seconds>(bScan - aScan).count() << " seconds" <<  endl;
  cout << "================== "<<  endl;
  cout << " "<<  endl;

  // output txt to plot with python
  std::ofstream outFile0((txtDir + "/qweighted_survElec.txt").c_str());
  std::ofstream outFile1((txtDir + "/qweighted_survPix.txt").c_str());
  std::ofstream outFile2((txtDir + "/qweighted_rates.txt").c_str());

  if (!outFile0) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile1) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile2) { cerr << "Error opening file for writing!" << endl; }

  for (size_t t = 0; t < thresholds.size(); ++t) {
    outFile0 << thresholds.at(t) << "\t"
             << n1EventsVec.at(t) << "\t"
             << n2EventsVec.at(t) << "\t"
             << n345EventsVec.at(t) << endl;

    outFile1 << thresholds.at(t) << "\t"
             << unmaskedPixelsVec.at(t) << endl;

    outFile2 << thresholds.at(t) << "\t"
             << rate1Vec.at(t) << "\t"
             << rate2Vec.at(t) << "\t"
             << rate345Vec.at(t) << endl;
  }

  outFile0.close();
  outFile1.close();
  outFile2.close();

  ////////////////////////////////////////////////////////////////////////////
  // Extra diagnostic run at threshold = 0.1.
  // Not written to txts: produces pdf plots and a root file with the
  // histograms instead.
  ////////////////////////////////////////////////////////////////////////////

  cout << " "<<  endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"Extra run at threshold = 0.1 (pdfs + root file only, no txt output)"<<endl;
  cout<<"==================================================================================================================================================================================="<<endl;

  float extraThreshold = 0.1f;

  string tag = "thr0p1_r" + to_string(radius)
               + "_w2-" + to_string(static_cast<int>(w2))
               + "_w3-" + to_string(static_cast<int>(w3));

  TH1D* hN    = new TH1D(("hN_"   +tag).c_str(), "Distribution of N per window",           100, 0,    3000);
  TH1D* hK1   = new TH1D(("hK1_"  +tag).c_str(), "Distribution of k_{1} per window",       20,  0,    20);
  TH1D* hK2   = new TH1D(("hK2_"  +tag).c_str(), "Distribution of k_{2} per window",       20,  0,    20);
  TH1D* hK345 = new TH1D(("hK345_"+tag).c_str(), "Distribution of k_{3+4+5} per window",   20,  0,    20);
  TH1D* hQ    = new TH1D(("hQ_"   +tag).c_str(), "Distribution of q = #frac{k_{1}+w_{2}k_{2}+w_{3}(k_{3}+k_{4}+k_{5})}{N}",
                                                                                            200, 0,    0.10);

  // Combined statistic K = k1 + w2*k2 + w3*(k3+k4+k5) (i.e. q*N), plus its
  // three per-channel weighted contributions, all on the same binning so
  // they can be stacked in one plot.
  TH1D* hK      = new TH1D(("hK_"     +tag).c_str(), "Combined K = k_{1}+w_{2}k_{2}+w_{3}(k_{3}+k_{4}+k_{5}) per window", 40, 0, 20);
  TH1D* hK1w    = new TH1D(("hK1w_"   +tag).c_str(), "Contribution of k_{1}",                                            40, 0, 20);
  TH1D* hK2w    = new TH1D(("hK2w_"   +tag).c_str(), "Contribution of w_{2}k_{2}",                                       40, 0, 20);
  TH1D* hK345w  = new TH1D(("hK345w_" +tag).c_str(), "Contribution of w_{3}(k_{3}+k_{4}+k_{5})",                        40, 0, 20);

  // Clock
  auto aDiag = high_resolution_clock::now();

  for (const auto& rf : inputFiles) {

    int lta = extractLTA(rf);
    auto ltaIt = GOOD_QUADS.find(lta);
    if (ltaIt == GOOD_QUADS.end()) {
      cout<<"                   *skipping LTA "<< lta <<" (not in GOOD_QUADS)"<<endl;
      continue;
    }
    const vector<int>& fileOhdus = ltaIt->second;

    cout<<"  READING FILE: "<< rf << endl;

    for (int oh=0; oh<fileOhdus.size(); oh++){
      int ohdu = fileOhdus.at(oh);

      cout<<"                   *ohdu "<< ohdu << endl;

      analyzeWeighted(extraThreshold, w2, w3, rf, ohdu, hQ, hN, hK1, hK2, hK345,
                       hK, hK1w, hK2w, hK345w, radius);
    }
  }

  // chrono
  auto bDiag = high_resolution_clock::now();
  cout << " "<<  endl;
  cout << "================== "<<  endl;
  cout << "Diagnostic pass took " << duration_cast<seconds>(bDiag - aDiag).count() << " seconds" <<  endl;
  cout << "================== "<<  endl;
  cout << " "<<  endl;

  canvasMaker(hN,    "Distribution of N",         (pdfDir + "/hN_thr0p1.pdf").c_str());
  canvasMaker(hK1,   "Distribution of k_{1}",      (pdfDir + "/k1_thr0p1.pdf").c_str());
  canvasMaker(hK2,   "Distribution of k_{2}",      (pdfDir + "/k2_thr0p1.pdf").c_str());
  canvasMaker(hK345, "Distribution of k_{3+4+5}",  (pdfDir + "/k345_thr0p1.pdf").c_str());
  canvasMaker(hQ,    "q statistic",                (pdfDir + "/q_statistic_thr0p1.pdf").c_str());

  ////////////////////////////////////////////////////////////////////////////
  // Combined K's three weighted contributions, stacked.
  ////////////////////////////////////////////////////////////////////////////

  hK1w  ->SetFillColor(kAzure+1);   hK1w  ->SetLineColor(kAzure+1);
  hK2w  ->SetFillColor(kGreen+2);   hK2w  ->SetLineColor(kGreen+2);
  hK345w->SetFillColor(kRed+1);     hK345w->SetLineColor(kRed+1);
  hK1w->SetFillStyle(1001);   hK2w->SetFillStyle(1001);   hK345w->SetFillStyle(1001);

  THStack* stackK = new THStack(("stackK_"+tag).c_str(),
                                 "Combined K = k_{1}+w_{2}k_{2}+w_{3}(k_{3}+k_{4}+k_{5})");
  stackK->Add(hK1w);
  stackK->Add(hK2w);
  stackK->Add(hK345w);

  TCanvas* canvK = new TCanvas("canvK");
  canvK->SetFillColor(kBlack);
  canvK->SetBatch(kTRUE);
  canvK->SetLogy();

  stackK->Draw("HIST");
  stackK->GetXaxis()->SetTitle("K");
  stackK->GetYaxis()->SetTitle("Number of entries");
  stackK->GetXaxis()->SetAxisColor(kWhite);
  stackK->GetYaxis()->SetAxisColor(kWhite);
  stackK->GetXaxis()->SetLabelColor(kWhite);
  stackK->GetYaxis()->SetLabelColor(kWhite);
  stackK->GetXaxis()->SetTitleColor(kWhite);
  stackK->GetYaxis()->SetTitleColor(kWhite);

  TLegend* legK = new TLegend(0.6, 0.7, 0.88, 0.88);
  legK->SetFillColor(kBlack);
  legK->SetTextColor(kWhite);
  legK->SetLineColor(kWhite);
  legK->AddEntry(hK1w,   "k_{1}",                      "f");
  legK->AddEntry(hK2w,   "w_{2}k_{2}",                 "f");
  legK->AddEntry(hK345w, "w_{3}(k_{3}+k_{4}+k_{5})",   "f");
  legK->Draw();

  canvK->SaveAs((pdfDir + "/K_combined_stack_thr0p1.pdf").c_str());
  canvK->Close();

  TFile* outRootFile = new TFile((rootDir + "/qWeighted_thr0p1.root").c_str(), "RECREATE");
  hN     ->Write();
  hK1    ->Write();
  hK2    ->Write();
  hK345  ->Write();
  hQ     ->Write();
  hK     ->Write();
  hK1w   ->Write();
  hK2w   ->Write();
  hK345w ->Write();
  stackK ->Write();
  outRootFile->Close();

  return 0;
}
