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
#include "TStyle.h"
#include "TFile.h"
#include "TSystem.h"
#include <chrono>

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
static constexpr int N_CHANNELS            = 3;  // 0: 1e, 1: 2e, 2: 3e+4e+5e

// Same weighted channels as qWeighted.C (q_1 = k1/N, q_2 = w2*k2/N,
// q_3 = w3*(k3+k4+k5)/N), but -- like multichannel.C, and unlike
// qWeighted.C -- each channel keeps its OWN evolving mask. A channel's
// survival therefore depends only on its own weighted ratio crossing the
// threshold, never on how abundant (or masked) the other channels are.
//
// Return vector layout:
//   [0] surviving 1e events     [3] unmasked pixels, channel 1 (1e)
//   [1] surviving 2e events     [4] unmasked pixels, channel 2 (2e)
//   [2] surviving (3+4+5)e ev.  [5] unmasked pixels, channel 3 (3+4+5e)
vector<int> analyzeSeparate(float threshold,
                             float w2,
                             float w3,
                             string fileName,
                             int ohdu,
                             TH1D* hq1,
                             TH1D* hq2,
                             TH1D* hq3,
                             TH1D* hN,
                             TH1D* h_k1,
                             TH1D* h_k2,
                             TH1D* h_k345,
                             int radius) {

  int n   = 2*radius;
  float r = radius;

  auto [matrix_data, baseMask] = loadImage3(fileName, ohdu);

  // Each channel starts from the same original image mask, but evolves
  // independently as its own weighted ratio crosses threshold.
  vector<vector<vector<int>>> masks(N_CHANNELS, baseMask);

  int Nx = 3200;
  int Ny = 520;

  double weight[N_CHANNELS] = {1.0, static_cast<double>(w2), static_cast<double>(w3)};

  for (int xcoord = 0; xcoord <= Nx - n; xcoord = xcoord + n/2) {
    for (int ycoord = 0; ycoord <= Ny - n; ycoord = ycoord + n/2) {

      int x_center = xcoord + r;
      int y_center = ycoord + r;

      int k[N_CHANNELS]      = {0, 0, 0};
      int trials[N_CHANNELS] = {0, 0, 0};

      for (int i = ycoord; i < ycoord + n; ++i) {
        for (int j = xcoord; j < xcoord + n; ++j) {

          bool isInside = isInsideCircle(r, x_center, y_center, j, i);
          if (!isInside) continue;

          int eventInPixel = matrix_data.at(i).at(j);

          for (int channel = 0; channel < N_CHANNELS; ++channel) {

            int maskInPixel = masks.at(channel).at(i).at(j);

            // Reject pixels with any standard SENSEI mask.
            if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;

            // Reject pixels already masked by a previous window, in
            // THIS channel's own mask.
            if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;

            trials[channel] += 1;

            if (channel == 0 && eventInPixel == 1) k[0] += 1;
            if (channel == 1 && eventInPixel == 2) k[1] += 1;
            if (channel == 2 && (eventInPixel == 3 || eventInPixel == 4 || eventInPixel == 5)) k[2] += 1;
          }
        }
      }

      double q[N_CHANNELS] = {0.0, 0.0, 0.0};
      for (int channel = 0; channel < N_CHANNELS; ++channel) {
        if (trials[channel] == 0) continue;
        q[channel] = weight[channel] * k[channel] / static_cast<double>(trials[channel]);
      }

      if (trials[0] > 0) { hN->Fill(trials[0]); h_k1->Fill(k[0]);   hq1->Fill(q[0]); }
      if (trials[1] > 0) {                      h_k2->Fill(k[1]);   hq2->Fill(q[1]); }
      if (trials[2] > 0) {                      h_k345->Fill(k[2]); hq3->Fill(q[2]); }

      // Apply threshold independently in each channel's own mask.
      for (int channel = 0; channel < N_CHANNELS; ++channel) {

        if (trials[channel] == 0) continue;

        if (q[channel] > threshold) {
          for (int i = ycoord; i < ycoord + n; ++i) {
            for (int j = xcoord; j < xcoord + n; ++j) {

              bool insideMask = isInsideCircle(r, x_center, y_center, j, i);
              if (!insideMask) continue;

              masks.at(channel).at(i).at(j) |= BINOMIAL_MASK_BIT;
            }
          }
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////////
  ////////////////// Count unmasked pixels by channel //////////////////
  /////////////////////////////////////////////////////////////////////

  int nUnmaskedPixels[N_CHANNELS] = {0, 0, 0};
  for (int channel = 0; channel < N_CHANNELS; ++channel) {
    nUnmaskedPixels[channel] = countUnmaskedPixels2(masks.at(channel), true);
  }

  /////////////////////////////////////////////////////////////////////
  ///////// Counting n-electron events that survived each mask /////////
  /////////////////////////////////////////////////////////////////////

  int nEventsUnmasked[N_CHANNELS] = {0, 0, 0};

  for (int x = 0; x < Ny; ++x) {
    for (int y = 0; y < Nx; ++y) {

      int event = matrix_data.at(x).at(y);
      if (event < 1 || event > 5) continue;

      int channel = (event == 1) ? 0 : (event == 2) ? 1 : 2;

      int maskInPixel = masks.at(channel).at(x).at(y);

      if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;
      if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;

      nEventsUnmasked[channel] += 1;
    }
  }

  return {
    nEventsUnmasked[0], nEventsUnmasked[1], nEventsUnmasked[2],
    nUnmaskedPixels[0], nUnmaskedPixels[1], nUnmaskedPixels[2]
  };
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////                Main           /////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Same weights/thresholds/inputs as qWeighted.C, but each channel (1e, 2e,
// 3+4+5e) keeps its own separate mask instead of one combined mask -- see
// analyzeSeparate() above.
//
// Config: ../configs/qweighted_cfg.yaml (w2, w3 -- used only when w2Arg/w3Arg
//         are left at their default, negative sentinel)
//         ../configs/config.yaml        (input_files)
int qWeightedSeparate(int radius = 30, float w2Arg = -1.0f, float w3Arg = -1.0f){

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

  // Same threshold list as qWeighted.C: finer spacing near 0, out to 0.005.
  vector<float> thresholds = {};
  for (float i=0;     i<0.002; i+=0.0001) thresholds.push_back(static_cast<float>(i));
  for (float i=0.002; i<0.005; i+=0.0005) thresholds.push_back(static_cast<float>(i));

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
  cout<<"q_1 = k1/N, q_2 = w2*k2/N, q_3 = w3*(k3+k4+k5)/N -- separate masks per channel"<<endl;
  cout<<"             " <<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"             " <<endl;
  cout<<"             " <<endl;

  string outTag = "qweighted_separate/r" + to_string(radius)
                  + "_w2-" + to_string(static_cast<int>(w2))
                  + "_w3-" + to_string(static_cast<int>(w3));
  string txtDir  = "../txts/"    + outTag;
  string pdfDir  = "../results/" + outTag;
  string rootDir = "../root/"    + outTag;
  gSystem->mkdir(txtDir.c_str(), true);
  gSystem->mkdir(pdfDir.c_str(), true);
  gSystem->mkdir(rootDir.c_str(), true);

  ////////////////////////////////////////////////////////////////////////////
  // Threshold scan: masking uses a separate mask per channel.
  ////////////////////////////////////////////////////////////////////////////

  // Histograms below are only needed to satisfy analyzeSeparate()'s
  // signature during the scan; their content mixes every threshold/file
  // together and is not written out. The threshold=0.1 diagnostic pass
  // further down makes fresh histograms of its own for the actual plots.
  TH1D* hN_scan   = new TH1D("","Distribution of N per window",     50, 0, 3000);
  TH1D* hK1_scan  = new TH1D("","Distribution of k_{1} per window", 20, 0, 20);
  TH1D* hK2_scan  = new TH1D("","Distribution of k_{2} per window", 20, 0, 20);
  TH1D* hK345_scan= new TH1D("","Distribution of k_{3+4+5} per window", 20, 0, 20);
  TH1D* hQ1_scan  = new TH1D("","Distribution of q_1", 200, 0, 0.10);
  TH1D* hQ2_scan  = new TH1D("","Distribution of q_2", 200, 0, 0.10);
  TH1D* hQ3_scan  = new TH1D("","Distribution of q_3", 200, 0, 0.10);

  // Clock
  auto aScan = high_resolution_clock::now();

  vector<float> n1EventsVec   = {};
  vector<float> n2EventsVec   = {};
  vector<float> n345EventsVec = {};
  vector<float> pix1Vec = {};
  vector<float> pix2Vec = {};
  vector<float> pix345Vec = {};
  vector<float> rate1Vec   = {};
  vector<float> rate2Vec   = {};
  vector<float> rate345Vec = {};

  for (auto wp : thresholds) {
    cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Threshold = "<< wp <<endl;

    float n1Events   = 0;
    float n2Events   = 0;
    float n345Events = 0;
    float pix1   = 0;
    float pix2   = 0;
    float pix345 = 0;

    for (const auto& rf : inputFiles) {

      int lta = extractLTA(rf);
      auto ltaIt = GOOD_QUADS.find(lta);
      if (ltaIt == GOOD_QUADS.end()) continue;
      const vector<int>& fileOhdus = ltaIt->second;

      for (int oh=0; oh<fileOhdus.size(); oh++){
        int ohdu = fileOhdus.at(oh);

        vector<int> qs = analyzeSeparate(wp, w2, w3, rf, ohdu,
                                          hQ1_scan, hQ2_scan, hQ3_scan,
                                          hN_scan, hK1_scan, hK2_scan, hK345_scan, radius);

        n1Events   += qs.at(0);
        n2Events   += qs.at(1);
        n345Events += qs.at(2);
        pix1       += qs.at(3);
        pix2       += qs.at(4);
        pix345     += qs.at(5);
      }
    }

    n1EventsVec.push_back(n1Events);
    n2EventsVec.push_back(n2Events);
    n345EventsVec.push_back(n345Events);
    pix1Vec.push_back(pix1);
    pix2Vec.push_back(pix2);
    pix345Vec.push_back(pix345);

    rate1Vec.push_back(  pix1   > 0 ? n1Events/pix1     : 0.0f);
    rate2Vec.push_back(  pix2   > 0 ? n2Events/pix2     : 0.0f);
    rate345Vec.push_back(pix345 > 0 ? n345Events/pix345 : 0.0f);
  }

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

    // Unlike qWeighted.C (one shared mask -> one pixel count), each
    // channel here has its own mask, so its own pixel count.
    outFile1 << thresholds.at(t) << "\t"
             << pix1Vec.at(t) << "\t"
             << pix2Vec.at(t) << "\t"
             << pix345Vec.at(t) << endl;

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

  TH1D* hN    = new TH1D(("hN_"   +tag).c_str(), "Distribution of N per window",         100, 0,    3000);
  TH1D* hK1   = new TH1D(("hK1_"  +tag).c_str(), "Distribution of k_{1} per window",     20,  0,    20);
  TH1D* hK2   = new TH1D(("hK2_"  +tag).c_str(), "Distribution of k_{2} per window",     20,  0,    20);
  TH1D* hK345 = new TH1D(("hK345_"+tag).c_str(), "Distribution of k_{3+4+5} per window", 20,  0,    20);
  TH1D* hQ1   = new TH1D(("hQ1_"  +tag).c_str(), "Distribution of q_{1} = #frac{k_{1}}{N}",             200, 0, 0.10);
  TH1D* hQ2   = new TH1D(("hQ2_"  +tag).c_str(), "Distribution of q_{2} = #frac{w_{2}k_{2}}{N}",        200, 0, 0.10);
  TH1D* hQ3   = new TH1D(("hQ3_"  +tag).c_str(), "Distribution of q_{3} = #frac{w_{3}(k_{3}+k_{4}+k_{5})}{N}", 200, 0, 0.10);

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

      analyzeSeparate(extraThreshold, w2, w3, rf, ohdu, hQ1, hQ2, hQ3, hN, hK1, hK2, hK345, radius);
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
  canvasMaker(hQ1,   "q_{1} = #frac{k_{1}}{N}",                        (pdfDir + "/q1_statistic_thr0p1.pdf").c_str());
  canvasMaker(hQ2,   "q_{2} = #frac{w_{2}k_{2}}{N}",                   (pdfDir + "/q2_statistic_thr0p1.pdf").c_str());
  canvasMaker(hQ3,   "q_{3} = #frac{w_{3}(k_{3}+k_{4}+k_{5})}{N}",     (pdfDir + "/q3_statistic_thr0p1.pdf").c_str());

  TFile* outRootFile = new TFile((rootDir + "/qWeightedSeparate_thr0p1.root").c_str(), "RECREATE");
  hN     ->Write();
  hK1    ->Write();
  hK2    ->Write();
  hK345  ->Write();
  hQ1    ->Write();
  hQ2    ->Write();
  hQ3    ->Write();
  outRootFile->Close();

  return 0;
}
