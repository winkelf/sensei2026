#include <yaml-cpp/yaml.h>
#include "helperFunctions_fixed.C"
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
#include <nlohmann/json.hpp>
#include <chrono>

using namespace chrono;
using namespace std;

// Dynamic mask bit used by this algorithm.
// It is OR-ed with the pre-existing mask so previous mask bits are preserved.
static constexpr int BINOMIAL_MASK_BIT = 65536;
static constexpr int N_CHANNELS = 4;
static constexpr int PREEXISTING_MASK_BITS = 0x56FC;

// Return vector layout:
//   [0] surviving 1e events
//   [1] surviving 2e events
//   [2] surviving 3e events
//   [3] surviving 4e events
//   [4] unmasked pixels for the 1e mask
//   [5] unmasked pixels for the 2e mask
//   [6] unmasked pixels for the 3e mask
//   [7] unmasked pixels for the 4e mask
vector<int> analyze(float threshold,
                    string fileName,
                    int ohdu,
                    TH1D* hq_1e,
                    TH1D* hq_2e,
                    TH1D* hq_3e,
                    TH1D* hq_4e,
                    TH1D* hN,
                    TH1D* h_k1,
                    TH1D* h_k2,
                    TH1D* h_k3) {

  int n   = 60;
  float r = n/2;

  ////////////////////////////////////////////////////////////////////
  ////////////////////// Get previous masks (only for data) //////////
  ////////////////////////////////////////////////////////////////////

  auto [matrix_data, baseMask] = loadImage3(fileName, ohdu);
  //auto [matrix_data, baseMask] = loadImage2(fileName, ohdu, {});

  // Each electron channel starts from the same original image mask,
  // but evolves independently during the threshold scan.
  vector<vector<vector<int>>> previousMasks(N_CHANNELS, baseMask);

  ////////////////////////////////////////////////////////////////////
  ////////////////////// Scan grid ///////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  
  int Nx = 3200;
  int Ny = 520;
  
  for (int xcoord = 0; xcoord <= Nx - n; xcoord = xcoord + n/2) {
    for (int ycoord = 0; ycoord <= Ny - n; ycoord = ycoord + n/2) {
  
      // One k and one number of trials for each electron channel.
      // Once the masks diverge, the effective N can be different in each channel.
      int k[N_CHANNELS]      = {0, 0, 0, 0};
      int trials[N_CHANNELS] = {0, 0, 0, 0};
  
      // Coordinates of the circle center
      int x_center = xcoord + r;
      int y_center = ycoord + r;
  
      // Count hits in sliding window independently for every channel mask
      for (int i = ycoord; i < ycoord + n; ++i) {
        for (int j = xcoord; j < xcoord + n; ++j) {
  
          int eventInPixel = matrix_data.at(i).at(j);
  
          // Check if point is inside circle
          bool isInside = isInsideCircle(r, x_center, y_center, j, i);
          if (!isInside) continue;
  
          for (int channel = 0; channel < N_CHANNELS; ++channel) {
  
            int maskInPixel =
                int(previousMasks.at(channel).at(i).at(j));
  
            // Reject pixels with any standard SENSEI mask.
            if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;
  
            // Reject pixels already masked by a previous binomial window.
            if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;
  
            // Count number of valid binomial trials.
            trials[channel] += 1;
  
            // Count events belonging to this electron channel.
            if (eventInPixel == channel + 1) {
              k[channel] += 1;
            }
          }
        }
      }
  
      long double P[N_CHANNELS] = {0.0L, 0.0L, 0.0L, 0.0L};
  
      for (int channel = 0; channel < N_CHANNELS; ++channel) {
  
        if (trials[channel] == 0) continue;
  
        P[channel] = static_cast<long double>(k[channel]) / static_cast<long double>(trials[channel]);
      }
  
      // Diagnostic histograms
      if (trials[0] > 0) hN->Fill(trials[0]);
  
      h_k1->Fill(k[0]);
      h_k2->Fill(k[1]);
      h_k3->Fill(k[2]);
  
      if (trials[0] > 0) hq_1e->Fill(P[0]);
      if (trials[1] > 0) hq_2e->Fill(P[1]);
      if (trials[2] > 0) hq_3e->Fill(P[2]);
      if (trials[3] > 0) hq_4e->Fill(P[3]);
  
      // Apply threshold independently in each channel.
      for (int channel = 0; channel < N_CHANNELS; ++channel) {
  
        // No usable pixels in this window.
        if (trials[channel] == 0) continue;
  
        if (P[channel] > threshold) {
  
          // Set BINOMIAL_MASK_BIT in this channel only.
          for (int i = ycoord; i < ycoord + n; ++i) {
            for (int j = xcoord; j < xcoord + n; ++j) {
  
              bool isInside = isInsideCircle(r, x_center, y_center, j, i);
  
              if (!isInside) continue;
  
              // Preserve all existing mask bits and add the binomial bit.
              previousMasks.at(channel).at(i).at(j) |= BINOMIAL_MASK_BIT;
            }
          }
        }
      }
    }
  }
  
  /////////////////////////////////////////////////////////////////////
  ////////////////// Count unmasked pixels by channel //////////////////
  /////////////////////////////////////////////////////////////////////
  
  int nUnmaskedPixels[N_CHANNELS] = {0, 0, 0, 0};
  
  for (int channel = 0; channel < N_CHANNELS; ++channel) {
    nUnmaskedPixels[channel] = countUnmaskedPixels2(previousMasks.at(channel), true);
  }
  
  /////////////////////////////////////////////////////////////////////
  ///////// Counting n-electron events that survived each mask /////////
  /////////////////////////////////////////////////////////////////////
  
  int nEventsUnmasked[N_CHANNELS] = {0, 0, 0, 0};
  
  for (int x = 0; x < 520; ++x) {
    for (int y = 0; y < 3200; ++y) {
  
      int event = matrix_data.at(x).at(y);
  
      // An ne event is tested only against the mask constructed
      // for that same ne channel.
      for (int channel = 0; channel < N_CHANNELS; ++channel) {
  
        if (event != channel + 1) continue;
  
        int maskInPixel = previousMasks.at(channel).at(x).at(y);
  
        // Standard SENSEI masks.
        if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;
  
        // Binomial mask.
        if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;
  
        nEventsUnmasked[channel] += 1;
      }
    }
  }

  vector<int> rVector = {
    nEventsUnmasked[0],
    nEventsUnmasked[1],
    nEventsUnmasked[2],
    nEventsUnmasked[3],
    nUnmaskedPixels[0],
    nUnmaskedPixels[1],
    nUnmaskedPixels[2],
    nUnmaskedPixels[3]
  };

  return rVector;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////                Main           /////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int multichannel(){

  // Make list of thresholds/cut values to loop over
  //vector<float> thresholds = {0.1};
  vector<float> thresholds = {};
  //for (float i=0; i<0.01; i+=0.0005) thresholds.push_back(static_cast<float>(i));
  for (float i=0; i<0.01; i+=0.001) thresholds.push_back(static_cast<float>(i));

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

  // input files
  YAML::Node configFile = YAML::LoadFile("../config/config.yaml");
  vector<string> inputFiles = configFile["input_files"].as<vector<string>>();
  vector<int>    ohdus      = configFile["ohdus"]      .as<vector<int>>();

  cout<<"             " <<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"             " <<endl;
  cout<<"Input files: " <<endl;
  for (const auto& file : inputFiles) cout <<"  " <<file << endl;
  cout<<"             " <<endl;
  cout<<"OHDUs: " <<endl;
  for (const auto& o : ohdus) cout <<"  " <<o;
  cout<<"             " <<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"             " <<endl;
  cout<<"             " <<endl;

  // Histo
  TH1D* hN    = new TH1D("","Distribution of N per window",        50, 0, 3000);
  TH1D* h_k1  = new TH1D("","Distribution of k_{1} per window",    10, 0, 10);
  TH1D* h_k2  = new TH1D("","Distribution of k_{2} per window",    10, 0, 10);
  TH1D* h_k3  = new TH1D("","Distribution of k_{3} per window",    10, 0, 10);
  TH1D* hq_1e = new TH1D("","Distribution of q_{1} = #frac{k}{N}", 50, 0, .1);
  TH1D* hq_2e = new TH1D("","Distribution of q_{2} = #frac{k}{N}", 50, 0, .1);
  TH1D* hq_3e = new TH1D("","Distribution of q_{3} = #frac{k}{N}", 50, 0, .1);
  TH1D* hq_4e = new TH1D("","Distribution of q_{4} = #frac{k}{N}", 50, 0, .1);
  TH1D* hq_5e = new TH1D("","Distribution of q_{5} = #frac{k}{N}", 50, 0, .1);

  // Clock
  auto a = high_resolution_clock::now();

  vector<float> n1ElectronEventsVec = {};
  vector<float> n2ElectronEventsVec = {};
  vector<float> n3ElectronEventsVec = {};
  vector<float> n4ElectronEventsVec = {};

  vector<float> unmaskedPixels1eVec = {};
  vector<float> unmaskedPixels2eVec = {};
  vector<float> unmaskedPixels3eVec = {};
  vector<float> unmaskedPixels4eVec = {};

  vector<float> oneElectronRateVec   = {};
  vector<float> twoElectronRateVec   = {};
  vector<float> threeElectronRateVec = {};
  vector<float> fourElectronRateVec  = {};

  for (auto wp : thresholds) {
    cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Threshold = "<< wp <<endl;

    int nImages = 0;
    int nImagesInConfig = inputFiles.size();

    float n1ElectronEvents = 0;
    float n2ElectronEvents = 0;
    float n3ElectronEvents = 0;
    float n4ElectronEvents = 0;

    float unmaskedPixels1e = 0;
    float unmaskedPixels2e = 0;
    float unmaskedPixels3e = 0;
    float unmaskedPixels4e = 0;

    for (const auto& rf : inputFiles) {
      nImages+=1;

      cout<<"                                                    "<<endl;
      cout<<"  READING FILE: "<< nImages << "/"<< nImagesInConfig <<endl;
      cout<<"                "<< rf                               <<endl;

      for (int oh=0; oh<ohdus.size(); oh++){
        int ohdu = ohdus.at(oh);

        cout<<"                   *ohdu "<< ohdu << endl;

        vector<int> qs = analyze(wp, rf, ohdu, hq_1e, hq_2e, hq_3e, hq_4e, hN, h_k1, h_k2, h_k3);

        cout<<"                      - Number of 1 electron events " << qs.at(0) << endl;
        cout<<"                      - Number of 2 electron events " << qs.at(1) << endl;
        cout<<"                      - Number of 3 electron events " << qs.at(2) << endl;
        cout<<"                      - Number of 4 electron events " << qs.at(3) << endl;
        cout<<"                      - Number of unmasked pixels 1e "<< qs.at(4) << endl;
        cout<<"                      - Number of unmasked pixels 2e "<< qs.at(5) << endl;
        cout<<"                      - Number of unmasked pixels 3e "<< qs.at(6) << endl;
        cout<<"                      - Number of unmasked pixels 4e "<< qs.at(7) << endl;

        n1ElectronEvents += qs.at(0);
        n2ElectronEvents += qs.at(1);
        n3ElectronEvents += qs.at(2);
        n4ElectronEvents += qs.at(3);

        unmaskedPixels1e += qs.at(4);
        unmaskedPixels2e += qs.at(5);
        unmaskedPixels3e += qs.at(6);
        unmaskedPixels4e += qs.at(7);
      }
      //if (nImages>120) break;
    }

    n1ElectronEventsVec.push_back(n1ElectronEvents);
    n2ElectronEventsVec.push_back(n2ElectronEvents);
    n3ElectronEventsVec.push_back(n3ElectronEvents);
    n4ElectronEventsVec.push_back(n4ElectronEvents);

    unmaskedPixels1eVec.push_back(unmaskedPixels1e);
    unmaskedPixels2eVec.push_back(unmaskedPixels2e);
    unmaskedPixels3eVec.push_back(unmaskedPixels3e);
    unmaskedPixels4eVec.push_back(unmaskedPixels4e);

    oneElectronRateVec.push_back(   unmaskedPixels1e > 0 ? n1ElectronEvents/unmaskedPixels1e : 0.0f);
    twoElectronRateVec.push_back(   unmaskedPixels2e > 0 ? n2ElectronEvents/unmaskedPixels2e : 0.0f);
    threeElectronRateVec.push_back( unmaskedPixels3e > 0 ? n3ElectronEvents/unmaskedPixels3e : 0.0f);
    fourElectronRateVec.push_back(  unmaskedPixels4e > 0 ? n4ElectronEvents/unmaskedPixels4e : 0.0f);
  }

  cout << " "<<  endl;
  cout << " "<<  endl;
//  canvasMaker(hN,    "Distribution of N",     "./pdfs/hN.pdf");
//  canvasMaker(h_k1,  "Distribution of k_{1}", "./pdfs/k1.pdf");
//  canvasMaker(h_k2,  "Distribution of k_{2}", "./pdfs/k2.pdf");
//  canvasMaker(h_k3,  "Distribution of k_{3}", "./pdfs/k3.pdf");
//  canvasMaker(hq_1e, "q = #frac{k_{1}}{N}",   "./pdfs/q_statistic_1e.pdf");
//  canvasMaker(hq_2e, "q = #frac{k_{2}}{N}",   "./pdfs/q_statistic_2e.pdf");
//  canvasMaker(hq_3e, "q = #frac{k_{3}}{N}",   "./pdfs/q_statistic_3e.pdf");
//  canvasMaker(hq_4e, "q = #frac{k_{4}}{N}",   "./pdfs/q_statistic_4e.pdf");

  // chrono
  auto b = high_resolution_clock::now();
  cout << " "<<  endl;
  cout << " "<<  endl;
  cout << "================== "<<  endl;
  cout << "Took " << duration_cast<seconds>(b - a).count() << " seconds" <<  endl;
  cout << "================== "<<  endl;
  cout << " "<<  endl;
  cout << " "<<  endl;

  // output txt to plot with python
  std::ofstream outFile0("./txts/more_survElec.txt");
  std::ofstream outFile1("./txts/more_survPix.txt");
  std::ofstream outFile2("./txts/more_rates.txt");

  if (!outFile0) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile1) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile2) { cerr << "Error opening file for writing!" << endl; }

  for (size_t t = 0; t < thresholds.size(); ++t) {
    outFile0 << thresholds.at(t) << "\t"
             << n1ElectronEventsVec.at(t) << "\t"
             << n2ElectronEventsVec.at(t) << "\t"
             << n3ElectronEventsVec.at(t) << "\t"
             << n4ElectronEventsVec.at(t) << endl;

    outFile1 << thresholds.at(t) << "\t"
             << unmaskedPixels1eVec.at(t) << "\t"
             << unmaskedPixels2eVec.at(t) << "\t"
             << unmaskedPixels3eVec.at(t) << "\t"
             << unmaskedPixels4eVec.at(t) << endl;

    outFile2 << thresholds.at(t) << "\t"
             << oneElectronRateVec.at(t) << "\t"
             << twoElectronRateVec.at(t) << "\t"
             << threeElectronRateVec.at(t) << "\t"
             << fourElectronRateVec.at(t) << endl;
  }

  outFile0.close();
  outFile1.close();
  outFile2.close();

  return 0;
}
