#include <yaml-cpp/yaml.h>
#include "helperFunctions.C"
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

vector<int> analyze(vector<int> maskList, float threshold, string fileName, int ohdu, TH1D* hq) {

  int n   = 30;
  float r = n/2;
  
  ////////////////////////////////////////////////////////////////////
  ////////////////////// Get previous masks (only for data) //////////
  ////////////////////////////////////////////////////////////////////

  auto [matrix_data, previousMask] = loadImage2(fileName, ohdu, maskList);

  // count unmasked pixels
  int nUnmaskedPixels = countUnmaskedPixels(previousMask, maskList, true);
  //cout << "  Number of unmasked pixels: "<< nUnmaskedPixels << endl;

  /////////////////////////////////////////////////////////////////////
  ///////// Counting n-electron events that survived mask /////////////
  /////////////////////////////////////////////////////////////////////

  //exit(1);

  int n1eEventsUnmasked = 0;
  int n2eEventsUnmasked = 0;
  int n3eEventsUnmasked = 0;

  for (int x = 0; x < 520; ++x) {
    for (int y = 0; y < 3200; ++y) {
      int event       = matrix_data.at(x).at(y);
      int maskInPixel = previousMask.at(x).at(y);
      //if (maskInPixel==8194) continue;

      int isMasked = 0;

      for (int m=0; m<maskList.size(); m++) {
        if ((maskInPixel&maskList.at(m))==maskList.at(m) ) {
          isMasked += 1;
        }
      }

      if (isMasked!=0){ continue; }

      if (event==1) n1eEventsUnmasked+=1;
      if (event==2) n2eEventsUnmasked+=1;
      if (event==3) n3eEventsUnmasked+=1;
    }
  }

  //exit(1);
  
  vector<int> rVector = { n1eEventsUnmasked, n2eEventsUnmasked, n3eEventsUnmasked, nUnmaskedPixels};

  return rVector;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////                Main           /////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int LEC(){

  // Make list of thresholds/cut values to loop over	
  //vector<float> thresholds = {0.1};
  vector<float> thresholds = {};
  //for (float i=0; i<0.01; i+=0.0005) thresholds.push_back(static_cast<float>(i));
  //for (float i=0; i<0.1; i+=0.005) thresholds.push_back(static_cast<float>(i));
  //for (float i=0; i<0.05; i+=0.0005) thresholds.push_back(static_cast<float>(i));
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
  vector<int>    ohdus      = configFile["ohdus"].as<vector<int>>();
  vector<int>    maskList   = configFile["mask_list_LEC"].as<vector<int>>();
 
  cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
  cout<<"             " <<endl; 
  cout<<"Input files: " <<endl; 
  for (const auto& file : inputFiles) cout <<"  " <<file << endl;
  cout<<"             " <<endl; 
  cout<<"OHDUs: " <<endl; 
  for (const auto& o : ohdus) cout <<"  " <<o;
  cout<<"             " <<endl; 

  // Histo
  TH1D* hq = new TH1D("","Distribution of q = #frac{k}{N}",50,0,.1);
 
  // Clock 
  auto a = high_resolution_clock::now();

  vector<float> n1ElectronEventsVec = {};
  vector<float> n2ElectronEventsVec = {};
  vector<float> n3ElectronEventsVec = {};
  vector<float> unmaskedPixelsVec   = {};
  vector<float> oneElectronRateVec  = {};
  vector<float> twoElectronRateVec  = {};
  vector<float> threeElectronRateVec  = {};
  
  for (auto wp : thresholds) {
    cout<<"-------------------------------------------------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Threshold = "<< wp <<endl;
  
    int nImages = 0;
    int nImagesInConfig = inputFiles.size();

    float n1ElectronEvents = 0;
    float n2ElectronEvents = 0;
    float n3ElectronEvents = 0;
    float unmaskedPixels   = 0;


    for (const auto& rf : inputFiles) {
      nImages+=1;

      cout<<"                                                    "<<endl;
      cout<<"  READING FILE: "<< nImages << "/"<< nImagesInConfig <<endl;
      cout<<"                "<< rf                               <<endl;

      for (int oh=0; oh<ohdus.size(); oh++){    
        int ohdu = ohdus.at(oh);

        cout<<"                   *ohdu "<< ohdu << endl;

        vector<int> qs  = analyze(maskList, wp, rf, ohdu, hq);
        
        cout<<"                      - Number of 1 electron events "<< qs.at(0) << endl;
        cout<<"                      - Number of 2 electron events "<< qs.at(1) << endl;
        cout<<"                      - Number of unmasked pixels   "<< qs.at(2) << endl;
        n1ElectronEvents += qs.at(0);
        n2ElectronEvents += qs.at(1);
        n2ElectronEvents += qs.at(2);
        unmaskedPixels   += qs.at(3);
      }	      
    //  if (nImages>120) break;
    }

    cout<<"                      - Number of 1 electron events "<< n1ElectronEvents << endl;
    cout<<"                      - Number of 2 electron events "<< n2ElectronEvents << endl;
    cout<<"                      - Number of unmasked pixels   "<< unmaskedPixels   << endl;
    //exit(1);

    n1ElectronEventsVec.push_back( n1ElectronEvents                );
    n2ElectronEventsVec.push_back( n2ElectronEvents                );
    n3ElectronEventsVec.push_back( n3ElectronEvents                );
    unmaskedPixelsVec  .push_back( unmaskedPixels                  );
    oneElectronRateVec .push_back( n1ElectronEvents/unmaskedPixels );
    twoElectronRateVec .push_back( n2ElectronEvents/unmaskedPixels );
    threeElectronRateVec .push_back( n3ElectronEvents/unmaskedPixels );
  }


  // output txt to plot with python
  std::ofstream outFile0("./txts/LEC_survElec.txt");
  std::ofstream outFile1("./txts/LEC_survPix.txt");
  std::ofstream outFile2("./txts/LEC_rates.txt");

  if (!outFile0) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile1) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile2) { cerr << "Error opening file for writing!" << endl; }


  for (size_t t = 0; t < thresholds.size(); ++t) {
    outFile0 << thresholds.at(t) <<"\t" << n1ElectronEventsVec.at(t) << "\t" << n2ElectronEventsVec.at(t) << "\t" << n3ElectronEventsVec.at(t) <<endl;
    outFile1 << thresholds.at(t) <<"\t" << unmaskedPixelsVec.  at(t) << "\t" << unmaskedPixelsVec.  at(t) << "\t" << unmaskedPixelsVec.  at(t) <<endl;
    outFile2 << thresholds.at(t) <<"\t" << oneElectronRateVec. at(t) << "\t" << twoElectronRateVec. at(t) << "\t" << threeElectronRateVec. at(t) <<endl;
  }

  outFile0.close();
  outFile1.close();
  outFile2.close();

  //TCanvas* c0 = new TCanvas();

  //c0->SetFillColor(kBlack);
  //c0->SetBatch(kTRUE);
  //hq->SetStats(0);
  //hq->GetXaxis()->SetTitle("q = #frac{k}{N}");
  //hq->GetYaxis()->SetTitle("Number of entries");
  //hq->GetYaxis()->SetAxisColor(kWhite);
  //hq->GetXaxis()->SetAxisColor(kWhite);
  ////hq->Scale(1/hq->Integral("width"));
  //hq->Draw("HIST");
  //c0->SetLogy();
  //c0->SaveAs("q_statistic.pdf");

  // chrono
  auto b = high_resolution_clock::now();
  cout << " "<<  endl;
  cout << " "<<  endl;
  cout << "================== "<<  endl;
  cout << "Took " << duration_cast<seconds>(b - a).count() << " seconds" <<  endl;
  cout << "================== "<<  endl;
  cout << " "<<  endl;
  cout << " "<<  endl;

  return 0;
}
