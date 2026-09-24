#include <yaml-cpp/yaml.h>
#include "helperFunctions.C"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
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

// Same standard SENSEI mask bits as qWeighted.C, plus the LEC bit on top.
// The LEC mask replaces qWeighted's dynamic binomial mask; everything else
// (image loading, event definition, base mask) is identical so the two
// jobs are directly comparable.
static constexpr int PREEXISTING_MASK_BITS = 0x56FC;
static constexpr int LEC_MASK_BIT          = 8192;

// Counts events/pixels of one quadrant that survive
//   [0..3] the base mask only        (n1, n2, n345, nPix)
//   [4..7] base mask + LEC mask      (n1, n2, n345, nPix)
// Events are placed at their barycenter pixel with e truncated to an
// integer (loadImage3), exactly as in qWeighted.C.
vector<int> analyzeLEC(string fileName, int ohdu) {

  auto [matrix_data, mask] = loadImage3(fileName, ohdu, PREEXISTING_MASK_BITS | LEC_MASK_BIT);

  int Nx = 3200;
  int Ny = 520;

  vector<int> out(8, 0);

  for (int x = 0; x < Ny; ++x) {
    for (int y = 0; y < Nx; ++y) {

      int maskInPixel = mask.at(x).at(y);

      // Standard SENSEI masks.
      if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;

      bool lecMasked = (maskInPixel & LEC_MASK_BIT) != 0;

      out[3] += 1;
      if (!lecMasked) out[7] += 1;

      int event = matrix_data.at(x).at(y);
      if (event < 1 || event > 5) continue;

      int idx = (event == 1) ? 0 : (event == 2) ? 1 : 2;
      out[idx] += 1;
      if (!lecMasked) out[4 + idx] += 1;
    }
  }

  return out;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////                Main           /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Config: ../configs/config.yaml (input_files)
//
// LEC is a fixed static mask, so there is no threshold scan: one pass over
// the GOOD_QUADS images. Output in ../txts/lec/:
//   LEC_perImage.txt  one row per (file, ohdu)
//   LEC_survElec.txt / LEC_survPix.txt / LEC_rates.txt  totals, same
//   column layout as qWeighted.C's (single row, threshold column = 0)
int LEC(){

  YAML::Node configFile = YAML::LoadFile("../configs/config.yaml");
  vector<string> inputFiles = configFile["input_files"].as<vector<string>>();

  string outDir = "../txts/lec";
  gSystem->mkdir(outDir.c_str(), true);

  std::ofstream outPerImage((outDir + "/LEC_perImage.txt").c_str());
  if (!outPerImage) { cerr << "Error opening file for writing!" << endl; }
  outPerImage << "#lta\tohdu\trunId"
              << "\tn1_base\tn2_base\tn345_base\tnPix_base"
              << "\tn1_LEC\tn2_LEC\tn345_LEC\tnPix_LEC" << endl;

  auto a = high_resolution_clock::now();

  vector<double> tot(8, 0.0);
  int nImages = 0;

  for (const auto& rf : inputFiles) {

    int lta = extractLTA(rf);
    auto ltaIt = GOOD_QUADS.find(lta);
    if (ltaIt == GOOD_QUADS.end()) continue;
    const vector<int>& fileOhdus = ltaIt->second;

    int runId = extractRunId(rf);
    cout<<"  READING FILE: "<< rf << endl;

    for (int ohdu : fileOhdus) {
      vector<int> c = analyzeLEC(rf, ohdu);
      nImages += 1;

      outPerImage << lta << "\t" << ohdu << "\t" << runId;
      for (int v : c) outPerImage << "\t" << v;
      outPerImage << endl;

      for (int k = 0; k < 8; ++k) tot[k] += c[k];

      cout<<"                   *ohdu "<< ohdu
          <<"  base: "<< c[0] <<"/"<< c[1] <<"/"<< c[2] <<" pix "<< c[3]
          <<"   LEC: "<< c[4] <<"/"<< c[5] <<"/"<< c[6] <<" pix "<< c[7] << endl;
    }
  }
  outPerImage.close();

  std::ofstream outFile0((outDir + "/LEC_survElec.txt").c_str());
  std::ofstream outFile1((outDir + "/LEC_survPix.txt").c_str());
  std::ofstream outFile2((outDir + "/LEC_rates.txt").c_str());

  if (!outFile0) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile1) { cerr << "Error opening file for writing!" << endl; }
  if (!outFile2) { cerr << "Error opening file for writing!" << endl; }

  outFile0 << 0 << "\t" << tot[4] << "\t" << tot[5] << "\t" << tot[6] << endl;
  outFile1 << 0 << "\t" << tot[7] << endl;
  outFile2 << 0 << "\t" << tot[4]/tot[7] << "\t" << tot[5]/tot[7] << "\t" << tot[6]/tot[7] << endl;

  outFile0.close();
  outFile1.close();
  outFile2.close();

  auto b = high_resolution_clock::now();
  cout << " "<<  endl;
  cout << "================== "<<  endl;
  cout << nImages << " images" << endl;
  cout << "Base mask:       n1 "<< tot[0] <<"  n2 "<< tot[1] <<"  n345 "<< tot[2] <<"  pix "<< tot[3] << endl;
  cout << "Base + LEC mask: n1 "<< tot[4] <<"  n2 "<< tot[5] <<"  n345 "<< tot[6] <<"  pix "<< tot[7] << endl;
  cout << "Took " << duration_cast<seconds>(b - a).count() << " seconds" <<  endl;
  cout << "================== "<<  endl;

  return 0;
}
