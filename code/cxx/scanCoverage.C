#include "helperFunctions.C"
#include <iostream>
#include <vector>
#include <string>
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

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

// Paints, in red, every pixel touched by at least one of the sliding
// circular binomial windows (same n=60/r=30 grid as analyze() in
// multichannel.C), overlaid on the quadrant's event image. The coverage
// is purely geometric (it does not depend on masks or hit data), so gaps
// left uncovered between circles show up in black.
int scanCoverage(string fileName, int ohdu){

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

  auto [matrix_data, baseMask] = loadImage3(fileName, ohdu);

  int Nx  = 3200;
  int Ny  = 520;
  int n   = 60;
  float r = n/2;

  vector<vector<int>> coverage(Ny, vector<int>(Nx, 0));

  for (int xcoord = 0; xcoord <= Nx - n; xcoord = xcoord + n/2) {
    for (int ycoord = 0; ycoord <= Ny - n; ycoord = ycoord + n/2) {

      int x_center = xcoord + r;
      int y_center = ycoord + r;

      for (int i = ycoord; i < ycoord + n; ++i) {
        for (int j = xcoord; j < xcoord + n; ++j) {
          if (isInsideCircle(r, x_center, y_center, j, i)) {
            coverage.at(i).at(j) = 1;
          }
        }
      }
    }
  }

  int nScanned = 0;
  for (int i = 0; i < Ny; i++)
    for (int j = 0; j < Nx; j++)
      if (coverage.at(i).at(j)) nScanned++;

  cout << "Scanned pixels: " << nScanned << " / " << Nx*Ny
       << " (" << 100.0*nScanned/(Nx*Ny) << "%)" << endl;

  int lta = extractLTA(fileName);

  TH2F* h1 = new TH2F("h1", Form("LTA %d, OHDU %d;X;Y", lta, ohdu), Nx, 0, Nx, Ny, 0, Ny);
  TH2F* h2 = new TH2F("h2", "Window coverage",                     Nx, 0, Nx, Ny, 0, Ny);

  for (int ix = 0; ix < Nx; ++ix) {
    for (int iy = 0; iy < Ny; ++iy) {
      h1->SetBinContent(ix + 1, iy + 1, matrix_data.at(iy).at(ix));
      if (coverage.at(iy).at(ix)) h2->SetBinContent(ix + 1, iy + 1, 1);
    }
  }

  h1->SetStats(false);
  h2->SetStats(false);
  h2->SetFillColorAlpha(kRed, 0.5);
  h2->SetLineColorAlpha(kRed, 0.5);

  TCanvas* c = new TCanvas("c", "Scan coverage", 1600, 400);
  c->SetFillColor(kBlack);
  c->cd();
  h1->Draw("COLZ");
  h2->Draw("BOX SAME");

  gSystem->mkdir("./results/pdfs", true);
  string outName = "./results/pdfs/scanCoverage_lta" + to_string(lta) + "_ohdu" + to_string(ohdu) + ".png";
  c->SaveAs(outName.c_str());

  return 0;
}
