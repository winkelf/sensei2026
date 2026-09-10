#include <yaml-cpp/yaml.h>
#include "helperFunctions.C"
#include <iostream>
#include <vector>
#include <string>
#include "TH2F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"

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

static constexpr int BINOMIAL_MASK_BIT     = 65536;
static constexpr int PREEXISTING_MASK_BITS = 0x56FC;

// Animates, frame by frame over a sweep of q1 cut values, which sliding
// windows of one quadrant get masked as the q1 threshold rises, next to
// a reference q1 = k1/N distribution (bottom panel, read from a
// previously generated root file) with a vertical line marking the
// current cut. Same n=60/r=30 window grid and mask bookkeeping as
// analyze() in multichannel.C, restricted to channel 1.
//
// Top panel layering:
//   - event image (COLZ)
//   - pre-existing SENSEI mask (PREEXISTING_MASK_BITS), static, in red
//   - windows newly cut by the current q1 threshold, in violet
//
// The quadrant (input_file/ohdu) and the reference histogram source
// (q1_root_file/q1_hist_name) are read from ../config/config_gif.yaml.
int q1MaskGif(){

  YAML::Node configFile = YAML::LoadFile("../config/gif_cfg.yaml");
  string fileName    = configFile["input_file"]  .as<string>();
  int    ohdu        = configFile["ohdu"]        .as<int>();
  string q1RootFile  = configFile["q1_root_file"].as<string>();
  string q1HistName  = configFile["q1_hist_name"].as<string>();

  cout << "Quadrant: " << fileName << "  (ohdu " << ohdu << ")" << endl;

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
  gStyle->SetTitleFontSize(0.04);

  auto [matrix_data, baseMask] = loadImage3(fileName, ohdu);

  int Nx  = 3200;
  int Ny  = 520;
  int n   = 60;
  float r = n/2;

  // Reference q1 distribution: read once, detach from its file so it
  // survives the TFile::Close() below.
  TFile* qFile = TFile::Open(q1RootFile.c_str(), "READ");
  if (!qFile || qFile->IsZombie()) {
    cerr << "Could not open " << q1RootFile << endl;
    return 1;
  }
  TH1D* hQ1ref = (TH1D*)qFile->Get(q1HistName.c_str());
  if (!hQ1ref) {
    cerr << "Histogram " << q1HistName << " not found in " << q1RootFile << endl;
    return 1;
  }
  hQ1ref = (TH1D*)hQ1ref->Clone("hQ1ref");
  hQ1ref->SetDirectory(0);
  qFile->Close();

  double qMax = hQ1ref->GetXaxis()->GetXmax();

  int lta = extractLTA(fileName);

  // Same threshold grid used for the production runs in multichannel.C.
  vector<float> thresholds = {
    0.000f, 0.001f, 0.002f, 0.003f, 0.004f,
    0.005f, 0.006f, 0.007f, 0.008f, 0.009f
  };

  gSystem->mkdir("./gif", true);
  string outGif = "./gif/q1_mask_lta" + to_string(lta) + "_ohdu" + to_string(ohdu) + ".gif";
  gSystem->Unlink(outGif.c_str());

  // Pre-existing SENSEI mask overlay: same for every frame, so build it
  // once outside the threshold loop.
  TH2F* hPre = new TH2F("hPre", "pre-existing mask", Nx, 0, Nx, Ny, 0, Ny);
  for (int ix = 0; ix < Nx; ++ix) {
    for (int iy = 0; iy < Ny; ++iy) {
      if ((baseMask.at(iy).at(ix) & PREEXISTING_MASK_BITS) != 0) {
        hPre->SetBinContent(ix + 1, iy + 1, 1);
      }
    }
  }
  hPre->SetStats(false);
  hPre->SetFillColorAlpha(kRed, 0.6);
  hPre->SetLineColorAlpha(kRed, 0.6);

  TCanvas* c = new TCanvas("c", "q1 masking", 1600, 1000);
  c->SetFillColor(kBlack);

  TPad* padTop = new TPad("padTop", "padTop", 0, 0.50, 1, 1);
  TPad* padBot = new TPad("padBot", "padBot", 0, 0,    1, 0.50);
  padTop->SetFillColor(kBlack);
  padBot->SetFillColor(kBlack);
  padTop->SetBottomMargin(0.15);
  padBot->SetBottomMargin(0.15);
  padTop->Draw();
  padBot->Draw();

  for (size_t f = 0; f < thresholds.size(); ++f) {

    float threshold = thresholds.at(f);

    // Fresh per-frame copy: the binomial mask accumulates window by
    // window only within a single threshold pass, exactly like analyze().
    vector<vector<int>> mask = baseMask;
    vector<vector<int>> maskedCoverage(Ny, vector<int>(Nx, 0));

    for (int xcoord = 0; xcoord <= Nx - n; xcoord = xcoord + n/2) {
      for (int ycoord = 0; ycoord <= Ny - n; ycoord = ycoord + n/2) {

        int x_center = xcoord + r;
        int y_center = ycoord + r;

        int k1     = 0;
        int trials = 0;

        for (int i = ycoord; i < ycoord + n; ++i) {
          for (int j = xcoord; j < xcoord + n; ++j) {

            if (!isInsideCircle(r, x_center, y_center, j, i)) continue;

            int maskInPixel = mask.at(i).at(j);
            if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;
            if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;

            trials += 1;
            if (matrix_data.at(i).at(j) == 1) k1 += 1;
          }
        }

        if (trials == 0) continue;

        double q1 = static_cast<double>(k1) / static_cast<double>(trials);

        if (q1 > threshold) {
          for (int i = ycoord; i < ycoord + n; ++i) {
            for (int j = xcoord; j < xcoord + n; ++j) {
              if (!isInsideCircle(r, x_center, y_center, j, i)) continue;
              mask.at(i).at(j) |= BINOMIAL_MASK_BIT;
              maskedCoverage.at(i).at(j) = 1;
            }
          }
        }
      }
    }

    int nMasked = 0;
    for (int i = 0; i < Ny; i++)
      for (int j = 0; j < Nx; j++)
        if (maskedCoverage.at(i).at(j)) nMasked++;

    cout << "threshold=" << threshold
         << "  masked pixels=" << nMasked << "/" << Nx*Ny
         << " (" << 100.0*nMasked/(Nx*Ny) << "%)" << endl;

    ////////////////////////////////////////////////////////////////
    // Top panel: image + pre-existing mask (red) + q1 cut (violet)
    ////////////////////////////////////////////////////////////////

    TH2F* h1 = new TH2F(Form("h1_%zu", f), Form("LTA %d, OHDU %d  (q_{1} cut: %.3f);X;Y", lta, ohdu, threshold),
                         Nx, 0, Nx, Ny, 0, Ny);
    TH2F* h2 = new TH2F(Form("h2_%zu", f), "masked windows", Nx, 0, Nx, Ny, 0, Ny);

    for (int ix = 0; ix < Nx; ++ix) {
      for (int iy = 0; iy < Ny; ++iy) {
        h1->SetBinContent(ix + 1, iy + 1, matrix_data.at(iy).at(ix));
        if (maskedCoverage.at(iy).at(ix)) h2->SetBinContent(ix + 1, iy + 1, 1);
      }
    }

    h1->SetStats(false);
    h2->SetStats(false);
    h2->SetFillColorAlpha(kViolet, 0.6);
    h2->SetLineColorAlpha(kViolet, 0.6);

    padTop->cd();
    h1->Draw("COLZ");
    hPre->Draw("BOX SAME");
    h2->Draw("BOX SAME");

    TLegend* leg = new TLegend(0.70, 0.90, 0.98, 0.99);
    leg->SetFillColor(kBlack);
    leg->SetTextColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetNColumns(2);
    leg->AddEntry(hPre, "pre-existing mask", "f");
    leg->AddEntry(h2,   "q_{1} cut",         "f");
    leg->Draw();

    ////////////////////////////////////////////////////////////////
    // Bottom panel: reference q1 distribution + cut line
    ////////////////////////////////////////////////////////////////

    padBot->cd();
    padBot->SetLogy();

    TH1D* hQ1frame = (TH1D*)hQ1ref->Clone(Form("hQ1frame_%zu", f));
    hQ1frame->SetStats(false);
    hQ1frame->SetLineColor(kWhite);
    hQ1frame->GetXaxis()->SetTitle("q_{1} = k_{1}/N");
    hQ1frame->GetYaxis()->SetTitle("Number of entries");
    hQ1frame->Draw("HIST");

    double yMin = 0.5; // logy floor
    double yMax = hQ1frame->GetMaximum() * 2;

    TLine* line = new TLine(threshold, yMin, threshold, yMax);
    line->SetLineColor(kViolet);
    line->SetLineWidth(2);
    line->Draw("SAME");

    c->cd();
    c->Modified();
    c->Update();
    c->Print((outGif + "+70").c_str());

    delete leg;
    delete line;
    delete hQ1frame;
    delete h1;
    delete h2;
  }

  // ROOT's gif writer doesn't set the "loop forever" (NETSCAPE2.0) block,
  // so strict viewers play the animation once and stop. Add it in place.
  gSystem->Exec(("convert -loop 0 " + outGif + " " + outGif).c_str());

  cout << "Wrote " << outGif << endl;

  return 0;
}
