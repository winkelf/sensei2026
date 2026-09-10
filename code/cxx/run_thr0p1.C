#include "multichannel.C"

// Standalone re-run of the threshold = 0.1 diagnostic pass from multichannel(),
// reusing its analyze()/GOOD_QUADS/extractLTA, but with:
//  - a 4th k-histogram (k_{4+}) that analyze() already computes internally
//    but the original pass never plotted, and
//  - per-channel q binning (0-0.10 / 0-0.05 / 0-0.02 / 0-0.02) instead of a
//    single narrow 0-0.002 axis shared by every channel.
vector<int> analyzeK4(float threshold,
                       string fileName,
                       int ohdu,
                       TH1D* hq_1e, TH1D* hq_2e, TH1D* hq_3e, TH1D* hq_4e,
                       TH1D* hN,
                       TH1D* h_k1, TH1D* h_k2, TH1D* h_k3, TH1D* h_k4) {

  int n   = 60;
  float r = n/2;

  auto [matrix_data, baseMask] = loadImage3(fileName, ohdu);
  vector<vector<vector<int>>> previousMasks(N_CHANNELS, baseMask);

  int Nx = 3200;
  int Ny = 520;

  for (int xcoord = 0; xcoord <= Nx - n; xcoord = xcoord + n/2) {
    for (int ycoord = 0; ycoord <= Ny - n; ycoord = ycoord + n/2) {

      int k[N_CHANNELS]      = {0, 0, 0, 0};
      int trials[N_CHANNELS] = {0, 0, 0, 0};

      int x_center = xcoord + r;
      int y_center = ycoord + r;

      for (int i = ycoord; i < ycoord + n; ++i) {
        for (int j = xcoord; j < xcoord + n; ++j) {

          int eventInPixel = matrix_data.at(i).at(j);

          bool isInside = isInsideCircle(r, x_center, y_center, j, i);
          if (!isInside) continue;

          for (int channel = 0; channel < N_CHANNELS; ++channel) {

            int maskInPixel = int(previousMasks.at(channel).at(i).at(j));

            if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;
            if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;

            trials[channel] += 1;

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

      if (trials[0] > 0) hN->Fill(trials[0]);

      h_k1->Fill(k[0]);
      h_k2->Fill(k[1]);
      h_k3->Fill(k[2]);
      h_k4->Fill(k[3]);

      if (trials[0] > 0) hq_1e->Fill(P[0]);
      if (trials[1] > 0) hq_2e->Fill(P[1]);
      if (trials[2] > 0) hq_3e->Fill(P[2]);
      if (trials[3] > 0) hq_4e->Fill(P[3]);

      for (int channel = 0; channel < N_CHANNELS; ++channel) {
        if (trials[channel] == 0) continue;

        if (P[channel] > threshold) {
          for (int i = ycoord; i < ycoord + n; ++i) {
            for (int j = xcoord; j < xcoord + n; ++j) {
              bool isInside = isInsideCircle(r, x_center, y_center, j, i);
              if (!isInside) continue;
              previousMasks.at(channel).at(i).at(j) |= BINOMIAL_MASK_BIT;
            }
          }
        }
      }
    }
  }

  int nUnmaskedPixels[N_CHANNELS] = {0, 0, 0, 0};
  for (int channel = 0; channel < N_CHANNELS; ++channel) {
    nUnmaskedPixels[channel] = countUnmaskedPixels2(previousMasks.at(channel), true);
  }

  int nEventsUnmasked[N_CHANNELS] = {0, 0, 0, 0};

  for (int x = 0; x < 520; ++x) {
    for (int y = 0; y < 3200; ++y) {
      int event = matrix_data.at(x).at(y);
      for (int channel = 0; channel < N_CHANNELS; ++channel) {
        if (event != channel + 1) continue;
        int maskInPixel = previousMasks.at(channel).at(x).at(y);
        if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;
        if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;
        nEventsUnmasked[channel] += 1;
      }
    }
  }

  vector<int> rVector = {
    nEventsUnmasked[0], nEventsUnmasked[1], nEventsUnmasked[2], nEventsUnmasked[3],
    nUnmaskedPixels[0], nUnmaskedPixels[1], nUnmaskedPixels[2], nUnmaskedPixels[3]
  };

  return rVector;
}

int run_thr0p1(){

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

  YAML::Node configFile = YAML::LoadFile("../config/config.yaml");
  vector<string> inputFiles = configFile["input_files"].as<vector<string>>();

  cout<<"==================================================================================================================================================================================="<<endl;
  cout<<"Threshold = 0.1 (pdfs + root file only, no txt output, rebinned + k4)"<<endl;
  cout<<"==================================================================================================================================================================================="<<endl;

  float extraThreshold = 0.1f;

  string tag      = "thr0p1";
  string fullTitle = "threshold = 0.1";

  TH1D* hN  = new TH1D(("hN_" +tag).c_str(), ("N per window ("      +fullTitle+")").c_str(), 100, 0,    3000);
  TH1D* hK1 = new TH1D(("hK1_"+tag).c_str(), ("k_{1} per window ("  +fullTitle+")").c_str(), 11, -0.5, 10.5);
  TH1D* hK2 = new TH1D(("hK2_"+tag).c_str(), ("k_{2} per window ("  +fullTitle+")").c_str(), 11, -0.5, 10.5);
  TH1D* hK3 = new TH1D(("hK3_"+tag).c_str(), ("k_{3} per window ("  +fullTitle+")").c_str(), 11, -0.5, 10.5);
  TH1D* hK4 = new TH1D(("hK4_"+tag).c_str(), ("k_{4+} per window (" +fullTitle+")").c_str(), 11, -0.5, 10.5);
  // Same bin width (0.0005) across all four channels, rather than a
  // fixed bin count, so the q resolution is consistent regardless of
  // each channel's axis range.
  TH1D* hQ1 = new TH1D(("hQ1_"+tag).c_str(), ("q_{1} = k_{1}/N ("   +fullTitle+")").c_str(), 200, 0,    0.10);
  TH1D* hQ2 = new TH1D(("hQ2_"+tag).c_str(), ("q_{2} = k_{2}/N ("   +fullTitle+")").c_str(), 100, 0,    0.05);
  TH1D* hQ3 = new TH1D(("hQ3_"+tag).c_str(), ("q_{3} = k_{3}/N ("   +fullTitle+")").c_str(), 40,  0,    0.02);
  TH1D* hQ4 = new TH1D(("hQ4_"+tag).c_str(), ("q_{4+} = k_{4+}/N (" +fullTitle+")").c_str(), 40,  0,    0.02);

  for (const auto& rf : inputFiles) {

    int lta = extractLTA(rf);
    auto ltaIt = GOOD_QUADS.find(lta);
    if (ltaIt == GOOD_QUADS.end()) continue;
    const vector<int>& fileOhdus = ltaIt->second;

    cout<<"  READING FILE: "<< rf << endl;

    for (int oh=0; oh<fileOhdus.size(); oh++){
      int ohdu = fileOhdus.at(oh);

      cout<<"                   *ohdu "<< ohdu << endl;

      analyzeK4(extraThreshold, rf, ohdu,
                hQ1, hQ2, hQ3, hQ4,
                hN, hK1, hK2, hK3, hK4);
    }
  }

  canvasMaker(hN,  "Distribution of N",      "./results/pdfs/hN_thr0p1.pdf");
  canvasMaker(hK1, "Distribution of k_{1}",  "./results/pdfs/k1_thr0p1.pdf");
  canvasMaker(hK2, "Distribution of k_{2}",  "./results/pdfs/k2_thr0p1.pdf");
  canvasMaker(hK3, "Distribution of k_{3}",  "./results/pdfs/k3_thr0p1.pdf");
  canvasMaker(hK4, "Distribution of k_{4+}", "./results/pdfs/k4_thr0p1.pdf");
  canvasMaker(hQ1, "q = #frac{k_{1}}{N}",    "./results/pdfs/q_statistic_1e_thr0p1.pdf");
  canvasMaker(hQ2, "q = #frac{k_{2}}{N}",    "./results/pdfs/q_statistic_2e_thr0p1.pdf");
  canvasMaker(hQ3, "q = #frac{k_{3}}{N}",    "./results/pdfs/q_statistic_3e_thr0p1.pdf");
  canvasMaker(hQ4, "q = #frac{k_{4+}}{N}",   "./results/pdfs/q_statistic_4e_thr0p1.pdf");

  TFile* outRootFile = new TFile("./results/root/multichannel_thr0p1.root", "RECREATE");
  hN ->Write();
  hK1->Write();
  hK2->Write();
  hK3->Write();
  hK4->Write();
  hQ1->Write();
  hQ2->Write();
  hQ3->Write();
  hQ4->Write();
  outRootFile->Close();

  return 0;
}
