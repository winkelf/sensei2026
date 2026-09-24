#include <yaml-cpp/yaml.h>
#include "helperFunctions.C"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TSystem.h"

using namespace std;

static constexpr int BINOMIAL_MASK_BIT     = 65536;
static constexpr int PREEXISTING_MASK_BITS = 0x56FC;

// Exports everything the Python q/mask applet (qApplet.py) needs for one
// quadrant: the raw event image, the pre-existing SENSEI mask (as a
// pixel-level boolean grid), and a per-window table of the q1 = k1/N
// statistic used by the binomial mask.
//
// Same n=60/r=30 window grid and sequential mask bookkeeping as
// q1MaskGif.C: a window's trial count excludes both the pre-existing
// mask and any BINOMIAL_MASK_BIT set by an earlier (lower xcoord/ycoord)
// window at the same threshold, so results.masked here matches what the
// production binomial mask would flag at this threshold.
//
// Config: ../configs/qapplet_cfg.yaml (input_file, ohdu, threshold, output_dir)
int qAppletExport(){

  YAML::Node configFile = YAML::LoadFile("../configs/qapplet_cfg.yaml");
  string fileName  = configFile["input_file"] .as<string>();
  int    ohdu      = configFile["ohdu"]       .as<int>();
  float  threshold = configFile["threshold"]  .as<float>();
  string outDir    = configFile["output_dir"] .as<string>();

  cout << "Quadrant: " << fileName << "  (ohdu " << ohdu << ")" << endl;
  cout << "Threshold: " << threshold << endl;

  gSystem->mkdir(outDir.c_str(), true);

  auto [matrix_data, baseMask] = loadImage3(fileName, ohdu);

  int Nx  = 3200;
  int Ny  = 520;
  int n   = 60;
  float r = n/2;

  ////////////////////////////////////////////////////////////////////
  // Event image (Ny rows x Nx cols)
  ////////////////////////////////////////////////////////////////////

  ofstream imgOut(outDir + "/image.txt");
  for (int iy = 0; iy < Ny; ++iy) {
    for (int ix = 0; ix < Nx; ++ix) {
      if (ix > 0) imgOut << " ";
      imgOut << matrix_data.at(iy).at(ix);
    }
    imgOut << "\n";
  }
  imgOut.close();

  ////////////////////////////////////////////////////////////////////
  // Pre-existing SENSEI mask, as a 0/1 grid (Ny rows x Nx cols)
  ////////////////////////////////////////////////////////////////////

  ofstream maskOut(outDir + "/premask.txt");
  for (int iy = 0; iy < Ny; ++iy) {
    for (int ix = 0; ix < Nx; ++ix) {
      if (ix > 0) maskOut << " ";
      maskOut << (((baseMask.at(iy).at(ix) & PREEXISTING_MASK_BITS) != 0) ? 1 : 0);
    }
    maskOut << "\n";
  }
  maskOut.close();

  ////////////////////////////////////////////////////////////////////
  // Per-window q1 table
  ////////////////////////////////////////////////////////////////////

  vector<vector<int>> mask = baseMask;

  ofstream winOut(outDir + "/windows.csv");
  winOut << "xcoord,ycoord,xcenter,ycenter,radius,k1,k2,k3plus,trials,q1,masked\n";

  int nWindows = 0;
  int nMasked  = 0;

  for (int xcoord = 0; xcoord <= Nx - n; xcoord += n/2) {
    for (int ycoord = 0; ycoord <= Ny - n; ycoord += n/2) {

      int x_center = xcoord + r;
      int y_center = ycoord + r;

      int k1     = 0;
      int k2     = 0;
      int k3plus = 0;
      int trials = 0;

      for (int i = ycoord; i < ycoord + n; ++i) {
        for (int j = xcoord; j < xcoord + n; ++j) {

          if (!isInsideCircle(r, x_center, y_center, j, i)) continue;

          int maskInPixel = mask.at(i).at(j);
          if ((maskInPixel & PREEXISTING_MASK_BITS) != 0) continue;
          if ((maskInPixel & BINOMIAL_MASK_BIT) != 0) continue;

          trials += 1;
          int val = matrix_data.at(i).at(j);
          if (val == 1) k1 += 1;
          else if (val == 2) k2 += 1;
          else if (val >= 3) k3plus += 1;
        }
      }

      if (trials == 0) continue;

      double q1 = static_cast<double>(k1) / static_cast<double>(trials);
      int masked = (q1 > threshold) ? 1 : 0;

      if (masked) {
        for (int i = ycoord; i < ycoord + n; ++i) {
          for (int j = xcoord; j < xcoord + n; ++j) {
            if (!isInsideCircle(r, x_center, y_center, j, i)) continue;
            mask.at(i).at(j) |= BINOMIAL_MASK_BIT;
          }
        }
        nMasked += 1;
      }

      winOut << xcoord << "," << ycoord << "," << x_center << "," << y_center << ","
             << r << "," << k1 << "," << k2 << "," << k3plus << "," << trials << "," << q1 << "," << masked << "\n";

      nWindows += 1;
    }
  }
  winOut.close();

  cout << "Wrote " << outDir << "/{image.txt,premask.txt,windows.csv}" << endl;
  cout << nWindows << " windows, " << nMasked << " masked by q1 > " << threshold << endl;

  return 0;
}
