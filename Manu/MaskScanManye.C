#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <string>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <yaml-cpp/yaml.h>

using namespace std;

// key = (RUNID, LTA, OHDU) -- a single detector image, not just a detector
typedef tuple<int,int,int> DetKey;

// LTA -> set of allowed OHDU quadrants
typedef map<int,set<int>> QuadMap;

QuadMap loadQuadMap(const YAML::Node &node){

    QuadMap quads;

    for(auto it = node.begin(); it!=node.end(); ++it){

        int lta = it->first.as<int>();
        vector<int> ohdus = it->second.as<vector<int>>();

        quads[lta] = set<int>(ohdus.begin(),ohdus.end());
    }

    return quads;
}

bool quadAllowed(const QuadMap &quads, int lta, int ohdu){

    auto it = quads.find(lta);

    if(it==quads.end())
        return false;

    return it->second.count(ohdu)>0;
}

struct WindowInfo{

    int xCenter;
    int yCenter;

    int N;

    int k1;
    int k2;
    int k3;
    int k4;

    WindowInfo(){

        N = 0;

        k1 = 0;
        k2 = 0;
        k3 = 0;
        k4 = 0;
    }

};

// full-grid dimensions used for the packed bitmask
static const int NX = 3200;
static const int NY = 520;

inline int pixIdx(int x,int y){
    return x*NY + y;
}

// returns (creating/resizing if needed) the packed mask for a given key
vector<bool>& getMask(map<DetKey,vector<bool>> &m, const DetKey &k){

    auto &v = m[k];

    if((int)v.size() != NX*NY)
        v.assign(NX*NY,false);

    return v;
}

int windowStudy(const char *fileName, int nImages, int radius,
                 const QuadMap &quadMap, const char *label){

    gROOT->SetBatch(kTRUE);

    //////////////////////////////////////////////////////////
    // masks
    //////////////////////////////////////////////////////////

    // combined bad-pixel mask (excludes the "Neighbor" and
    // "Cluster shape" bits, which get set on multi-electron
    // hits themselves rather than on defective pixels)
    const int maskBits = 0x56FC;

    //////////////////////////////////////////////////////////
    // input
    //////////////////////////////////////////////////////////

    //const char *fileName =
    //    "data/hadded_EXPOSURE7200_Integrated.root";
        //"data/hits_corr_proc_skp_sensei_2026-03-17_run8_135K_sciencerun_commissioning_NROW520_NBINROW1_NCOL3200_NBINCOL1_EXPOSURE72000_11_519.root";

    TFile *file = TFile::Open(fileName);

    TTree *hitTree =
        (TTree*)file->Get("hitSumm");

    TTree *calPixTree =
        (TTree*)file->Get("calPixTree");

    //////////////////////////////////////////////////////////
    // hit tree branches
    //////////////////////////////////////////////////////////

    const int maxPix = 100000;

    int runID;
    int lta;
    int ohdu;
    int flag;

    float e;
    float xBary;
    float yBary;

    int nSavedPix;

    int xPix[maxPix];
    int yPix[maxPix];
    float ePix[maxPix];

    hitTree->SetBranchAddress("runID",&runID);
    hitTree->SetBranchAddress("LTANAME",&lta);
    hitTree->SetBranchAddress("ohdu",&ohdu);
    hitTree->SetBranchAddress("flag",&flag);

    hitTree->SetBranchAddress("e",&e);

    hitTree->SetBranchAddress("xBary",&xBary);
    hitTree->SetBranchAddress("yBary",&yBary);

    hitTree->SetBranchAddress("nSavedPix",&nSavedPix);

    hitTree->SetBranchAddress("xPix",xPix);
    hitTree->SetBranchAddress("yPix",yPix);
    hitTree->SetBranchAddress("ePix",ePix);

    //////////////////////////////////////////////////////////
    // calPixTree branches
    //////////////////////////////////////////////////////////

    int cal_x;
    int cal_y;
    double cal_e;
    int cal_mask;
    int cal_lta;
    int cal_ohdu;
    int cal_runID;

    calPixTree->SetBranchAddress("x",&cal_x);
    calPixTree->SetBranchAddress("y",&cal_y);
    calPixTree->SetBranchAddress("ePix",&cal_e);
    calPixTree->SetBranchAddress("mask",&cal_mask);
    calPixTree->SetBranchAddress("LTANAME",&cal_lta);
    calPixTree->SetBranchAddress("ohdu",&cal_ohdu);
    calPixTree->SetBranchAddress("RUNID",&cal_runID);

    //////////////////////////////////////////////////////////
    // discover detector *images*: (runID, LTA, OHDU)
    //////////////////////////////////////////////////////////

    set<DetKey> detectors;

    Long64_t nEvents = hitTree->GetEntries();

    for(Long64_t i=0;i<nEvents;i++){

        hitTree->GetEntry(i);

        if(!quadAllowed(quadMap,lta,ohdu))
            continue;

        detectors.insert(make_tuple(runID,lta,ohdu));
    }

    //////////////////////////////////////////////////////////
    // build masked-pixel lookup, per (runID, LTA, OHDU)
    // stored as a packed bitmask (not a set<pair>) to keep
    // memory usage bounded when there are many runs
    //////////////////////////////////////////////////////////

    map<DetKey,vector<bool>> maskedPixels;

    Long64_t nPixels = calPixTree->GetEntries();

    cout << "Building mask lookup..." << endl;

    for(Long64_t i=0;i<nPixels;i++){

        calPixTree->GetEntry(i);

        if(cal_x<0 || cal_x>=NX) continue;
        if(cal_y<0 || cal_y>=NY) continue;

        if((cal_mask & maskBits)==0) continue;

        DetKey key = make_tuple(cal_runID,cal_lta,cal_ohdu);

        auto &detMask = getMask(maskedPixels,key);

        detMask[pixIdx(cal_x,cal_y)] = true;
    }

    cout << "Done." << endl;

    //////////////////////////////////////////////////////////
    // build windows (geometric template, per run/detector image)
    //////////////////////////////////////////////////////////

    map<DetKey, vector<WindowInfo>> windows;

    for(auto &det : detectors){

        // ensures a properly-sized (all-false) mask exists even if
        // this run/detector had no entries in calPixTree
        auto &mask = getMask(maskedPixels,det);

        for(int xc=0+radius; xc<NX-radius; xc+=radius){

            for(int yc=0+radius; yc<NY-radius; yc+=radius){

                WindowInfo w;

                w.xCenter = xc;
                w.yCenter = yc;

                int N=0;

                for(int dx=-radius;dx<=radius;dx++){

                    for(int dy=-radius;dy<=radius;dy++){

                        if(dx*dx+dy*dy>radius*radius)
                            continue;

                        int xx=xc+dx;
                        int yy=yc+dy;

                        if(xx<0 || xx>=NX)
                            continue;

                        if(yy<0 || yy>=NY)
                            continue;

                        if(mask[pixIdx(xx,yy)]){
                            continue;
			}

                        N++;
                    }
                }

                w.N=N;

                windows[det].push_back(w);

            }
        }
    }

    cout << "Total detector-image groups: "
         << windows.size()
         << endl;

    //////////////////////////////////////////////////////////
    // output file -- one subfolder per quad category, then one
    // per radius, so different radii never overwrite each other
    //////////////////////////////////////////////////////////

    string outDir = string(label)+"/r"+to_string(radius);

    gSystem->mkdir(outDir.c_str(),kTRUE);

    string outName =
        outDir+"/windowStatisticsManye_"+label+"_"+to_string(nImages)+".root";

    TFile *out =
        new TFile(outName.c_str(),"RECREATE");

    TTree *tree =
        new TTree("windows","window statistics");

    int outRunID;

    int outLTA;
    int outOHDU;

    int outX;
    int outY;

    int outN;
    int outRadius = radius;

    int outK1;
    int outK2;
    int outK3;
    int outK4;

    double outR1;
    double outR2;
    double outR3;
    double outR4;

    tree->Branch("runID",&outRunID,"runID/I");
    tree->Branch("radius",&outRadius,"radius/I");

    tree->Branch("LTA",&outLTA,"LTA/I");
    tree->Branch("OHDU",&outOHDU,"OHDU/I");

    tree->Branch("x",&outX,"x/I");
    tree->Branch("y",&outY,"y/I");

    tree->Branch("N",&outN,"N/I");

    tree->Branch("k1",&outK1,"k1/I");
    tree->Branch("k2",&outK2,"k2/I");
    tree->Branch("k3",&outK3,"k3/I");
    tree->Branch("k4",&outK4,"k4/I");

    tree->Branch("r1",&outR1,"r1/D");
    tree->Branch("r2",&outR2,"r2/D");
    tree->Branch("r3",&outR3,"r3/D");
    tree->Branch("r4",&outR4,"r4/D");

    //////////////////////////////////////////////////////////
    // Histograms
    //////////////////////////////////////////////////////////

    TH1D hR1("hR1","k_{1}/N;k_{1}/N;Windows",100,0,0.10);
    TH1D hR2("hR2","k_{2}/N;k_{2}/N;Windows",100,0,0.05);
    TH1D hR3("hR3","k_{3}/N;k_{3}/N;Windows",100,0,0.02);
    TH1D hR4("hR4","k_{4+}/N;k_{4+}/N;Windows",100,0,0.02);

    //////////////////////////////////////////////////////////
    // Loop over hitSumm, counting each (runID, LTA, OHDU) image
    // separately, skipping events whose barycenter lands on a
    // masked pixel for THAT image
    //////////////////////////////////////////////////////////

    // counts are accumulated directly into "windows" -- no second
    // copy of the template is kept, which was doubling memory use

    cout << "Looping over hitSumm..." << endl;

    for(Long64_t i=0;i<nEvents;i++){

        if(i%100000==0)
            cout << i << " / " << nEvents << endl;

        hitTree->GetEntry(i);

        DetKey detKey = make_tuple(runID,lta,ohdu);

        auto tmplIt = windows.find(detKey);

        if(tmplIt==windows.end())
            continue;

        // mask-aware check: skip events whose barycenter falls
        // on a masked pixel for this specific run/detector image,
        // consistent with how N was built

        auto &mask = getMask(maskedPixels,detKey);

        int ix = (int)round(xBary);
        int iy = (int)round(yBary);

        if(ix<0 || ix>=NX || iy<0 || iy>=NY)
            continue;

        if(mask[pixIdx(ix,iy)])
            continue;

        auto &vec = tmplIt->second;

        for(auto &w : vec){

            double dx = xBary - w.xCenter;
            double dy = yBary - w.yCenter;

            if(dx*dx + dy*dy > radius*radius)
                continue;

            int ne = (int)e;

            if(ne==1)
                w.k1++;

            else if(ne==2)
                w.k2++;

            else if(ne==3)
                w.k3++;

            else if(ne>=4)
                w.k4++;
        }
    }

    cout << "Done counting events." << endl;

    // done reading from the input file -- close it now so its
    // TTree baskets don't stay cached in memory for the rest of
    // the process (this was leaking ~1 open file's worth of
    // basket buffers per call, across every file processed)
    file->Close();
    delete file;

    //////////////////////////////////////////////////////////
    // Fill output tree (one entry per run/detector-image/window)
    //////////////////////////////////////////////////////////

    for(auto &entry : windows){

        outRunID = get<0>(entry.first);
        outLTA   = get<1>(entry.first);
        outOHDU  = get<2>(entry.first);

        for(auto &w : entry.second){

            outX = w.xCenter;
            outY = w.yCenter;

            outN = w.N;

            outK1 = w.k1;
            outK2 = w.k2;
            outK3 = w.k3;
            outK4 = w.k4;

            if(w.N>0){

                outR1 = (double)w.k1/w.N;
                outR2 = (double)w.k2/w.N;
                outR3 = (double)w.k3/w.N;
                outR4 = (double)w.k4/w.N;

            }
            else continue;


            tree->Fill();

            hR1.Fill(outR1);
            hR2.Fill(outR2);
            hR3.Fill(outR3);
            hR4.Fill(outR4);

        }
    }

    //////////////////////////////////////////////////////////
    // Write everything
    //////////////////////////////////////////////////////////

    out->cd();

    tree->Write();

    hR1.Write();
    hR2.Write();
    hR3.Write();
    hR4.Write();

    out->Close();
    delete out;

    cout << "Output written to " << outName << endl;

    return 0;
}

int MaskScanManye(){

    YAML::Node configFile = YAML::LoadFile("config.yaml");
    vector<string> inputFiles = configFile["input_files"].as<vector<string>>();
    vector<int> radii = configFile["radii"].as<vector<int>>();

    QuadMap goodQuads   = loadQuadMap(configFile["good_quads"]);
    QuadMap notBadQuads = loadQuadMap(configFile["notbad_quads"]);

    struct Category{
        string label;
        QuadMap quads;
    };

    vector<Category> categories = {
        { "goodquads",   goodQuads   },
        { "notbadquads", notBadQuads }
    };

    for(auto &cat : categories){

        for(int radius : radii){

            cout << "=== " << cat.label << ", radius=" << radius << " ===" << endl;

            int nImages = 0;

            for (const auto& rf : inputFiles) {
              nImages+=1;
              int qs = windowStudy(rf.c_str(),nImages,radius,cat.quads,cat.label.c_str());
            }
        }
    }

    cout << "Done. Run PlotWindowStats.C to produce the summary plots." << endl;

    return 0;
}
