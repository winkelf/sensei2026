// How much the window-flag masks (BuildFlagMasks.C) remove, per dataset,
// quad selection and exposure, from the flagSummary trees of the flagged
// copies. Also lists, per dataset, a typical image and the most affected one
// (good quads) as candidates for FlagDisplay.C.
//
// usage (from the project folder): root -l -b -q FlagSummary.C

#include <yaml-cpp/yaml.h>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include <map>
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

using namespace std;

int FlagSummary(const char *flagDir="data/March_flagged"){

    struct DS{ string config, name; };
    vector<DS> datasets = {
        {"config.yaml",                               "Old"},
        {"results/hits_skp_2025-09-28/config.yaml",   "Sept"},
        {"results/hits_skp_2026-03-17_run8/config.yaml","March"},
        {"results/hits_skp_2026-07-03/config.yaml",   "July"},
        {"results/hits_skp_all/config.yaml",          "Inclusive"},
    };
    const char *flagName[3] = {"A (k1>=2)","B (k1>=3)","C (k1>=5)"};

    printf("%-10s %-8s %-7s %6s %9s | %-26s %-26s %-26s\n","dataset","quads","exp","images","windows",
           "flag A: windows, live px","flag B: windows, live px","flag C: windows, live px");

    for(auto &ds : datasets){

        YAML::Node cfg = YAML::LoadFile(ds.config);
        map<string,map<int,set<int>>> quads;
        for(const char *q : {"good_quads","notbad_quads"})
            for(auto it=cfg[q].begin(); it!=cfg[q].end(); ++it){
                auto v = it->second.as<vector<int>>();
                quads[q][it->first.as<int>()] = set<int>(v.begin(),v.end());
            }

        struct Cand{ double frac; string file; int run, lta, ohdu; };
        vector<Cand> cands;

        for(auto ex=cfg["exposures"].begin(); ex!=cfg["exposures"].end(); ++ex){

            string exp = ex->first.as<string>();
            auto files = ex->second.as<vector<string>>();

            for(const char *q : {"good_quads","notbad_quads"}){

                long nImg=0, nWin=0, fWin[3]={0,0,0};
                double live=0, fLive[3]={0,0,0};

                for(auto &path : files){
                    string f = string(flagDir)+"/"+gSystem->BaseName(path.c_str());
                    TFile F(f.c_str());
                    TTree *t = F.IsZombie() ? nullptr : (TTree*)F.Get("flagSummary");
                    if(!t){ cout << "missing flagSummary: " << f << endl; continue; }
                    int r,l,o,nw; Long64_t lp, fw[3], fl[3];
                    t->SetBranchAddress("runID",&r); t->SetBranchAddress("LTA",&l); t->SetBranchAddress("OHDU",&o);
                    t->SetBranchAddress("livePix",&lp); t->SetBranchAddress("nWindows",&nw);
                    t->SetBranchAddress("flagWindows",fw); t->SetBranchAddress("flagLivePix",fl);
                    for(Long64_t i=0;i<t->GetEntries();i++){
                        t->GetEntry(i);
                        auto qi = quads[q].find(l);
                        if(qi==quads[q].end() || !qi->second.count(o)) continue;
                        nImg++; nWin+=nw; live+=lp;
                        for(int k=0;k<3;k++){ fWin[k]+=fw[k]; fLive[k]+=fl[k]; }
                        if(string(q)=="good_quads" && lp>0 && fw[0]>0)
                            cands.push_back({(double)fl[0]/lp,f,r,l,o});
                    }
                }

                printf("%-10s %-8s %-7s %6ld %9ld |",ds.name.c_str(),string(q)=="good_quads"?"good":"not-bad",exp.c_str(),nImg,nWin);
                for(int k=0;k<3;k++) printf(" %7ld (%5.2f%%), %5.2f%% |",fWin[k],100.0*fWin[k]/std::max(1L,nWin),100.0*fLive[k]/std::max(1.0,live));
                printf("\n");
            }
        }

        if(!cands.empty()){
            sort(cands.begin(),cands.end(),[](const Cand&a,const Cand&b){ return a.frac<b.frac; });
            const Cand &med = cands[cands.size()/2], &mx = cands.back();
            printf("   %s display candidates (good quads, flag A live-pixel fraction): typical %.2f%% run %d LTA %d OHDU %d; most affected %.2f%% run %d LTA %d OHDU %d\n",
                   ds.name.c_str(),100*med.frac,med.run,med.lta,med.ohdu,100*mx.frac,mx.run,mx.lta,mx.ohdu);
            printf("     typical: %s\n     most:    %s\n",med.file.c_str(),mx.file.c_str());
        }
    }
    return 0;
}
