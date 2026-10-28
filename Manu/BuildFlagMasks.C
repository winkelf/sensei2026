// Window-flag masks written into copies of the input files.
//
// For every detector image (runID, LTA, OHDU) the r = 30 window grid is
// rebuilt exactly as in MaskScanManye.C (standard mask 0x56FC, hits whose
// rounded barycentre lies on a masked pixel are vetoed, counts by (int)e).
// A window is flagged by
//   flag A (bit 20, 0x100000): k2 >= 2, or k2 == 1 and k1 >= 2
//   flag B (bit 21, 0x200000): k2 >= 2, or k2 == 1 and k1 >= 3
//   flag C (bit 22, 0x400000): k2 >= 2, or k2 == 1 and k1 >= 5
// and every pixel inside a flagged window's circle gets that flag's bit ORed
// into calPixTree's "mask". The flags are nested (C implies B implies A).
// The bits are outside every bit used by the original masks; a file whose
// original masks already use any of them is not written.
//
// hitSumm keeps its "flag" (the OR of the original masks of the event's
// pixels) and gets three new branches, one per flag:
//   flagA = flag | 0x100000 if any pixel of the event carries flag A
//   flagB = flag | 0x200000 if any pixel of the event carries flag B
//   flagC = flag | 0x400000 if any pixel of the event carries flag C
//
// The input file is never modified: a full copy is written to the output
// folder under the same name, with every tree copied unchanged except
// calPixTree (mask updated) and hitSumm (flagA/flagB/flagC added), plus a
// new "flagWindows" tree listing the flagged windows and a "flagSummary"
// tree with per-image counts.
//
// usage (from the project folder):
//   root -l -b -q 'BuildFlagMasks.C("data/March","data/March_flagged")'
//   root -l -b -q 'BuildFlagMasks.C("data/March","data/March_flagged",0,2)'  // part 0 of 2
// Existing output files are skipped, so an interrupted run can be restarted.

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "TROOT.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TNamed.h"

using namespace std;

namespace flagmask {

const int NX = 3200, NY = 520;
const int RADIUS = 30;
const int STD_MASK = 0x56FC;

const int NFLAG = 3;
const int FLAG_BIT[NFLAG]  = { 0x100000, 0x200000, 0x400000 };
const int FLAG_K1[NFLAG]   = { 2, 3, 5 };          // k2 == 1 needs k1 >= this
const char *FLAG_NAME[NFLAG] = { "A", "B", "C" };
const int ALL_FLAGS = 0x700000;

typedef tuple<int,int,int> Img;   // (runID, LTA, OHDU)

inline int pix(int x, int y){ return x*NY+y; }

int windowFlags(int k1, int k2){
    int f = 0;
    for(int i=0;i<NFLAG;i++)
        if(k2>=2 || (k2==1 && k1>=FLAG_K1[i])) f |= FLAG_BIT[i];
    return f;
}

}

using namespace flagmask;

int buildFlagMask(const string &inPath, const string &outPath){

    if(!gSystem->AccessPathName(outPath.c_str())){
        cout << "exists, skipped: " << outPath << endl;
        return 0;
    }

    TFile *in = TFile::Open(inPath.c_str());
    if(!in || in->IsZombie()){ cout << "cannot open " << inPath << endl; return 1; }

    TTree *hit = (TTree*)in->Get("hitSumm");
    TTree *cal = (TTree*)in->Get("calPixTree");
    if(!hit || !cal){ cout << "missing hitSumm/calPixTree in " << inPath << endl; return 1; }

    // ---- standard mask per image, and a check that the flag bits are free
    int cx, cy, cm, cl, co, cr;
    cal->SetBranchStatus("*",0);
    for(const char *b : {"x","y","mask","LTANAME","ohdu","RUNID"}) cal->SetBranchStatus(b,1);
    cal->SetBranchAddress("x",&cx);       cal->SetBranchAddress("y",&cy);
    cal->SetBranchAddress("mask",&cm);    cal->SetBranchAddress("LTANAME",&cl);
    cal->SetBranchAddress("ohdu",&co);    cal->SetBranchAddress("RUNID",&cr);

    map<Img,vector<bool>> stdMask;
    unsigned int orMask = 0;
    Long64_t nCal = cal->GetEntries();
    for(Long64_t i=0;i<nCal;i++){
        cal->GetEntry(i);
        orMask |= (unsigned int)cm;
        if(cx<0||cx>=NX||cy<0||cy>=NY) continue;
        auto &m = stdMask[Img(cr,cl,co)];
        if(m.empty()) m.assign(NX*NY,false);
        if(cm & STD_MASK) m[pix(cx,cy)] = true;
    }
    if(orMask & ALL_FLAGS){
        cout << "ERROR: original masks already use flag bits (OR = 0x" << hex << orMask << dec
             << ") in " << inPath << " -- not written" << endl;
        in->Close();
        return 1;
    }

    // ---- hits, vetoed exactly as in MaskScanManye.C
    int hr, hl, ho; float he, hx, hy;
    hit->SetBranchStatus("*",0);
    for(const char *b : {"runID","LTANAME","ohdu","e","xBary","yBary"}) hit->SetBranchStatus(b,1);
    hit->SetBranchAddress("runID",&hr);  hit->SetBranchAddress("LTANAME",&hl); hit->SetBranchAddress("ohdu",&ho);
    hit->SetBranchAddress("e",&he);      hit->SetBranchAddress("xBary",&hx);   hit->SetBranchAddress("yBary",&hy);

    struct H{ float x, y; int ne; };
    map<Img,vector<H>> hits;
    for(Long64_t i=0;i<hit->GetEntries();i++){
        hit->GetEntry(i);
        Img key(hr,hl,ho);
        hits[key];                                    // images are discovered from hitSumm
        int ix=(int)round(hx), iy=(int)round(hy);
        if(ix<0||ix>=NX||iy<0||iy>=NY) continue;
        auto mit = stdMask.find(key);
        if(mit!=stdMask.end() && mit->second[pix(ix,iy)]) continue;
        hits[key].push_back({hx,hy,(int)he});
    }

    // ---- windows, flags and per-pixel flag bits
    vector<int> offX, offY;
    for(int dx=-RADIUS;dx<=RADIUS;dx++) for(int dy=-RADIUS;dy<=RADIUS;dy++)
        if(dx*dx+dy*dy<=RADIUS*RADIUS){ offX.push_back(dx); offY.push_back(dy); }

    map<Img,vector<uint8_t>> flagPix;   // bit i set = flag i

    // written under a temporary name and renamed when complete, so an
    // interrupted run never leaves a partial file that looks finished
    string tmpPath = outPath+".part";
    TFile *out = new TFile(tmpPath.c_str(),"RECREATE","",in->GetCompressionSettings());

    TTree *fw = new TTree("flagWindows","r = 30 windows flagged by A/B/C (flags = OR of mask bits)");
    int wRun, wLTA, wOHDU, wX, wY, wN, wK1, wK2, wFlags;
    fw->Branch("runID",&wRun,"runID/I");   fw->Branch("LTA",&wLTA,"LTA/I");   fw->Branch("OHDU",&wOHDU,"OHDU/I");
    fw->Branch("x",&wX,"x/I");             fw->Branch("y",&wY,"y/I");         fw->Branch("N",&wN,"N/I");
    fw->Branch("k1",&wK1,"k1/I");          fw->Branch("k2",&wK2,"k2/I");      fw->Branch("flags",&wFlags,"flags/I");

    TTree *fs = new TTree("flagSummary","per image: windows and pixels flagged by A/B/C");
    int sRun, sLTA, sOHDU, sNWin;
    Long64_t sLive, sFlagWin[NFLAG], sFlagPix[NFLAG], sFlagLive[NFLAG];
    fs->Branch("runID",&sRun,"runID/I");  fs->Branch("LTA",&sLTA,"LTA/I");  fs->Branch("OHDU",&sOHDU,"OHDU/I");
    fs->Branch("livePix",&sLive,"livePix/L");  fs->Branch("nWindows",&sNWin,"nWindows/I");
    fs->Branch("flagWindows",sFlagWin,"flagWindows[3]/L");
    fs->Branch("flagPix",sFlagPix,"flagPix[3]/L");
    fs->Branch("flagLivePix",sFlagLive,"flagLivePix[3]/L");

    for(auto &kv : hits){

        const Img &key = kv.first;
        auto mit = stdMask.find(key);
        static vector<bool> noMask(NX*NY,false);
        const vector<bool> &m = mit!=stdMask.end() ? mit->second : noMask;

        auto &fp = flagPix[key];
        fp.assign(NX*NY,0);

        sRun = get<0>(key); sLTA = get<1>(key); sOHDU = get<2>(key);
        sNWin = 0;
        for(int i=0;i<NFLAG;i++) sFlagWin[i]=sFlagPix[i]=sFlagLive[i]=0;
        sLive = 0;
        for(bool b : m) if(!b) sLive++;

        for(int xc=RADIUS; xc<NX-RADIUS; xc+=RADIUS){
            for(int yc=RADIUS; yc<NY-RADIUS; yc+=RADIUS){

                int N = 0;
                for(size_t o=0;o<offX.size();o++) if(!m[pix(xc+offX[o],yc+offY[o])]) N++;

                int k1=0, k2=0;
                for(auto &h : kv.second){
                    double dx=h.x-xc, dy=h.y-yc;
                    if(dx*dx+dy*dy > RADIUS*RADIUS) continue;
                    if(h.ne==1) k1++;
                    else if(h.ne==2) k2++;
                }

                sNWin++;
                int f = windowFlags(k1,k2);
                if(!f) continue;

                uint8_t fb = 0;
                for(int i=0;i<NFLAG;i++) if(f & FLAG_BIT[i]){ fb |= (1<<i); sFlagWin[i]++; }
                for(size_t o=0;o<offX.size();o++) fp[pix(xc+offX[o],yc+offY[o])] |= fb;

                wRun=sRun; wLTA=sLTA; wOHDU=sOHDU; wX=xc; wY=yc; wN=N; wK1=k1; wK2=k2; wFlags=f;
                fw->Fill();
            }
        }

        for(int p=0;p<NX*NY;p++){
            if(!fp[p]) continue;
            for(int i=0;i<NFLAG;i++) if(fp[p]>>i & 1){ sFlagPix[i]++; if(!m[p]) sFlagLive[i]++; }
        }
        fs->Fill();
    }

    // ---- copy every tree; calPixTree with the flag bits ORed into mask
    out->cd();
    TIter next(in->GetListOfKeys());
    TKey *key;
    set<string> done;
    while((key=(TKey*)next())){
        string name = key->GetName();
        if(done.count(name)) continue;           // only the highest cycle of each object
        done.insert(name);

        TObject *obj = in->Get(name.c_str());
        if(name=="calPixTree" || name=="hitSumm") continue;   // rewritten below
        if(obj->InheritsFrom(TTree::Class())){
            TTree *t = (TTree*)obj;
            t->SetBranchStatus("*",1);
            out->cd();
            TTree *c = t->CloneTree(-1,"fast");
            c->Write();
            delete c;
        }
        else{
            out->cd();
            obj->Write(name.c_str());
        }
    }

    cal->SetBranchStatus("*",1);
    cal->SetBranchAddress("x",&cx);       cal->SetBranchAddress("y",&cy);
    cal->SetBranchAddress("mask",&cm);    cal->SetBranchAddress("LTANAME",&cl);
    cal->SetBranchAddress("ohdu",&co);    cal->SetBranchAddress("RUNID",&cr);
    out->cd();
    TTree *ncal = cal->CloneTree(0);
    for(Long64_t i=0;i<nCal;i++){
        cal->GetEntry(i);
        if(cx>=0 && cx<NX && cy>=0 && cy<NY){
            auto it = flagPix.find(Img(cr,cl,co));
            if(it!=flagPix.end()){
                uint8_t fb = it->second[pix(cx,cy)];
                for(int b=0;b<NFLAG;b++) if(fb>>b & 1) cm |= FLAG_BIT[b];
            }
        }
        ncal->Fill();
    }
    ncal->Write();

    // hitSumm with flagA/flagB/flagC: the event's flag plus each flag's bit
    // if any of its pixels (xPix, yPix) carries that flag
    int hflag, hn;
    hit->SetBranchStatus("*",1);                  // before sizing: nSavedPix was disabled above
    int maxPix = std::max(1,(int)hit->GetMaximum("nSavedPix"));
    vector<int> xp(maxPix), yp(maxPix);
    hit->SetBranchAddress("runID",&hr);    hit->SetBranchAddress("LTANAME",&hl);
    hit->SetBranchAddress("ohdu",&ho);     hit->SetBranchAddress("flag",&hflag);
    hit->SetBranchAddress("nSavedPix",&hn);
    hit->SetBranchAddress("xPix",xp.data()); hit->SetBranchAddress("yPix",yp.data());
    out->cd();
    TTree *nhit = hit->CloneTree(0);
    int fOut[NFLAG];
    for(int b=0;b<NFLAG;b++)
        nhit->Branch(Form("flag%s",FLAG_NAME[b]),&fOut[b],Form("flag%s/I",FLAG_NAME[b]));
    for(Long64_t i=0;i<hit->GetEntries();i++){
        hit->GetEntry(i);
        uint8_t bits = 0;
        auto it = flagPix.find(Img(hr,hl,ho));
        if(it!=flagPix.end())
            for(int k=0;k<hn;k++)
                if(xp[k]>=0 && xp[k]<NX && yp[k]>=0 && yp[k]<NY) bits |= it->second[pix(xp[k],yp[k])];
        for(int b=0;b<NFLAG;b++) fOut[b] = hflag | ((bits>>b & 1) ? FLAG_BIT[b] : 0);
        nhit->Fill();
    }
    nhit->Write();

    fw->Write();
    fs->Write();
    TNamed info("flagMaskInfo",
        "calPixTree mask bits added by BuildFlagMasks.C on r = 30 windows (standard mask 0x56FC): "
        "0x100000 flag A: k2>=2 or (k2==1 and k1>=2); 0x200000 flag B: k2>=2 or (k2==1 and k1>=3); "
        "0x400000 flag C: k2>=2 or (k2==1 and k1>=5); every pixel of a flagged window's circle gets the bit. "
        "hitSumm: flagA/flagB/flagC = flag plus that flag's bit if any pixel of the event carries it");
    info.Write();

    Long64_t nWinTot = fw->GetEntries();
    out->Close();
    delete out;
    in->Close();
    delete in;
    gSystem->Rename(tmpPath.c_str(),outPath.c_str());

    cout << "wrote " << outPath << "  (" << hits.size() << " images, " << nWinTot << " flagged windows)" << endl;
    return 0;
}

int BuildFlagMasks(const char *inDir="data/March", const char *outDir="data/March_flagged",
                   int part=0, int nParts=1){

    gROOT->SetBatch(kTRUE);
    gSystem->mkdir(outDir,kTRUE);

    vector<string> files;
    TSystemDirectory d(inDir,inDir);
    TList *l = d.GetListOfFiles();
    if(l){
        for(auto o : *l){
            string n = o->GetName();
            if(n.size()>5 && n.substr(n.size()-5)==".root") files.push_back(n);
        }
        delete l;
    }
    sort(files.begin(),files.end());

    int nDone = 0, nFail = 0;
    for(size_t i=0;i<files.size();i++){
        if((int)(i % nParts)!=part) continue;
        int rc = buildFlagMask(string(inDir)+"/"+files[i], string(outDir)+"/"+files[i]);
        if(rc) nFail++; else nDone++;
    }

    cout << "part " << part << "/" << nParts << ": " << nDone << " ok, " << nFail << " failed" << endl;
    return nFail;
}
