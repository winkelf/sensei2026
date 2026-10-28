// Charge map of one OHDU (quadrant) with the standard mask and the
// window-flag masks written by BuildFlagMasks.C.
//
// Every pixel is coloured by its charge ePix (calPixTree): white at 0 e,
// yellow at 1 e, then orange to deep red at 5 e or more. The standard mask (0x56FC) is drawn in translucent grey
// and the window-flag masks in translucent blue over the charge, so masked
// hits stay visible underneath; the
// flags are nested (C inside B inside A): every flag masks the windows with
// k2 >= 2, and A / B / C add k2 = 1 windows with k1 >= 2 / 3 / 5. Each pixel
// shows the strictest flag among the windows covering it: light = A only
// (k2 = 1, k1 = 2), medium = A and B (k2 = 1, k1 = 3-4), dark = A, B and C
// (k2 >= 2, or k2 = 1 and k1 >= 5).
// The 3200 x 520 quadrant is split into two rows (x 0-1600 and 1600-3200) so
// a CCD pixel is about one screen pixel.
//
// With lec = true the blue layer shows the LEC mask instead (bit 0x2000,
// already set in the reconstruction; the mask of the LEC reference in the
// q-cut scans), for comparison with the flags.
//
// Written as PNG: the maps are raster data (as a vector PDF every pixel is a
// separate shape, tens of MB and very heavy to render).
//
// usage (from the project folder), reading a flagged copy:
//   root -l -b -q 'FlagDisplay.C("data/March_flagged/<file>.root",522,14,3,"flagdisplays/x.png","Old")'

#include "PlotStyle.C"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH2C.h"
#include "TExec.h"
#include "TBox.h"
#include "TH1F.h"

namespace fdisp {
const int NX = 3200, NY = 520;
const int STD_MASK = 0x56FC;
const int FLAG_BIT[3] = { 0x100000, 0x200000, 0x400000 };
const int LEC_BIT = 0x2000;
const double QMAX = 5.0;          // charge at which the scale is fully red
int maskCols[4];
}

// palettes switched inside the pads with TExec (charge, then masks)
void fdChargePalette(){
    // stops at 0, 1, 2, 3, 4, 5 e: white, yellow, then orange to deep red
    const int n = 6;
    double st[n] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    const char *hx[n] = {"#ffffff","#fde74c","#fdb03a","#f37b2c","#e0452a","#b8141c"};
    double r[n], g[n], b[n];
    for(int i=0;i<n;i++){ TColor *c = gROOT->GetColor(TColor::GetColor(hx[i])); r[i]=c->GetRed(); g[i]=c->GetGreen(); b[i]=c->GetBlue(); }
    TColor::CreateGradientColorTable(n,st,r,g,b,50);
    gStyle->SetNumberContours(50);
}

void fdMaskPalette(){
    gStyle->SetPalette(4,fdisp::maskCols);
    gStyle->SetNumberContours(4);
}

int flagDisplay(const char *file, int run, int lta, int ohdu, const char *outFile, const char *datasetLabel="",
                bool lec=false){

    using namespace fdisp;
    gROOT->SetBatch(kTRUE);
    applyPlotStyle();

    TFile f(file);
    if(f.IsZombie()){ cout << "cannot open " << file << endl; return 1; }
    TTree *cal = (TTree*)f.Get("calPixTree");
    TTree *fw  = (TTree*)f.Get("flagWindows");
    if(!cal || !fw){ cout << "not a flagged copy (needs calPixTree and flagWindows): " << file << endl; return 1; }

    // charge and masks of this OHDU
    vector<float> charge(NX*NY,0);
    vector<int>   mask(NX*NY,0);
    int cx, cy, cm, cl, co, cr;
    double ce;
    cal->SetBranchStatus("*",0);
    for(const char *b : {"x","y","ePix","mask","LTANAME","ohdu","RUNID"}) cal->SetBranchStatus(b,1);
    cal->SetBranchAddress("x",&cx); cal->SetBranchAddress("y",&cy); cal->SetBranchAddress("ePix",&ce);
    cal->SetBranchAddress("mask",&cm); cal->SetBranchAddress("LTANAME",&cl);
    cal->SetBranchAddress("ohdu",&co); cal->SetBranchAddress("RUNID",&cr);
    for(Long64_t i=0;i<cal->GetEntries();i++){
        cal->GetEntry(i);
        if(cr!=run || cl!=lta || co!=ohdu) continue;
        if(cx<0||cx>=NX||cy<0||cy>=NY) continue;
        charge[cx*NY+cy] = ce;
        mask[cx*NY+cy]   = cm;
    }

    int nWinFlag[3] = {0,0,0};
    int wr, wl, wo, wf;
    fw->SetBranchAddress("runID",&wr); fw->SetBranchAddress("LTA",&wl);
    fw->SetBranchAddress("OHDU",&wo);  fw->SetBranchAddress("flags",&wf);
    for(Long64_t i=0;i<fw->GetEntries();i++){
        fw->GetEntry(i);
        if(wr!=run || wl!=lta || wo!=ohdu) continue;
        for(int k=0;k<3;k++) if(wf & FLAG_BIT[k]) nWinFlag[k]++;
    }

    // mask layer: 1 = standard mask, 2 = flag A only, 3 = A and B, 4 = A, B and C
    // (LEC mode: 1 = standard mask, 4 = LEC)
    long nPix = (long)NX*NY, nStd = 0, nFlag[3] = {0,0,0}, nLec = 0;
    vector<char> layer(NX*NY,0);
    for(int p=0;p<NX*NY;p++){
        int m = mask[p];
        if(m & STD_MASK){ layer[p] = 1; nStd++; continue; }
        if(m & LEC_BIT) nLec++;
        for(int k=0;k<3;k++) if(m & FLAG_BIT[k]) nFlag[k]++;
        if(lec){ if(m & LEC_BIT) layer[p] = 4; continue; }
        if(m & FLAG_BIT[2]) layer[p] = 4;
        else if(m & FLAG_BIT[1]) layer[p] = 3;
        else if(m & FLAG_BIT[0]) layer[p] = 2;
    }

    maskCols[0] = TColor::GetColorTransparent(TColor::GetColor("#8f8e87"),0.55);
    maskCols[1] = TColor::GetColorTransparent(gStBlue,0.25);
    maskCols[2] = TColor::GetColorTransparent(gStBlue,0.45);
    maskCols[3] = TColor::GetColorTransparent(gStBlue,0.70);

    // canvas: two rows, ~1 screen pixel per CCD pixel
    const int CW = 1850;
    const double lm = 0.045, rm = 0.09, tm = 0.05, bm = 0.10;
    const double frameW = CW*(1-lm-rm);
    const int    half   = NX/2;
    const double frameH = frameW*NY/half;
    const double padH   = frameH/(1-tm-bm);
    const int headH = 145;
    const int CH = headH + (int)std::ceil(2*padH) + 10;

    TCanvas *c = new TCanvas("cflag","",CW,CH);
    c->SetBatch(kTRUE);

    string ctx = string(datasetLabel)+(strlen(datasetLabel) ? "  #upoint  " : "")+Form("run %d  #upoint  LTA %d  #upoint  OHDU %d",run,lta,ohdu);
    if(lec)
        canvasHeader(c,CH,"Pixel charge with the standard mask and the LEC mask",ctx,
            Form("live pixels: %.1f%% with the standard mask  #upoint  LEC removes %.1f%% more (%.1f%% live)  #upoint  for comparison, flag A removes %.1f%%, B %.1f%%, C %.1f%%",
                 100.0*(nPix-nStd)/nPix,100.0*nLec/nPix,100.0*(nPix-nStd-nLec)/nPix,100.0*nFlag[0]/nPix,100.0*nFlag[1]/nPix,100.0*nFlag[2]/nPix));
    else
        canvasHeader(c,CH,"Pixel charge with the standard mask and the window-flag masks",ctx,
            Form("live pixels: %.1f%% with the standard mask  #upoint  flag A removes %.1f%% (%d windows)  #upoint  flag B %.1f%% (%d)  #upoint  flag C %.1f%% (%d)",
                 100.0*(nPix-nStd)/nPix,100.0*nFlag[0]/nPix,nWinFlag[0],100.0*nFlag[1]/nPix,nWinFlag[1],100.0*nFlag[2]/nPix,nWinFlag[2]));

    // mask legend
    {
        c->cd();
        struct L{ int col; const char *txt; };
        vector<L> items;
        if(lec) items = {
            {maskCols[0],"standard mask 0x56FC"},
            {maskCols[3],"LEC mask (bit 0x2000), where not already masked"} };
        else items = {
            {maskCols[0],"standard mask 0x56FC"},
            {maskCols[1],"flag A only:  windows with k_{2} = 1, k_{1} = 2"},
            {maskCols[2],"flags A and B:  k_{2} = 1, k_{1} = 3 or 4"},
            {maskCols[3],"flags A, B and C:  k_{2} #geq 2, or k_{2} = 1, k_{1} #geq 5"} };
        double x = 0.015, y = 1.0-(headH-44.0)/CH;
        double bw = 22.0/CW, bh = 14.0/CH;
        TLatex t;
        t.SetNDC();
        t.SetTextFont(43);
        t.SetTextSize(16);
        t.SetTextColor(gStInk2);
        t.SetTextAlign(12);
        for(auto &it : items){
            TBox *b = new TBox(x,y-bh/2,x+bw,y+bh/2);
            b->SetFillColor(it.col);
            b->SetLineColor(gStAxis);
            b->Draw();
            TBox *o = new TBox(x,y-bh/2,x+bw,y+bh/2);
            o->SetFillStyle(0); o->SetLineColor(gStAxis); o->Draw();
            t.DrawLatex(x+bw+0.006,y,it.txt);
            x += 0.235;
        }
        t.SetTextColor(gStMuted);
        t.SetTextSize(15);
        t.DrawLatex(0.015,1.0-(headH-18.0)/CH, lec ?
            "LEC = Local Event Cut: mask bit 0x2000 set by the reconstruction, the mask behind the LEC reference lines of the q-cut scans." :
            "Every flag masks the windows with k_{2} #geq 2; A, B and C also mask k_{2} = 1 windows with k_{1} #geq 2, 3 and 5. "
            "A pixel shows the strictest flag among the r = 30 windows covering it.");
    }

    for(int row=0;row<2;row++){

        int x0 = row*half, x1 = x0+half;
        double y1 = CH-headH-row*padH, y0 = y1-padH;

        TPad *p = new TPad(Form("prow%d",row),"",0,y0/CH,1,y1/CH);
        p->SetFillColor(TColor::GetColor("#ffffff"));
        p->SetFrameFillColor(TColor::GetColor("#ffffff"));
        p->SetLeftMargin(lm); p->SetRightMargin(rm); p->SetTopMargin(tm); p->SetBottomMargin(bm);
        c->cd();
        p->Draw();
        p->cd();

        TH2F *hq = new TH2F(Form("hq%d",row),"",half,x0-0.5,x1-0.5,NY,-0.5,NY-0.5);
        TH2C *hm = new TH2C(Form("hm%d",row),"",half,x0-0.5,x1-0.5,NY,-0.5,NY-0.5);
        hq->SetDirectory(nullptr);
        hm->SetDirectory(nullptr);
        for(int x=x0;x<x1;x++) for(int y=0;y<NY;y++){
            int pidx = x*NY+y;
            double q = charge[pidx];
            if(q>0) hq->SetBinContent(x-x0+1,y+1,std::min(q,QMAX*0.9999));
            if(layer[pidx]) hm->SetBinContent(x-x0+1,y+1,layer[pidx]);
        }
        hq->SetMinimum(0); hq->SetMaximum(QMAX);
        hm->SetMinimum(0.5); hm->SetMaximum(4.5);

        TH1F *fr = p->DrawFrame(x0-0.5,-0.5,x1-0.5,NY-0.5);
        fr->GetXaxis()->SetTitle(row==1 ? "x [px]" : "");
        fr->GetYaxis()->SetTitle("y [px]");
        fr->GetYaxis()->CenterTitle();
        fr->GetYaxis()->SetNdivisions(505);
        fr->GetXaxis()->SetNdivisions(516);
        hq->GetZaxis()->SetTitle("charge [e]");
        hq->GetZaxis()->SetNdivisions(505);

        TExec *e1 = new TExec(Form("exq%d",row),"fdChargePalette();");
        e1->Draw();
        hq->Draw("COLZ SAME");
        TExec *e2 = new TExec(Form("exm%d",row),"fdMaskPalette();");
        e2->Draw();
        hm->SetContour(4);
        hm->Draw("COL SAME");
        p->RedrawAxis();

        padHeader(Form("x = %d to %d",x0,x1),"",0.012);
    }

    // PNG at twice the canvas resolution: crisp text, each CCD pixel an
    // exact 2 x 2 block
    gStyle->SetImageScaling(2.);
    c->SaveAs(outFile);
    gStyle->SetImageScaling(1.);
    delete c;
    string o = outFile;
    if(o.size()>4 && o.substr(o.size()-4)==".pdf") tightPdf(o);
    cout << "wrote " << outFile << endl;
    return 0;
}

int FlagDisplay(){ return 0; }
