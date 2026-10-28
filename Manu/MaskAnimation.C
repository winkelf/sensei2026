// Pictorial walk-through of how the circular-window scan (MaskScanManye.C)
// is applied to one detector image: raw hits -> bad-pixel mask (hits on
// masked pixels vetoed) -> grid of r-pixel windows -> per-window live-pixel
// count N and hit counts k1..k4 -> the resulting per-window maps.
//
// Counts are recomputed here with exactly the rules of MaskScanManye.C and
// checked against the stored windows tree for that image, so every number
// in the frames is the real one.
//
// Output: animation/<category>_e<exposure>_r<radius>_run<run>_LTA<lta>_OHDU<ohdu>/
//   frame_NN.png (one per step) -- frames2gif.py turns them into a GIF.

#include "PlotK1Independence.C"
#include "TEllipse.h"
#include "TH2C.h"
#include "TH2F.h"
#include "TMarker.h"
#include "TBox.h"
#include "TArrow.h"
#include <yaml-cpp/yaml.h>
#include <fstream>

namespace anim {

const int NXP = 3200, NYP = 520;
const int W = 1600, H = 900;

struct Hit{ float x, y; int ne; bool vetoed; };
struct Win{ int xc, yc, N; int k[4]; };

int cSurf, cInk, cInk2, cMuted, cGrid, cBlue, cOrange, cAqua;
int cLive, cMaskIn, cMaskOut, cCircle, cCircleFaint;

// discrete palette shared by every pixel map in a frame:
// 1 = live pixel inside the window, 2 = masked inside, 3 = masked (elsewhere)
void discretePalette(){
    int cols[3] = { cLive, cMaskIn, cMaskOut };
    gStyle->SetPalette(3,cols);
}

void setupColors(){
    applyPlotStyle();
    cSurf   = gStSurf;
    cInk    = gStInk;
    cInk2   = gStInk2;
    cMuted  = gStMuted;
    cGrid   = gStGrid;
    cBlue   = gStBlue;
    cOrange = gStOrange;
    cAqua   = TColor::GetColor("#1baf7a");
    cLive   = TColor::GetColor("#cde2fb");
    cMaskIn = TColor::GetColor("#a8a79f");
    cMaskOut= TColor::GetColor("#dcdbd4");
    cCircle = cBlue;
    cCircleFaint = TColor::GetColorTransparent(cBlue,0.22);
}

int hitColor(int ne){ return ne==1 ? cInk2 : ne==2 ? cBlue : ne==3 ? cOrange : cAqua; }
double hitSize(int ne, double base){ return ne==1 ? base : base*1.9; }

void text(double x, double y, const string &s, int size, int color, int align=11){
    TLatex t;
    t.SetNDC();
    t.SetTextFont(43);
    t.SetTextSize(size);
    t.SetTextColor(color);
    t.SetTextAlign(align);
    t.DrawLatex(x,y,s.c_str());
}

// NDC box of a pad given in canvas pixels (x from left, y from bottom)
TPad* pixelPad(const char *name, double x0, double y0, double x1, double y1){
    TPad *p = new TPad(name,"",x0/W,y0/H,x1/W,y1/H);
    p->SetFillColor(cSurf);
    p->SetBorderMode(0);
    p->SetFrameFillColor(cSurf);
    return p;
}

struct Scene{
    vector<bool> mask;          // NXP*NYP, true = masked
    vector<Hit>  hits;
    vector<Win>  wins;
    int radius;
    int nFull;                  // pixels in a full circle
    string label;               // e.g. "good quads  #upoint  exposure 216000  #upoint  r = 30"
    string imageName;           // e.g. "run 522  #upoint  LTA 14  #upoint  OHDU 3"
    TH2C *hMaskTop = nullptr;
};

bool masked(const Scene &s, int x, int y){ return s.mask[x*NYP+y]; }

// top pad: the whole image, isotropic scale so circles stay round
TPad* drawTop(Scene &s, bool showMask, bool showHits, int grid, const Win *hl){

    // frame width 1504 px -> frame height 1504*520/3200 = 244.4 px
    const double lm=0.045, rm=0.015, tm=0.03, bm=0.17;
    double frameW = W*(1-lm-rm);
    double frameH = frameW*NYP/NXP;
    double padH = frameH/(1-tm-bm);
    double y1 = H-62, y0 = y1-padH;

    TPad *p = pixelPad("top",0,y0,W,y1);
    p->SetLeftMargin(lm); p->SetRightMargin(rm);
    p->SetTopMargin(tm);  p->SetBottomMargin(bm);
    p->Draw();
    p->cd();

    TH1F *fr = p->DrawFrame(0,0,NXP,NYP);
    fr->GetXaxis()->SetTitle("x, y [px]");
    fr->GetYaxis()->SetTitle("");
    fr->GetXaxis()->SetTitleOffset(1.05);
    fr->GetYaxis()->SetTitleOffset(0.55);
    fr->GetYaxis()->CenterTitle();
    fr->GetYaxis()->SetNdivisions(503);
    fr->GetXaxis()->SetTickLength(0.02);
    fr->GetYaxis()->SetTickLength(0.005);

    if(showMask){
        s.hMaskTop->Draw("COL SAME");
    }

    if(grid>0){
        int c = grid==1 ? cCircle : cCircleFaint;
        for(auto &w : s.wins){
            TEllipse *e = new TEllipse(w.xc,w.yc,s.radius,s.radius);
            e->SetFillStyle(0);
            e->SetLineColor(c);
            e->SetLineWidth(1);
            e->Draw();
        }
    }

    if(showHits){
        for(int pass=0;pass<2;pass++)          // 1e first, multi-e on top
        for(auto &h : s.hits){
            if((pass==0)!=(h.ne==1)) continue;
            if(showMask && h.vetoed){
                TMarker *m = new TMarker(h.x,h.y,24);
                m->SetMarkerColor(cMuted);
                m->SetMarkerSize(hitSize(h.ne,0.45));
                m->Draw();
                continue;
            }
            TMarker *m = new TMarker(h.x,h.y,20);
            m->SetMarkerColor(hitColor(h.ne));
            m->SetMarkerSize(hitSize(h.ne,0.45));
            m->Draw();
        }
    }

    if(hl){
        TEllipse *e = new TEllipse(hl->xc,hl->yc,s.radius,s.radius);
        e->SetFillStyle(0);
        e->SetLineColor(cInk);
        e->SetLineWidth(4);
        e->Draw();
    }

    p->RedrawAxis();
    return p;
}

// pixel-level zoom around (cx,cy), half-size z; no margins so every
// primitive is clipped at the pad edge
void drawZoom(Scene &s, int cx, int cy, int z, bool showMask, bool showHits,
               bool showGrid, const Win *hl, double px0, double py0, double side){

    TPad *p = pixelPad("zoom",px0,py0,px0+side,py0+side);
    p->SetMargin(0,0,0,0);
    p->Draw();
    p->cd();

    int x0 = cx-z, y0 = cy-z;
    p->Range(x0-0.5, y0-0.5, x0+2*z+0.5, y0+2*z+0.5);

    // pixel states
    TH2C *h = new TH2C(Form("hz_%d_%d_%d",cx,cy,(int)(hl!=nullptr)),"",2*z+1,x0-0.5,x0+2*z+0.5,2*z+1,y0-0.5,y0+2*z+0.5);
    h->SetDirectory(nullptr);
    for(int x=x0;x<=x0+2*z;x++){
        for(int y=y0;y<=y0+2*z;y++){
            if(x<0||x>=NXP||y<0||y>=NYP){ h->Fill(x,y,3); continue; }
            bool inside = hl && (x-hl->xc)*(x-hl->xc)+(y-hl->yc)*(y-hl->yc) <= s.radius*s.radius;
            bool m = showMask && masked(s,x,y);
            if(inside) h->Fill(x,y, m ? 2 : 1);
            else if(m) h->Fill(x,y,3);
        }
    }
    h->SetMinimum(0.5);
    h->SetMaximum(3.5);
    h->SetContour(3);
    h->Draw("COL A");

    // faint pixel lattice every 10 px, for scale
    for(int x=((x0+9)/10)*10; x<=x0+2*z; x+=10){
        TLine *l = new TLine(x-0.5,y0-0.5,x-0.5,y0+2*z+0.5);
        l->SetLineColor(cGrid); l->Draw();
    }
    for(int y=((y0+9)/10)*10; y<=y0+2*z; y+=10){
        TLine *l = new TLine(x0-0.5,y-0.5,x0+2*z+0.5,y-0.5);
        l->SetLineColor(cGrid); l->Draw();
    }

    if(showGrid){
        for(auto &w : s.wins){
            if(std::abs(w.xc-cx)>z+s.radius || std::abs(w.yc-cy)>z+s.radius) continue;
            TEllipse *e = new TEllipse(w.xc,w.yc,s.radius,s.radius);
            e->SetFillStyle(0);
            e->SetLineColor(hl ? cCircleFaint : cCircle);
            e->SetLineWidth(2);
            e->Draw();
            if(!hl){
                TMarker *c = new TMarker(w.xc,w.yc,2);
                c->SetMarkerColor(cCircle);
                c->SetMarkerSize(1.2);
                c->Draw();
            }
        }
    }

    if(hl){
        TEllipse *e = new TEllipse(hl->xc,hl->yc,s.radius,s.radius);
        e->SetFillStyle(0);
        e->SetLineColor(cInk);
        e->SetLineWidth(4);
        e->Draw();
        TMarker *c = new TMarker(hl->xc,hl->yc,2);
        c->SetMarkerColor(cInk);
        c->SetMarkerSize(1.6);
        c->Draw();
    }

    if(showHits){
        for(auto &hh : s.hits){
            if(std::abs(hh.x-cx)>z+1 || std::abs(hh.y-cy)>z+1) continue;
            bool in = hl && (hh.x-hl->xc)*(hh.x-hl->xc)+(hh.y-hl->yc)*(hh.y-hl->yc) <= s.radius*s.radius;
            if(showMask && hh.vetoed){
                TMarker *m = new TMarker(hh.x,hh.y,24);
                m->SetMarkerColor(cInk2);
                m->SetMarkerSize(hitSize(hh.ne,1.6));
                m->Draw();
                continue;
            }
            TMarker *m = new TMarker(hh.x,hh.y,20);
            m->SetMarkerColor(hitColor(hh.ne));
            m->SetMarkerSize(hitSize(hh.ne, hl && !in ? 1.2 : 1.6));
            m->Draw();
            if(in){
                TMarker *r = new TMarker(hh.x,hh.y,24);
                r->SetMarkerColor(cInk);
                r->SetMarkerSize(hitSize(hh.ne,1.6)+0.9);
                r->Draw();
            }
        }
    }

    // scale bar: radius length
    double sbx = x0+4, sby = y0+4;
    TLine *sb = new TLine(sbx,sby,sbx+s.radius,sby);
    sb->SetLineColor(cInk); sb->SetLineWidth(3); sb->Draw();
    TLatex t;
    t.SetTextFont(43); t.SetTextSize(15); t.SetTextColor(cInk); t.SetTextAlign(11);
    t.DrawLatex(sbx,sby+1.5,Form("%d px",s.radius));

    // frame
    TBox *b = new TBox(x0-0.5,y0-0.5,x0+2*z+0.5,y0+2*z+0.5);
    b->SetFillStyle(0); b->SetLineColor(gStAxis); b->SetLineWidth(1); b->Draw();
}

// marker legend at (x,y) NDC of the current pad, one row
void hitLegend(double x, double y, bool withVeto){
    struct E{ int ne; const char *lab; };
    E es[4] = {{1,"1e"},{2,"2e"},{3,"3e"},{4,"4e+"}};
    double dx = 0.085;
    for(int i=0;i<4;i++){
        TMarker *m = new TMarker(x+i*dx,y,20);
        m->SetNDC();
        m->SetMarkerColor(hitColor(es[i].ne));
        m->SetMarkerSize(hitSize(es[i].ne,1.3));
        m->Draw();
        text(x+i*dx+0.018,y-0.012,es[i].lab,17,cInk2);
    }
    if(withVeto){
        TMarker *m = new TMarker(x+4*dx+0.01,y,24);
        m->SetNDC();
        m->SetMarkerColor(cInk2);
        m->SetMarkerSize(1.3);
        m->Draw();
        text(x+4*dx+0.028,y-0.012,"vetoed (on masked pixel)",17,cInk2);
    }
}

void swatch(double x, double y, int color, const string &lab){
    TBox *b = new TBox(x,y-0.018,x+0.03,y+0.018);
    b->SetFillColor(color);
    b->SetLineColor(gStAxis);
    b->Draw();
    TBox *o = new TBox(x,y-0.018,x+0.03,y+0.018);
    o->SetFillStyle(0);
    o->SetLineColor(gStAxis);
    o->Draw();
    text(x+0.042,y-0.012,lab,17,cInk2);
}

TCanvas* newFrame(const string &step, const string &caption, const Scene &s){
    TCanvas *c = new TCanvas("cf","",W,H);
    c->SetBatch(kTRUE);
    c->SetFillColor(cSurf);
    c->cd();
    text(0.02,0.962,step,26,cInk);
    text(0.98,0.962,s.imageName+"   #upoint   "+s.label,17,cInk2,31);
    text(0.02,0.933,caption,17,cInk2);
    return c;
}

void save(TCanvas *c, const string &dir, int &idx){
    discretePalette();
    string f = Form("%s/frame_%02d.png",dir.c_str(),idx++);
    c->SaveAs(f.c_str());
    delete c;
}

} // namespace anim

using namespace anim;

int maskAnimation(const char *category="goodquads", int exposure=216000, int fileIndex=11,
                   int run=522, int lta=14, int ohdu=3, int radius=30, int maskBits=0x56FC){

    gROOT->SetBatch(kTRUE);
    setupColors();

    YAML::Node cfg = YAML::LoadFile("config.yaml");
    vector<string> files = cfg["exposures"][to_string(exposure)].as<vector<string>>();
    string fileName = files.at(fileIndex-1);

    Scene s;
    s.radius = radius;
    s.label  = string(category)+"  #upoint  exposure "+to_string(exposure)+"  #upoint  r = "+to_string(radius);
    s.imageName = Form("run %d  #upoint  LTA %d  #upoint  OHDU %d",run,lta,ohdu);
    s.mask.assign(NXP*NYP,false);

    TFile *f = TFile::Open(fileName.c_str());
    TTree *cal = (TTree*)f->Get("calPixTree");
    TTree *hit = (TTree*)f->Get("hitSumm");

    int cx,cy,cm,cl,co,cr;
    cal->SetBranchStatus("*",0);
    for(const char *b : {"x","y","mask","LTANAME","ohdu","RUNID"}) cal->SetBranchStatus(b,1);
    cal->SetBranchAddress("x",&cx); cal->SetBranchAddress("y",&cy); cal->SetBranchAddress("mask",&cm);
    cal->SetBranchAddress("LTANAME",&cl); cal->SetBranchAddress("ohdu",&co); cal->SetBranchAddress("RUNID",&cr);
    for(Long64_t i=0;i<cal->GetEntries();i++){
        cal->GetEntry(i);
        if(cr!=run||cl!=lta||co!=ohdu) continue;
        if(cx<0||cx>=NXP||cy<0||cy>=NYP) continue;
        if(cm & maskBits) s.mask[cx*NYP+cy] = true;
    }

    int hr,hl,ho; float he,hx,hy;
    hit->SetBranchStatus("*",0);
    for(const char *b : {"runID","LTANAME","ohdu","e","xBary","yBary"}) hit->SetBranchStatus(b,1);
    hit->SetBranchAddress("runID",&hr); hit->SetBranchAddress("LTANAME",&hl); hit->SetBranchAddress("ohdu",&ho);
    hit->SetBranchAddress("e",&he); hit->SetBranchAddress("xBary",&hx); hit->SetBranchAddress("yBary",&hy);
    for(Long64_t i=0;i<hit->GetEntries();i++){
        hit->GetEntry(i);
        if(hr!=run||hl!=lta||ho!=ohdu) continue;
        int ne = (int)he;
        if(ne<1) continue;
        int ix=(int)round(hx), iy=(int)round(hy);
        if(ix<0||ix>=NXP||iy<0||iy>=NYP) continue;
        s.hits.push_back({hx,hy,std::min(ne,4),s.mask[ix*NYP+iy]});
    }
    f->Close();

    // windows, exactly as MaskScanManye.C
    s.nFull = 0;
    for(int dx=-radius;dx<=radius;dx++) for(int dy=-radius;dy<=radius;dy++) if(dx*dx+dy*dy<=radius*radius) s.nFull++;

    for(int xc=radius; xc<NXP-radius; xc+=radius){
        for(int yc=radius; yc<NYP-radius; yc+=radius){
            Win w{xc,yc,0,{0,0,0,0}};
            for(int dx=-radius;dx<=radius;dx++) for(int dy=-radius;dy<=radius;dy++){
                if(dx*dx+dy*dy>radius*radius) continue;
                if(!s.mask[(xc+dx)*NYP+(yc+dy)]) w.N++;
            }
            for(auto &h : s.hits){
                if(h.vetoed) continue;
                double ddx=h.x-xc, ddy=h.y-yc;
                if(ddx*ddx+ddy*ddy > radius*radius) continue;
                w.k[h.ne-1]++;
            }
            s.wins.push_back(w);
        }
    }

    // cross-check against the stored windows tree for this image
    {
        string wf = Form("%s/e%d/r%d/windowStatisticsManye_%s_%d.root",category,exposure,radius,category,fileIndex);
        TFile F(wf.c_str());
        TTree *t = F.IsZombie() ? nullptr : (TTree*)F.Get("windows");
        if(t){
            int R,L,O,X,Y,N,k1,k2,k3,k4;
            t->SetBranchAddress("runID",&R); t->SetBranchAddress("LTA",&L); t->SetBranchAddress("OHDU",&O);
            t->SetBranchAddress("x",&X); t->SetBranchAddress("y",&Y); t->SetBranchAddress("N",&N);
            t->SetBranchAddress("k1",&k1); t->SetBranchAddress("k2",&k2); t->SetBranchAddress("k3",&k3); t->SetBranchAddress("k4",&k4);
            map<pair<int,int>,const Win*> byPos;
            for(auto &w : s.wins) byPos[{w.xc,w.yc}] = &w;
            int nMatch=0, nDiff=0;
            for(Long64_t i=0;i<t->GetEntries();i++){
                t->GetEntry(i);
                if(R!=run||L!=lta||O!=ohdu) continue;
                auto it = byPos.find({X,Y});
                if(it==byPos.end()){ nDiff++; continue; }
                const Win *w = it->second;
                if(w->N==N && w->k[0]==k1 && w->k[1]==k2 && w->k[2]==k3 && w->k[3]==k4) nMatch++; else nDiff++;
            }
            cout << "Cross-check vs " << wf << ": " << nMatch << " windows identical, "
                 << nDiff << " different (of " << s.wins.size() << " recomputed)" << endl;
        }
        else cout << "Cross-check skipped: cannot read " << wf << endl;
    }

    s.hMaskTop = new TH2C("hMaskTop","",NXP,-0.5,NXP-0.5,NYP,-0.5,NYP-0.5);
    s.hMaskTop->SetDirectory(nullptr);
    for(int x=0;x<NXP;x++) for(int y=0;y<NYP;y++) if(s.mask[x*NYP+y]) s.hMaskTop->SetBinContent(x+1,y+1,3);
    s.hMaskTop->SetMinimum(0.5);
    s.hMaskTop->SetMaximum(3.5);
    s.hMaskTop->SetContour(3);

    int nHits=s.hits.size(), nVeto=0, nByE[4]={0,0,0,0};
    for(auto &h : s.hits){ if(h.vetoed) nVeto++; nByE[h.ne-1]++; }
    long nMasked=0; for(bool b : s.mask) if(b) nMasked++;
    double maskedPct = 100.0*nMasked/(NXP*NYP);

    // windows to walk through: the most centred window on each multi-e
    // hit, a few partly masked ones, a few clean ones, one fully masked
    vector<const Win*> pick;
    auto add = [&](const Win *w){
        for(auto *p : pick) if(p==w) return;
        for(auto *p : pick) if(std::abs(p->xc-w->xc)<2*radius && std::abs(p->yc-w->yc)<2*radius) return;
        pick.push_back(w);
    };
    for(auto &h : s.hits){
        if(h.ne<2 || h.vetoed) continue;
        const Win *best=nullptr; double bd=1e9;
        for(auto &w : s.wins){ double d=hypot(h.x-w.xc,h.y-w.yc); if(d<bd){bd=d;best=&w;} }
        if(best && bd<=radius) add(best);
        if(pick.size()>=5) break;
    }
    int nPart=0, nClean=0, nZero=0;
    for(auto &w : s.wins){
        double fr = (double)w.N/s.nFull;
        int nk = w.k[0]+w.k[1]+w.k[2]+w.k[3];
        if(nPart<3 && fr>0.25 && fr<0.65 && nk>=1){ size_t b=pick.size(); add(&w); if(pick.size()>b) nPart++; }
        else if(nClean<2 && fr>0.97 && w.k[0]>=2 && w.k[1]+w.k[2]+w.k[3]==0){ size_t b=pick.size(); add(&w); if(pick.size()>b) nClean++; }
        else if(nZero<1 && w.N==0){ size_t b=pick.size(); add(&w); if(pick.size()>b) nZero++; }
    }
    std::sort(pick.begin(),pick.end(),[](const Win*a,const Win*b){ return a->xc<b->xc || (a->xc==b->xc && a->yc<b->yc); });

    // zoom target for the intro frames: first window with a 2e hit
    const Win *focus = pick.empty() ? &s.wins[s.wins.size()/2] : pick[0];
    for(auto *w : pick) if(w->k[1]>0){ focus=w; break; }

    string dir = Form("animation/%s_e%d_r%d_run%d_LTA%d_OHDU%d",category,exposure,radius,run,lta,ohdu);
    gSystem->mkdir(dir.c_str(),kTRUE);
    gSystem->Exec(("rm -f "+dir+"/frame_*.png").c_str());

    int idx = 1;
    const double zx=40, zy=14, zs=505;       // zoom pad (canvas px)
    const double ix0=(zx+zs+50)/W;           // info column start (NDC)
    const int zIntro = 70, zWin = 45;

    // 1) raw image
    {
        TCanvas *c = newFrame("1.  The image: every reconstructed hit",
            Form("%d hits on the %d #times %d px CCD, coloured by charge (e = number of electrons)",nHits,NXP,NYP),s);
        drawTop(s,false,true,0,nullptr);
        c->cd();
        drawZoom(s,focus->xc,focus->yc,zIntro,false,true,false,nullptr,zx,zy,zs);
        c->cd();
        text(ix0,0.60,"Zoom (left) around one of the 2e hits.",19,cInk);
        text(ix0,0.55,Form("1e: %d    2e: %d    3e: %d    4e+: %d",nByE[0],nByE[1],nByE[2],nByE[3]),19,cInk2);
        hitLegend(ix0+0.01,0.47,false);
        save(c,dir,idx);
    }

    // 2) mask
    {
        TCanvas *c = newFrame(Form("2.  Apply the bad-pixel mask (bits 0x%X)",maskBits),
            Form("%.1f%% of the pixels are masked; a hit whose barycentre lands on a masked pixel is dropped (%d of %d here)",maskedPct,nVeto,nHits),s);
        drawTop(s,true,true,0,nullptr);
        c->cd();
        drawZoom(s,focus->xc,focus->yc,zIntro,true,true,false,nullptr,zx,zy,zs);
        c->cd();
        swatch(ix0,0.60,cMaskOut,"masked pixel");
        text(ix0,0.52,"Only live (unmasked) pixels and the hits on them are used.",19,cInk);
        hitLegend(ix0+0.01,0.45,true);
        save(c,dir,idx);
    }

    // 3) window grid
    {
        TCanvas *c = newFrame(Form("3.  Lay the window grid: circles of r = %d px, centres every %d px",radius,radius),
            Form("%zu windows on this image; neighbouring circles overlap, so each pixel belongs to several windows",s.wins.size()),s);
        drawTop(s,true,true,1,nullptr);
        c->cd();
        drawZoom(s,focus->xc,focus->yc,zIntro,true,true,true,nullptr,zx,zy,zs);
        c->cd();
        text(ix0,0.60,Form("Each circle holds %d pixels when nothing is masked.",s.nFull),19,cInk);
        text(ix0,0.55,Form("Centres (#times) run from %d to %d in x and %d to %d in y, step %d.",
                            radius,s.wins.back().xc,radius,s.wins.back().yc,radius),19,cInk2);
        hitLegend(ix0+0.01,0.47,true);
        save(c,dir,idx);
    }

    // 4) walk through selected windows
    int step=0;
    for(auto *w : pick){
        step++;
        TCanvas *c = newFrame(Form("4.  Count, window by window  (%d/%zu)",step,pick.size()),
            "N = live pixels inside the circle;   k1, k2, k3, k4+ = hits with 1, 2, 3, #geq4 electrons whose barycentre is inside the circle and on a live pixel",s);
        drawTop(s,true,true,2,w);
        c->cd();
        drawZoom(s,w->xc,w->yc,zWin,true,true,true,w,zx,zy,zs);
        c->cd();

        double fr = (double)w->N/s.nFull;
        text(ix0,0.60,Form("Window centred at (x, y) = (%d, %d)",w->xc,w->yc),22,cInk);
        swatch(ix0,0.535,cLive,"");
        text(ix0+0.042,0.523,Form("N = %d live pixels   (%.0f%% of %d)",w->N,100*fr,s.nFull),20,cInk);
        swatch(ix0,0.48,cMaskIn,"");
        text(ix0+0.042,0.468,Form("%d masked pixels inside the circle",s.nFull-w->N),20,cInk2);

        const char *kn[4] = {"k_{1}","k_{2}","k_{3}","k_{4+}"};
        for(int i=0;i<4;i++){
            double yy = 0.39-0.055*i;
            TMarker *m = new TMarker(ix0+0.012,yy+0.006,20);
            m->SetNDC();
            m->SetMarkerColor(hitColor(i+1));
            m->SetMarkerSize(hitSize(i+1,1.3));
            m->Draw();
            string q = w->N==0 ? "q undefined (N = 0)" : w->k[i]==0 ? "q = 0" : Form("q = k/N = %.2e",(double)w->k[i]/w->N);
            text(ix0+0.042,yy-0.006,Form("%s = %d",kn[i],w->k[i]),21,w->k[i]>0 ? cInk : cMuted);
            text(ix0+0.13,yy-0.006,q,18,cInk2);
        }
        text(ix0,0.13,"Ringed hits are the ones counted; open circles were vetoed by the mask.",16,cMuted);
        text(ix0,0.09, w->N>0 ? "This (N, k_{1}..k_{4+}) is one row of the windows tree."
                              : "N = 0: fully masked, so this window is not written to the windows tree.",16,cMuted);
        save(c,dir,idx);
    }

    // 5) result maps (continuous palette: separate canvas, own palette)
    {
        int nx=0, ny=0;
        for(int xc=radius; xc<NXP-radius; xc+=radius) nx++;
        for(int yc=radius; yc<NYP-radius; yc+=radius) ny++;
        double xlo = radius-radius/2.0, ylo = radius-radius/2.0;
        TH2F *hFr = new TH2F("hFr","",nx,xlo,xlo+nx*radius,ny,ylo,ylo+ny*radius);
        TH2F *hK  = new TH2F("hK","",nx,xlo,xlo+nx*radius,ny,ylo,ylo+ny*radius);
        hFr->SetDirectory(nullptr); hK->SetDirectory(nullptr);
        double kMax=0;
        for(auto &w : s.wins){
            hFr->Fill(w.xc,w.yc,std::max(1e-6,(double)w.N/s.nFull));
            if(w.N>0) hK->Fill(w.xc,w.yc,std::max(1e-6,(double)w.k[0]));
            kMax = std::max(kMax,(double)w.k[0]);
        }

        TCanvas *c = new TCanvas("cf","",W,H);
        c->SetBatch(kTRUE);
        c->SetFillColor(cSurf);
        text(0.02,0.962,"5.  Result: one (N, k_{1}..k_{4+}) per window",26,cInk);
        text(0.98,0.962,s.imageName+"   #upoint   "+s.label,17,cInk2,31);
        text(0.02,0.933,"Each cell is one window (centre at the cell centre); these per-window values feed the distributions, cut scans and independence test",17,cInk2);

        const double lm=0.07, rm=0.09, tm=0.10, bm=0.16;
        double frameW = W*(1-lm-rm);
        double frameH = frameW*NYP/NXP;
        double padH = frameH/(1-tm-bm);

        double yTop1 = H-70;
        TPad *p1 = pixelPad("m1",0,yTop1-padH,W,yTop1);
        TPad *p2 = pixelPad("m2",0,yTop1-2*padH-10,W,yTop1-padH-10);
        for(TPad *p : {p1,p2}){ p->SetLeftMargin(lm); p->SetRightMargin(rm); p->SetTopMargin(tm); p->SetBottomMargin(bm); p->Draw(); }

        // sequential blue ramp (light = low)
        const int nSt = 5;
        double st[nSt] = {0.0,0.25,0.5,0.75,1.0};
        double rr[nSt], gg[nSt], bb[nSt];
        const char *hx[nSt] = {"#e8f1fd","#9ec5f4","#5598e7","#256abf","#0d366b"};
        for(int i=0;i<nSt;i++){ TColor *col=gROOT->GetColor(TColor::GetColor(hx[i])); rr[i]=col->GetRed(); gg[i]=col->GetGreen(); bb[i]=col->GetBlue(); }
        TColor::CreateGradientColorTable(nSt,st,rr,gg,bb,99);
        gStyle->SetNumberContours(99);

        auto style2 = [&](TH2F *h, const char *zt, double zmax){
            h->SetMinimum(0);
            h->SetMaximum(zmax);
            h->GetXaxis()->SetTitle("x, y [px]");
            h->GetYaxis()->SetTitle("");
            h->GetXaxis()->SetTitleOffset(1.05);
            h->GetYaxis()->SetTitleOffset(1.25);
            h->GetYaxis()->CenterTitle();
            h->GetYaxis()->SetNdivisions(503);
            h->GetZaxis()->SetTitle(zt);
            h->GetZaxis()->SetTitleOffset(0.9);
            h->GetZaxis()->SetNdivisions(503);
            h->GetXaxis()->SetTickLength(0.02);
            h->GetYaxis()->SetTickLength(0.005);
        };

        p1->cd();
        style2(hFr,"N / N_{full}",1.0);
        hFr->Draw("COLZ");
        text(lm,1-tm+0.03,Form("Live fraction of each window,  N / %d",s.nFull),18,cInk);

        p2->cd();
        style2(hK,"k_{1}",std::max(1.0,kMax));
        hK->Draw("COLZ");
        for(auto &w : s.wins){
            if(w.k[1]+w.k[2]+w.k[3]==0) continue;
            TMarker *m = new TMarker(w.xc,w.yc,w.k[1]>0 ? 20 : 21);
            m->SetMarkerColor(w.k[2]+w.k[3]>0 ? cOrange : cInk);
            m->SetMarkerSize(0.7);
            m->Draw();
        }
        text(lm,1-tm+0.03,"k_{1} per window   (dots: windows with k_{2} #geq 1 in ink, with k_{3} or k_{4+} in orange)",18,cInk);

        string fn = Form("%s/frame_%02d.png",dir.c_str(),idx++);
        c->SaveAs(fn.c_str());
        delete c;
    }

    cout << "Wrote " << idx-1 << " frames to " << dir << "/" << endl;

    return 0;
}

int MaskAnimation(){ return maskAnimation(); }
