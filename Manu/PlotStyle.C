// Shared look for every plot: pixel-sized fonts (same size in every pad,
// whatever its size), no default title boxes, hairline grids, and colours
// fixed per quantity so the same thing looks the same everywhere:
//   charge channel  1e = dark grey, 2e = blue, 3e = orange, 4e+ = aqua
//   N (live pixels) = neutral grey
// Include this from any plotting macro and call applyPlotStyle() first.

#ifndef PLOTSTYLE_C
#define PLOTSTYLE_C

#include <string>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TH1.h"
#include "TLine.h"

int gStInk, gStInk2, gStMuted, gStGrid, gStAxis, gStSurf, gStBlue, gStBlueBand, gStOrange, gStAqua;

void applyPlotStyle(){

    gStInk      = TColor::GetColor("#0b0b0b");
    gStInk2     = TColor::GetColor("#52514e");
    gStMuted    = TColor::GetColor("#898781");
    gStGrid     = TColor::GetColor("#e1e0d9");
    gStAxis     = TColor::GetColor("#c3c2b7");
    gStSurf     = TColor::GetColor("#fcfcfb");
    gStBlue     = TColor::GetColor("#2a78d6");
    gStBlueBand = TColor::GetColor("#b7d3f6");
    gStOrange   = TColor::GetColor("#eb6834");
    gStAqua     = TColor::GetColor("#1baf7a");

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    gStyle->SetCanvasColor(gStSurf);
    gStyle->SetPadColor(gStSurf);
    gStyle->SetFrameFillColor(gStSurf);
    gStyle->SetFrameLineColor(gStAxis);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);

    gStyle->SetTextFont(43);
    gStyle->SetLabelFont(43,"XYZ");
    gStyle->SetTitleFont(43,"XYZ");
    gStyle->SetLegendFont(43);
    gStyle->SetLabelSize(17,"XYZ");
    gStyle->SetTitleSize(19,"XYZ");
    gStyle->SetLabelColor(gStInk2,"XYZ");
    gStyle->SetTitleColor(gStInk,"XYZ");
    gStyle->SetAxisColor(gStAxis,"XYZ");
    gStyle->SetTitleOffset(1.25,"X");
    gStyle->SetTitleOffset(1.7,"Y");

    gStyle->SetGridColor(gStGrid);
    gStyle->SetGridStyle(1);
    gStyle->SetGridWidth(1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetPadLeftMargin(0.11);
    gStyle->SetPadRightMargin(0.04);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.12);

    gStyle->SetEndErrorSize(0);
    gStyle->SetLineScalePS(1.6);
    gStyle->SetLegendBorderSize(0);

    gROOT->ForceStyle();
}

// channel i = 0..3 for 1e, 2e, 3e, 4e+
int chanColor(int i){
    switch(i){
        case 0:  return gStInk2;
        case 1:  return gStBlue;
        case 2:  return gStOrange;
        default: return gStAqua;
    }
}

int chanFill(int i){ return TColor::GetColorTransparent(chanColor(i), i==0 ? 0.22 : 0.30); }

int nColor(){ return gStInk2; }
int nFill(){  return TColor::GetColorTransparent(gStMuted,0.25); }

TLegend* styledLegend(double x1, double y1, double x2, double y2){
    TLegend *l = new TLegend(x1,y1,x2,y2);
    l->SetFillStyle(0);
    l->SetBorderSize(0);
    l->SetTextFont(43);
    l->SetTextSize(16);
    l->SetTextColor(gStInk2);
    return l;
}

// title (left, ink) and context (right, secondary ink) just above the frame
void padHeader(const std::string &left, const std::string &right, double dy=0.025){
    TLatex t;
    t.SetNDC();
    t.SetTextFont(43);
    double y = 1-gPad->GetTopMargin()+dy;
    t.SetTextSize(19);
    t.SetTextColor(gStInk);
    t.SetTextAlign(11);
    t.DrawLatex(gPad->GetLeftMargin(),y,left.c_str());
    if(!right.empty()){
        t.SetTextSize(16);
        t.SetTextColor(gStInk2);
        t.SetTextAlign(31);
        t.DrawLatex(1-gPad->GetRightMargin(),y,right.c_str());
    }
}

// canvas-level title strip: title and context on the first line, an
// optional summary line under it. Returns the NDC height it used.
double canvasHeader(TCanvas *c, int heightPx, const std::string &title,
                     const std::string &context, const std::string &summary=""){
    c->cd();
    double lineH = 30.0/heightPx;
    TLatex t;
    t.SetNDC();
    t.SetTextFont(43);
    t.SetTextSize(22);
    t.SetTextColor(gStInk);
    t.SetTextAlign(13);
    t.DrawLatex(0.015,1-8.0/heightPx,title.c_str());
    t.SetTextSize(16);
    t.SetTextColor(gStInk2);
    t.SetTextAlign(33);
    t.DrawLatex(0.985,1-11.0/heightPx,context.c_str());
    double used = lineH+10.0/heightPx;
    if(!summary.empty()){
        t.SetTextAlign(13);
        t.DrawLatex(0.015,1-used,summary.c_str());
        used += 24.0/heightPx;
    }
    return used;
}

// ncol x nrow pads filling [0,yTop] of the current canvas, row-major from
// the top-left
std::vector<TPad*> gridPads(int ncol, int nrow, double yTop, const std::string &prefix){
    std::vector<TPad*> pads;
    double h = yTop/nrow;
    for(int r=0;r<nrow;r++){
        for(int c=0;c<ncol;c++){
            TPad *p = new TPad((prefix+std::to_string(r*ncol+c)).c_str(),"",
                               (double)c/ncol, yTop-(r+1)*h, (double)(c+1)/ncol, yTop-r*h);
            p->SetFillColor(gStSurf);
            p->SetLeftMargin(0.17);
            p->SetRightMargin(0.04);
            p->SetTopMargin(0.11);
            p->SetBottomMargin(0.16);
            p->Draw();
            pads.push_back(p);
        }
    }
    return pads;
}

// a single distribution: filled, outlined histogram, log y, hairline grid
void drawDistribution(TH1 *h, int lineColor, int fillColor, const char *xTitle,
                       const char *yTitle="Windows", bool logy=true){
    h->SetStats(0);
    h->SetLineColor(lineColor);
    h->SetLineWidth(2);
    h->SetFillColor(fillColor);
    h->SetFillStyle(1001);
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle(yTitle);
    h->GetXaxis()->SetNdivisions(505);
    if(logy){
        gPad->SetLogy();
        double mx = h->GetMaximum();
        h->SetMinimum(0.5);
        h->SetMaximum(std::max(2.0,mx*4));
    }
    else gPad->SetGridy();
    h->Draw("HIST");
    gPad->RedrawAxis();
}

// ROOT writes landscape canvases as a rotated page whose box doesn't match
// the drawing, so viewers show large blank margins -- re-render through
// pdftocairo (poppler), which gives a page exactly the size of the canvas.
// Left untouched if pdftocairo isn't available.
void tightPdf(const std::string &path){
    std::string tmp = path+".tmp.pdf";
    if(gSystem->Exec(("pdftocairo -pdf '"+path+"' '"+tmp+"' 2>/dev/null").c_str())==0)
        gSystem->Rename(tmp.c_str(),path.c_str());
    else
        gSystem->Unlink(tmp.c_str());
}

#endif
