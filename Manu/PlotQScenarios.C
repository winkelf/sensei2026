// Two concrete selection scenarios instead of arbitrary weighted combos:
//   "no cut"  -- every window, unconditionally
//   "k2 < 2"  -- only windows with k2 == 0 or k2 == 1
//
// For each (category, exposure, radius) this produces:
//   - qdist_nocut_<tag>.pdf   -- q1,q2,q3,q4 distributions, no cut
//   - qdist_k2lt2_<tag>.pdf   -- q1,q2,q3,q4 distributions, k2<2 only
//   - qdist_compare_<tag>.pdf -- both scenarios overlaid on the same 4 panels
//   - qcut_scan_nocut_<tag>.pdf  -- independent per-channel q_i cut scan, no cut
//   - qcut_scan_k2lt2_<tag>.pdf  -- same scan, restricted to the k2<2 subset
//
// No LEC references, no channel weighting -- this replaces the older
// weighted-combo machinery (which lived in PlotQ1CutScan.C) now that the
// analysis direction has moved to these two concrete selections.
//
// Restricted to a single radius at a time -- r30 first, other radii can be
// added the same way. #include's PlotWindowStats.C to reuse
// buildChain()/discoverTaggedSubdirs()/styleAndDraw().

#include <numeric>

#include "PlotWindowStats.C"
#include "PlotStyle.C"
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveText.h"

const int kNScanSteps = 300;
const double kMaxQCut = 0.002;
const int   kChanColor[4] = { kBlue+1, kRed+1, kGreen+2, kMagenta+1 };
const char *kChanName[4]  = { "k_{1}", "k_{2}", "k_{3}", "k_{4+}" };
const char *kQName[4]     = { "q_{1}", "q_{2}", "q_{3}", "q_{4+}" };
const char *kQFormula[4]  = { "q_{1} = k_{1}/N", "q_{2} = k_{2}/N", "q_{3} = k_{3}/N", "q_{4+} = k_{4+}/N" };

// reference values from the LEC ("Local Event Cut") mask -- not affected by
// any q_i cut, so it's just one number per channel, drawn as a dashed
// horizontal line to help pick a working point. Requires
// MaskScanManyeLEC.C's output at <category>/lec/e<exposure>/r<radius>/.
struct LecTotals{
    bool     valid = false;
    Long64_t sumN  = 0;
    Long64_t sumK[4] = {0,0,0,0};

    double q(int i) const { return sumN>0 ? (double)sumK[i]/sumN : 0; }
};

LecTotals loadLecTotals(const string &label, int exposure, int radius){

    LecTotals t;

    string dir = label+"/lec/e"+to_string(exposure)+"/r"+to_string(radius);

    TChain *ch = buildChain(dir,label);

    Long64_t nFiles = ch->GetListOfFiles() ? ch->GetListOfFiles()->GetEntries() : 0;

    if(nFiles==0){
        cout << "  [LEC] no files found in '" << dir << "/' -- reference lines will be omitted"
             << " (run MaskScanManyeLEC.C first)" << endl;
        delete ch;
        return t;
    }

    int N,k1,k2,k3,k4;
    ch->SetBranchAddress("N", &N);
    ch->SetBranchAddress("k1",&k1);
    ch->SetBranchAddress("k2",&k2);
    ch->SetBranchAddress("k3",&k3);
    ch->SetBranchAddress("k4",&k4);

    Long64_t n = ch->GetEntries();

    for(Long64_t i=0;i<n;i++){
        ch->GetEntry(i);
        t.sumN    += N;
        t.sumK[0] += k1;
        t.sumK[1] += k2;
        t.sumK[2] += k3;
        t.sumK[3] += k4;
    }

    t.valid = true;

    delete ch;

    return t;
}

// draws a dashed horizontal reference line spanning [xmin,xmax] at y=value,
// in the given color, and adds it to the legend
TLine* referenceLine(double xmin, double xmax, double value, int color, TLegend *leg, const string &label){

    TLine *l = new TLine(xmin,value,xmax,value);
    l->SetLineColor(color);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();

    if(leg) leg->AddEntry(l,label.c_str(),"l");

    return l;
}

// all per-window quantities for one (category, exposure, radius) combo,
// loaded once and reused for both the distribution and the scan plots
struct WindowData{
    Long64_t n = 0;
    vector<int>    N;
    vector<int>    K[4];
    vector<double> Q[4]; // = K[c]/N, already stored as r1-r4 in the tree
};

WindowData loadWindowData(const string &label, const string &dir){

    WindowData wd;

    TChain *ch = buildChain(dir,label);

    Long64_t nFiles = ch->GetListOfFiles() ? ch->GetListOfFiles()->GetEntries() : 0;

    if(nFiles==0){
        delete ch;
        return wd;
    }

    int N,k1,k2,k3,k4;
    double r1,r2,r3,r4;

    ch->SetBranchAddress("N", &N);
    ch->SetBranchAddress("k1",&k1);
    ch->SetBranchAddress("k2",&k2);
    ch->SetBranchAddress("k3",&k3);
    ch->SetBranchAddress("k4",&k4);
    ch->SetBranchAddress("r1",&r1);
    ch->SetBranchAddress("r2",&r2);
    ch->SetBranchAddress("r3",&r3);
    ch->SetBranchAddress("r4",&r4);

    Long64_t n = ch->GetEntries();
    wd.n = n;
    wd.N.reserve(n);
    for(int c=0;c<4;c++){ wd.K[c].reserve(n); wd.Q[c].reserve(n); }

    for(Long64_t i=0;i<n;i++){

        ch->GetEntry(i);

        wd.N.push_back(N);

        int    kv[4] = {k1,k2,k3,k4};
        double rv[4] = {r1,r2,r3,r4};

        for(int c=0;c<4;c++){
            wd.K[c].push_back(kv[c]);
            wd.Q[c].push_back(rv[c]);
        }
    }

    delete ch;

    return wd;
}

// smallest non-zero and largest q of channel c over idx
void qRange(const WindowData &wd, int c, const vector<Long64_t> &idx, double &lo, double &hi){
    lo = 0; hi = 0;
    for(Long64_t i : idx){
        double q = wd.Q[c][i];
        if(q<=0) continue;
        if(lo==0 || q<lo) lo=q;
        if(q>hi) hi=q;
    }
    if(hi<=0){ lo=1e-6; hi=1e-5; }
    if(lo>=hi) lo=hi/10;
}

Long64_t qZeros(const WindowData &wd, int c, const vector<Long64_t> &idx){
    Long64_t n=0;
    for(Long64_t i : idx) if(wd.Q[c][i]<=0) n++;
    return n;
}

// log-spaced bins over [lo,hi]: q spans several decades (a few nearly
// fully masked windows reach q = 1), so a linear axis squashes everything
// into the first bin. q = 0 windows can't go on a log axis and are counted
// separately (qZeros).
TH1D* buildQHist(const string &name, const WindowData &wd, int c,
                  const vector<Long64_t> &idx, double lo, double hi){

    const int nb = 40;
    vector<double> edges(nb+1);
    double l0 = std::log10(lo*0.95), l1 = std::log10(hi*1.05);
    for(int i=0;i<=nb;i++) edges[i] = std::pow(10.0, l0+(l1-l0)*i/nb);

    TH1D *h = new TH1D(name.c_str(),"",nb,edges.data());
    h->SetDirectory(nullptr);
    for(Long64_t i : idx) if(wd.Q[c][i]>0) h->Fill(wd.Q[c][i]);
    return h;
}

void plotQDistributionScenarios(const string &label, const string &context, const string &dir,
                                  const string &tag, const WindowData &wd,
                                  const vector<Long64_t> &idxNoCut, const vector<Long64_t> &idxK2lt2){

    applyPlotStyle();

    const char *rateName[4] = { "1e rate", "2e rate", "3e rate", "#geq4e rate" };
    const int CW = 1600, CH = 470;

    auto drawSingle = [&](const vector<Long64_t> &idx, const string &scenarioTag, const string &scenarioTitle){
        size_t nWin = idx.size();
        TCanvas *c = new TCanvas(("cq_"+scenarioTag+"_"+tag).c_str(),"",CW,CH);
        c->SetBatch(kTRUE);
        double used = canvasHeader(c,CH,"q = k / N per window,  "+scenarioTitle,context,
            Form("%zu of %lld windows",nWin,(Long64_t)wd.n));
        vector<TPad*> pads = gridPads(4,1,1-used-0.01,"pq_"+scenarioTag+"_"+tag+"_");
        for(int i=0;i<4;i++){
            double lo, hi;
            qRange(wd,i,idx,lo,hi);
            TH1D *h = buildQHist(Form("hQ%d_%s_%s",i+1,scenarioTag.c_str(),tag.c_str()),wd,i,idx,lo,hi);
            pads[i]->cd();
            gPad->SetLogx();
            drawDistribution(h,chanColor(i),chanFill(i),kQFormula[i]);
            padHeader(rateName[i],Form("+ %lld at q = 0",qZeros(wd,i,idx)));
        }
        string out = dir+"/qdist_"+scenarioTag+"_"+tag+".pdf";
        c->SaveAs(out.c_str());
        delete c;
        tightPdf(out);
        cout << "  wrote " << out << endl;
    };

    drawSingle(idxNoCut, "nocut", "no cut");
    drawSingle(idxK2lt2, "k2lt2", "k_{2} < 2");

    // overlay: channel-coloured fill = no cut, ink outline = k2<2, on a
    // shared per-channel x-range so the shapes are directly comparable
    TCanvas *cc = new TCanvas(("cqcompare_"+tag).c_str(),"",CW,CH);
    cc->SetBatch(kTRUE);
    double used = canvasHeader(cc,CH,"q = k / N per window:  no cut vs k_{2} < 2",context,
        Form("no cut: %zu windows     k_{2} < 2: %zu windows (%.2f%%)",
             idxNoCut.size(),idxK2lt2.size(),100.0*idxK2lt2.size()/wd.n));
    vector<TPad*> pads = gridPads(4,1,1-used-0.01,"pqc_"+tag+"_");

    for(int c=0;c<4;c++){

        double loA, hiA, loB, hiB;
        qRange(wd,c,idxNoCut,loA,hiA);
        qRange(wd,c,idxK2lt2,loB,hiB);
        double lo = std::min(loA,loB), hi = std::max(hiA,hiB);

        TH1D *hA = buildQHist(Form("hQ%d_cmpA_%s",c+1,tag.c_str()),wd,c,idxNoCut,lo,hi);
        TH1D *hB = buildQHist(Form("hQ%d_cmpB_%s",c+1,tag.c_str()),wd,c,idxK2lt2,lo,hi);

        pads[c]->cd();
        gPad->SetLogx();
        drawDistribution(hA,chanColor(c),chanFill(c),kQFormula[c]);
        hB->SetLineColor(gStInk);
        hB->SetLineWidth(2);
        hB->SetLineStyle(1);
        hB->SetFillStyle(0);
        hB->Draw("HIST SAME");
        gPad->RedrawAxis();
        padHeader(rateName[c],"");

        TLegend *leg = styledLegend(0.40,0.72,0.95,0.87);
        leg->AddEntry(hA,Form("no cut  (+ %lld at q = 0)",qZeros(wd,c,idxNoCut)),"f");
        leg->AddEntry(hB,Form("k_{2} < 2  (+ %lld at q = 0)",qZeros(wd,c,idxK2lt2)),"l");
        leg->Draw();
    }

    string outCmp = dir+"/qdist_compare_"+tag+".pdf";
    cc->SaveAs(outCmp.c_str());
    delete cc;
    tightPdf(outCmp);

    cout << "  wrote " << outCmp << endl;
}

// per-channel independent q_i cut scan (no weighting, no LEC): cut q_c<cut,
// sum k_c and N over the surviving windows -- computed only within idx
struct ScanResult{
    vector<double> cutVals[4], sumK[4], sumN[4];
    double scanMax[4];
};

ScanResult computeScan(const WindowData &wd, const vector<Long64_t> &idx){

    ScanResult res;

    double maxQ[4] = {0,0,0,0};
    for(Long64_t i : idx)
        for(int c=0;c<4;c++)
            if(wd.Q[c][i]>maxQ[c]) maxQ[c]=wd.Q[c][i];

    for(int c=0;c<4;c++) res.scanMax[c] = std::min(maxQ[c], kMaxQCut);

    for(int c=0;c<4;c++){

        res.cutVals[c].assign(kNScanSteps+1,0);
        res.sumK[c].assign(kNScanSteps+1,0);
        res.sumN[c].assign(kNScanSteps+1,0);

        for(int s=0;s<=kNScanSteps;s++){

            double cut = res.scanMax[c] * s / kNScanSteps;
            res.cutVals[c][s] = cut;

            Long64_t ksum=0, Nsum=0;

            for(Long64_t i : idx){
                if(wd.Q[c][i] >= cut) continue;
                ksum += wd.K[c][i];
                Nsum += wd.N[i];
            }

            res.sumK[c][s] = (double)ksum;
            res.sumN[c][s] = (double)Nsum;
        }
    }

    return res;
}

void plotScanPanel(const ScanResult &res, const LecTotals &lec, const string &context, const string &scenarioLabel,
                     const string &out, Long64_t nSel, Long64_t nTotal){

    applyPlotStyle();

    double globalScanMax = *std::max_element(res.scanMax,res.scanMax+4);

    const int CW = 1600, CH = 520;
    TCanvas *c = new TCanvas(("cscan_"+out).c_str(),"",CW,CH);
    c->SetBatch(kTRUE);

    double used = canvasHeader(c,CH,"Independent q_{i} cut scan,  "+scenarioLabel,context,
        Form("%lld of %lld windows in this selection (%.2f%%)     each channel keeps the windows with q_{i} below the cut     dashed = LEC mask reference",
             nSel,nTotal,100.0*nSel/nTotal));

    vector<TPad*> pads = gridPads(4,1,1-used-0.01,"ps_"+out+"_");

    auto curve = [&](int i, const vector<double> &x, const vector<double> &y){
        TGraph *g = new TGraph((int)x.size(),x.data(),y.data());
        g->SetLineColor(chanColor(i));
        g->SetLineWidth(2);
        return g;
    };
    auto lecLine = [&](double xmax, double v, int color){
        if(v<=0) return;
        TLine *l = new TLine(0,v,xmax,v);
        l->SetLineColor(color);
        l->SetLineStyle(2);
        l->SetLineWidth(2);
        l->Draw();
    };
    auto frame = [&](double ylo, double yhi, const char *ytitle, bool logy=true){
        if(logy) gPad->SetLogy();
        TH1F *fr = gPad->DrawFrame(0,ylo,globalScanMax,yhi);
        fr->GetXaxis()->SetTitle("q_{i} cut");
        fr->GetYaxis()->SetTitle(ytitle);
        fr->GetXaxis()->SetNdivisions(505);
    };

    // 1: k_i kept
    pads[0]->cd();
    {
        double gmax=1;
        for(int i=0;i<4;i++) for(double v : res.sumK[i]) gmax=std::max(gmax,v);
        frame(0.5,gmax*4,"k_{i} kept");
        if(lec.valid) for(int i=0;i<4;i++) lecLine(res.scanMax[i],(double)lec.sumK[i],chanColor(i));
        for(int i=0;i<4;i++) curve(i,res.cutVals[i],res.sumK[i])->Draw("L");
        gPad->RedrawAxis();
        padHeader("Hits kept","");
    }

    // 2: N kept -- one curve per channel; the LEC N is common to all
    pads[1]->cd();
    {
        double gmax=1;
        for(int i=0;i<4;i++) for(double v : res.sumN[i]) gmax=std::max(gmax,v);
        if(lec.valid) gmax=std::max(gmax,(double)lec.sumN);
        frame(0,gmax*1.12,"N kept (live pixels)",false);
        gPad->SetGridy();
        if(lec.valid) lecLine(globalScanMax,(double)lec.sumN,gStInk2);
        for(int i=0;i<4;i++) curve(i,res.cutVals[i],res.sumN[i])->Draw("L");
        gPad->RedrawAxis();
        padHeader("Live pixels kept","");
    }

    // 3: q_i = k_i/N kept
    pads[2]->cd();
    {
        double gmin=1, gmax=1e-12;
        vector<double> rc[4], rq[4];
        for(int i=0;i<4;i++){
            for(int s=0;s<=kNScanSteps;s++){
                if(res.sumN[i][s]<=0 || res.sumK[i][s]<=0) continue;
                rc[i].push_back(res.cutVals[i][s]);
                rq[i].push_back(res.sumK[i][s]/res.sumN[i][s]);
                gmax=std::max(gmax,rq[i].back());
                gmin=std::min(gmin,rq[i].back());
            }
            if(lec.valid && lec.q(i)>0) gmin=std::min(gmin,lec.q(i));
        }
        frame(gmin*0.3,gmax*4,"q_{i} = k_{i} / N kept");
        if(lec.valid) for(int i=0;i<4;i++) lecLine(res.scanMax[i],lec.q(i),chanColor(i));
        for(int i=0;i<4;i++) curve(i,rc[i],rq[i])->Draw("L");
        gPad->RedrawAxis();
        padHeader("Rate kept, q_{i}","");
    }

    // 4: LEC reference values and scan ranges
    pads[3]->cd();
    {
        TLatex t;
        t.SetNDC();
        t.SetTextFont(43);
        t.SetTextAlign(11);
        t.SetTextSize(18);
        t.SetTextColor(gStInk);
        t.DrawLatex(0.06,0.84,"Scan (solid) and LEC mask reference (dashed)");
        t.SetTextSize(16);
        t.SetTextColor(gStInk2);
        if(!lec.valid){
            t.DrawLatex(0.06,0.72,"not available -- run MaskScanManyeLEC.C");
        }
        else{
            t.DrawLatex(0.06,0.73,Form("LEC  N = %.3g live pixels",(double)lec.sumN));
            t.DrawLatex(0.06,0.62,"channel");
            t.DrawLatex(0.36,0.62,"LEC k");
            t.DrawLatex(0.55,0.62,"LEC q");
            t.DrawLatex(0.78,0.62,"scan to");
            const char *nm[4] = {"k_{1}","k_{2}","k_{3}","k_{4+}"};
            for(int i=0;i<4;i++){
                double y = 0.51-0.11*i;
                TLine *sw = new TLine(0.06,y+0.012,0.11,y+0.012);
                sw->SetNDC(); sw->SetLineColor(chanColor(i)); sw->SetLineWidth(3); sw->Draw();
                TLine *sd = new TLine(0.13,y+0.012,0.18,y+0.012);
                sd->SetNDC(); sd->SetLineColor(chanColor(i)); sd->SetLineWidth(3); sd->SetLineStyle(2); sd->Draw();
                t.SetTextColor(gStInk);
                t.DrawLatex(0.21,y,nm[i]);
                t.DrawLatex(0.36,y,Form("%lld",(Long64_t)lec.sumK[i]));
                t.DrawLatex(0.55,y,Form("%.2e",lec.q(i)));
                t.SetTextColor(gStInk2);
                t.DrawLatex(0.78,y,Form("%.4f",res.scanMax[i]));
            }
        }
    }

    c->SaveAs(out.c_str());
    delete c;
    tightPdf(out);

    cout << "  wrote " << out << endl;
}

void plotQScenarios(const string &label, const string &title, int exposure, int radius){

    string dir = label+"/e"+to_string(exposure)+"/r"+to_string(radius);
    string tag = label+"_e"+to_string(exposure)+"_r"+to_string(radius);
    string fullTitle = title+"  #upoint  exposure "+to_string(exposure)+"  #upoint  r = "+to_string(radius);

    cout << "[" << tag << "] q scenarios (no cut / k2<2) from '" << dir << "/' ..." << endl;

    WindowData wd = loadWindowData(label,dir);

    if(wd.n==0){
        cout << "  no windowStatisticsManye_" << label << "_*.root files found -- skipping" << endl;
        return;
    }

    vector<Long64_t> idxNoCut(wd.n);
    std::iota(idxNoCut.begin(),idxNoCut.end(),0);

    vector<Long64_t> idxK2lt2;
    for(Long64_t i=0;i<wd.n;i++) if(wd.K[1][i]<2) idxK2lt2.push_back(i);

    cout << "  no cut: " << idxNoCut.size() << " windows; k2<2: " << idxK2lt2.size()
         << " windows (" << Form("%.3f",100.0*idxK2lt2.size()/wd.n) << "%)" << endl;

    if(idxK2lt2.empty()){
        cout << "  no windows with k2<2 -- skipping" << endl;
        return;
    }

    plotQDistributionScenarios(label,fullTitle,dir,tag,wd,idxNoCut,idxK2lt2);

    LecTotals lec = loadLecTotals(label,exposure,radius);

    ScanResult scanNoCut = computeScan(wd,idxNoCut);
    ScanResult scanK2lt2 = computeScan(wd,idxK2lt2);

    plotScanPanel(scanNoCut,lec,fullTitle,"no cut", dir+"/qcut_scan_nocut_"+tag+".pdf", (Long64_t)idxNoCut.size(),wd.n);
    plotScanPanel(scanK2lt2,lec,fullTitle,"k_{2} < 2", dir+"/qcut_scan_k2lt2_"+tag+".pdf", (Long64_t)idxK2lt2.size(),wd.n);
}

void plotQScenariosForRadius(const string &label, const string &title, int radius){

    vector<int> exposures = discoverTaggedSubdirs(label,'e');

    if(exposures.empty()){
        cout << "[" << label << "] no e<N> subfolders found under '" << label << "/' -- skipping" << endl;
        return;
    }

    for(int exposure : exposures){

        string expDir = label+"/e"+to_string(exposure);

        vector<int> radii = discoverTaggedSubdirs(expDir,'r');

        if(std::find(radii.begin(),radii.end(),radius)==radii.end()){
            cout << "[" << label << "/e" << exposure << "] no r" << radius << " subfolder -- skipping" << endl;
            continue;
        }

        plotQScenarios(label,title,exposure,radius);
    }
}

int PlotQScenarios(){

    gROOT->SetBatch(kTRUE);

    plotQScenariosForRadius("goodquads",  "good quads",   30);
    plotQScenariosForRadius("notbadquads","not-bad quads",30);

    return 0;
}
