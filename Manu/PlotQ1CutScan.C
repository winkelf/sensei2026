// For each (category, exposure, radius), and for each channel i in
// {1,2,3,4} INDEPENDENTLY (no weighting, no merging), scans a q_i = k_i/N
// cut from 0 up to min(kMaxQCut, highest observed q_i), keeping only
// windows with q_i < cut, and plots (all channels overlaid, one PDF per
// combination):
//   1) k_i (summed over kept windows) vs. its own q_i cut, for i=1..4
//   2) N (summed over kept windows) vs. its own q_i cut, for i=1..4 --
//      "living pixels" kept, one curve per channel since each channel
//      admits a different set of windows
//   3) q_i = (1)/(2) vs. its own q_i cut, for i=1..4
// Each curve also gets a single reference value computed once with the
// LEC ("Local Event Cut") mask bit added to the standard mask -- that
// selection isn't affected by the cut, so it's just one number per
// channel, drawn as a dashed horizontal line. kMaxQCut is intentionally
// small (tuned around where the curves cross their LEC reference) so that
// crossing region is clearly visible instead of being squashed by the
// full range of the data.
//
// Requires MaskScanManye.C's standard output (<category>/e<exp>/r<radius>/)
// and, for the reference lines, MaskScanManyeLEC.C's output
// (<category>/lec/e<exp>/r<radius>/).

#include "PlotWindowStats.C"
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveText.h"

const int kNCutSteps = 300;

// upper bound of each channel's scan -- kept small and tuned by hand
// around where the curves cross their LEC reference (currently ~0.001 for
// q1, goodquads/e216000/r30), so that crossing region fills the plot
// instead of being squashed. Each channel still uses its own true max if
// that happens to be smaller than this cap.
const double kMaxQCut = 0.002;

const int   kChannelColor[4] = { kBlue+1, kRed+1, kGreen+2, kMagenta+1 };
const char *kChannelName[4]  = { "k_{1}", "k_{2}", "k_{3}", "k_{4+}" };

// weights for the combined statistic shown (for reference only) in
// plotQDistribution -- q_weighted = (k1 + w2*k2 + w3*(k3+k4)) / N. The
// scan itself is per-channel and unweighted; this is just an extra
// distribution to compare against.
const double kWeight2 = 5.0;
const double kWeight3 = 10.0;
const int    kWeightedColor = kBlack;
const char  *kWeightedName  = "q (weighted: k_{1}+5k_{2}+10(k_{3}+k_{4}))";

struct LecTotals{
    bool     valid = false;
    Long64_t sumN  = 0;
    Long64_t sumK[4] = {0,0,0,0};

    double q(int i) const { return sumN>0 ? (double)sumK[i]/sumN : 0; }
    double weightedQ() const {
        return sumN>0 ? (sumK[0] + kWeight2*sumK[1] + kWeight3*(sumK[2]+sumK[3]))/sumN : 0;
    }
};

// sums k1-k4 and N over every window in the LEC-mask output for this
// (category, exposure, radius) -- no cut applied, computed once
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

void plotQCutScanForCombo(const string &label, const string &title, int exposure, int radius){

    string dir = label+"/e"+to_string(exposure)+"/r"+to_string(radius);
    string tag = label+"_e"+to_string(exposure)+"_r"+to_string(radius);
    string fullTitle = title+", exposure="+to_string(exposure)+", r="+to_string(radius);

    cout << "[" << tag << "] independent per-channel q cut scan from '" << dir << "/' ..." << endl;

    TChain *ch = buildChain(dir,label);

    Long64_t nFiles = ch->GetListOfFiles() ? ch->GetListOfFiles()->GetEntries() : 0;

    if(nFiles==0){
        cout << "  no windowStatisticsManye_" << label << "_*.root files found -- skipping" << endl;
        delete ch;
        return;
    }

    int N,k1,k2,k3,k4;

    ch->SetBranchAddress("N", &N);
    ch->SetBranchAddress("k1",&k1);
    ch->SetBranchAddress("k2",&k2);
    ch->SetBranchAddress("k3",&k3);
    ch->SetBranchAddress("k4",&k4);

    Long64_t n = ch->GetEntries();

    // pull the whole (small) per-window table into memory once
    vector<int> vN(n), vK[4];
    for(int i=0;i<4;i++) vK[i].assign(n,0);

    vector<double> vQ[4];
    for(int i=0;i<4;i++) vQ[i].assign(n,0.0);

    double maxQ[4] = {0,0,0,0};

    for(Long64_t i=0;i<n;i++){
        ch->GetEntry(i);

        vN[i] = N;
        int kv[4] = {k1,k2,k3,k4};

        for(int c=0;c<4;c++){
            vK[c][i] = kv[c];
            double q = (N>0) ? (double)kv[c]/N : 0;
            vQ[c][i] = q;
            if(q>maxQ[c]) maxQ[c]=q;
        }
    }

    delete ch;

    double scanMax[4];
    for(int c=0;c<4;c++) scanMax[c] = std::min(maxQ[c], kMaxQCut);

    cout << "  " << n << " windows; max q per channel = ("
         << maxQ[0] << ", " << maxQ[1] << ", " << maxQ[2] << ", " << maxQ[3]
         << "), scanning up to (" << scanMax[0] << ", " << scanMax[1] << ", "
         << scanMax[2] << ", " << scanMax[3] << ")" << endl;

    // scan each channel independently: cut on q_c < cut, sum k_c and N
    // over the windows that pass
    vector<double> cutVals[4], sumKv[4], sumNv[4];

    for(int c=0;c<4;c++){

        cutVals[c].assign(kNCutSteps+1,0);
        sumKv[c].assign(kNCutSteps+1,0);
        sumNv[c].assign(kNCutSteps+1,0);

        for(int s=0; s<=kNCutSteps; s++){

            double cut = scanMax[c] * s / kNCutSteps;
            cutVals[c][s] = cut;

            Long64_t ksum=0, Nsum=0;

            for(Long64_t i=0;i<n;i++){
                if(vQ[c][i] >= cut) continue;
                ksum += vK[c][i];
                Nsum += vN[i];
            }

            sumKv[c][s] = (double)ksum;
            sumNv[c][s] = (double)Nsum;
        }
    }

    LecTotals lec = loadLecTotals(label,exposure,radius);

    double globalScanMax = *std::max_element(scanMax,scanMax+4);

    //////////////////////////////////////////////////////////
    // One combined canvas: k_i, N_i, q_i (each vs its own cut), stats box
    //////////////////////////////////////////////////////////

    TCanvas *c = new TCanvas(("cqcut_"+tag).c_str(),"",1400,1000);
    c->SetBatch(kTRUE);
    c->Divide(2,2);

    // pad 1: k_i vs its own cut, all 4 channels overlaid
    c->cd(1);
    gPad->SetLogy();
    gPad->SetGrid();
    {
        TGraph *g[4];
        double gmax = 1;

        for(int i=0;i<4;i++){
            g[i] = new TGraph(kNCutSteps+1, cutVals[i].data(), sumKv[i].data());
            g[i]->SetLineColor(kChannelColor[i]);
            g[i]->SetLineWidth(2);
            for(double v : sumKv[i]) if(v>gmax) gmax=v;
        }

        g[0]->SetTitle((fullTitle+";q_{i} cut;k_{i} kept").c_str());
        g[0]->GetYaxis()->SetRangeUser(0.5,gmax*3);
        g[0]->GetXaxis()->SetLimits(0,globalScanMax);
        g[0]->Draw("AL");
        for(int i=1;i<4;i++) g[i]->Draw("L SAME");

        TLegend *leg = new TLegend(0.55,0.15,0.98,0.45);
        for(int i=0;i<4;i++) leg->AddEntry(g[i],kChannelName[i],"l");

        if(lec.valid)
            for(int i=0;i<4;i++)
                referenceLine(0,scanMax[i],(double)lec.sumK[i],kChannelColor[i],leg,
                              string(kChannelName[i])+" (LEC)");

        leg->Draw();
    }

    // pad 2: N_i (living pixels kept) vs its own cut, all 4 channels overlaid
    c->cd(2);
    gPad->SetLogy();
    gPad->SetGrid();
    {
        TGraph *g[4];
        double gmax = 1;

        for(int i=0;i<4;i++){
            g[i] = new TGraph(kNCutSteps+1, cutVals[i].data(), sumNv[i].data());
            g[i]->SetLineColor(kChannelColor[i]);
            g[i]->SetLineWidth(2);
            for(double v : sumNv[i]) if(v>gmax) gmax=v;
        }

        g[0]->SetTitle((fullTitle+";q_{i} cut;N kept").c_str());
        g[0]->GetYaxis()->SetRangeUser(1,gmax*3);
        g[0]->GetXaxis()->SetLimits(0,globalScanMax);
        g[0]->Draw("AL");
        for(int i=1;i<4;i++) g[i]->Draw("L SAME");

        TLegend *leg = new TLegend(0.55,0.15,0.98,0.45);
        for(int i=0;i<4;i++) leg->AddEntry(g[i],(string("N (")+kChannelName[i]+" cut)").c_str(),"l");

        if(lec.valid)
            for(int i=0;i<4;i++)
                referenceLine(0,scanMax[i],(double)lec.sumN,kChannelColor[i],leg,
                              string("N (LEC, ")+kChannelName[i]+")");

        leg->Draw();
    }

    // pad 3: q_i = k_i/N vs its own cut, all 4 channels overlaid
    c->cd(3);
    gPad->SetLogy();
    gPad->SetGrid();
    {
        TGraph *g[4];
        double gmax = 1e-9;
        vector<double> rc[4], rq[4];

        for(int i=0;i<4;i++){
            for(int s=0;s<=kNCutSteps;s++){
                if(sumNv[i][s]<=0) continue;
                rc[i].push_back(cutVals[i][s]);
                rq[i].push_back(sumKv[i][s]/sumNv[i][s]);
            }
            g[i] = new TGraph((int)rc[i].size(), rc[i].data(), rq[i].data());
            g[i]->SetLineColor(kChannelColor[i]);
            g[i]->SetLineWidth(2);
            for(double v : rq[i]) if(v>gmax) gmax=v;
        }

        g[0]->SetTitle((fullTitle+";q_{i} cut;q_{i} = k_{i}/N (kept windows)").c_str());
        g[0]->GetYaxis()->SetRangeUser(gmax*1e-6,gmax*3);
        g[0]->GetXaxis()->SetLimits(0,globalScanMax);
        g[0]->Draw("AL");
        for(int i=1;i<4;i++) g[i]->Draw("L SAME");

        TLegend *leg = new TLegend(0.55,0.65,0.98,0.95);
        for(int i=0;i<4;i++) leg->AddEntry(g[i],("q "+string(kChannelName[i])).c_str(),"l");

        if(lec.valid)
            for(int i=0;i<4;i++)
                referenceLine(0,scanMax[i],lec.q(i),kChannelColor[i],leg,
                              "q "+string(kChannelName[i])+" (LEC)");

        leg->Draw();
    }

    // pad 4: stats box
    c->cd(4);
    {
        TPaveText *pt = new TPaveText(0.02,0.02,0.98,0.98,"NDC");
        pt->SetTextAlign(12);
        pt->SetBorderSize(0);
        pt->SetFillColor(0);
        pt->AddText(Form("Category: %s",fullTitle.c_str()));
        pt->AddText("q_{i} = k_{i} / N, scanned independently per channel (no weighting)");
        pt->AddText(Form("Windows: %lld",(Long64_t)n));
        for(int i=0;i<4;i++)
            pt->AddText(Form("%s: max q = %.5f (scan capped at %.5f)",
                              kChannelName[i],maxQ[i],scanMax[i]));
        if(lec.valid)
            for(int i=0;i<4;i++)
                pt->AddText(Form("LEC %s: k = %lld, N = %lld, q = %.5f",
                                  kChannelName[i],(Long64_t)lec.sumK[i],(Long64_t)lec.sumN,lec.q(i)));
        else
            pt->AddText("LEC: not available (run MaskScanManyeLEC.C)");
        pt->Draw();
    }

    string out = dir+"/qcut_scan_summary_"+tag+".pdf";
    c->SaveAs(out.c_str());
    delete c;

    cout << "  wrote " << out << endl;
}

// safety-check plot: the full distribution of q1, q2, q3, q4 (each = k_i/N,
// independently, no weighting), overlaid (log-x so both the bulk and any
// outlier tail -- e.g. the near-fully-masked N=1 windows seen before -- are
// visible in one view), with each channel's LEC reference q and the scan's
// kMaxQCut boundary marked for context
void plotQDistribution(const string &label, const string &title, int exposure, int radius){

    string dir = label+"/e"+to_string(exposure)+"/r"+to_string(radius);
    string tag = label+"_e"+to_string(exposure)+"_r"+to_string(radius);
    string fullTitle = title+", exposure="+to_string(exposure)+", r="+to_string(radius);

    cout << "[" << tag << "] q distributions from '" << dir << "/' ..." << endl;

    TChain *ch = buildChain(dir,label);

    Long64_t nFiles = ch->GetListOfFiles() ? ch->GetListOfFiles()->GetEntries() : 0;

    if(nFiles==0){
        cout << "  no windowStatisticsManye_" << label << "_*.root files found -- skipping" << endl;
        delete ch;
        return;
    }

    int N,k1,k2,k3,k4;

    ch->SetBranchAddress("N", &N);
    ch->SetBranchAddress("k1",&k1);
    ch->SetBranchAddress("k2",&k2);
    ch->SetBranchAddress("k3",&k3);
    ch->SetBranchAddress("k4",&k4);

    Long64_t n = ch->GetEntries();

    // index 4 = the weighted-combined q, shown here for reference even
    // though the scan itself is per-channel and unweighted
    vector<vector<double>> vals(5, vector<double>(n));

    double maxQ = 0, minPosQ = -1;

    for(Long64_t i=0;i<n;i++){
        ch->GetEntry(i);

        int kv[4] = {k1,k2,k3,k4};

        for(int c=0;c<4;c++){
            double q = (N>0) ? (double)kv[c]/N : 0;
            vals[c][i] = q;
            if(q>maxQ) maxQ=q;
            if(q>0 && (minPosQ<0 || q<minPosQ)) minPosQ=q;
        }

        double qw = (N>0) ? (k1 + kWeight2*k2 + kWeight3*(k3+k4))/N : 0;
        vals[4][i] = qw;
        if(qw>maxQ) maxQ=qw;
        if(qw>0 && (minPosQ<0 || qw<minPosQ)) minPosQ=qw;
    }

    delete ch;

    if(minPosQ<=0 || minPosQ>=maxQ) minPosQ = maxQ>0 ? maxQ*1e-6 : 1e-9;

    cout << "  " << n << " windows, q range (nonzero, across all channels) ["
         << minPosQ << ", " << maxQ << "]" << endl;

    LecTotals lec = loadLecTotals(label,exposure,radius);

    // extend the low edge to include the smallest LEC reference if it
    // falls below the smallest observed window q, so its line is visible
    // only extend for LEC references within ~2 orders of magnitude of the
    // observed data -- a handful of LEC channels can have only 0-2 raw
    // counts over hundreds of millions of pixels, so their q is orders of
    // magnitude below anything meaningful and would otherwise squash the
    // whole plot into a sliver just to include that one far-off line
    double rangeMin = minPosQ;
    if(lec.valid){
        for(int i=0;i<4;i++)
            if(lec.q(i) > minPosQ*1e-2 && lec.q(i) < rangeMin) rangeMin=lec.q(i);
        if(lec.weightedQ() > minPosQ*1e-2 && lec.weightedQ() < rangeMin) rangeMin=lec.weightedQ();
    }

    // shared log-spaced bins across all four channels, so shapes are
    // directly comparable
    const int nBins = 100;
    vector<double> edges(nBins+1);
    double logMin = std::log10(rangeMin*0.9);
    double logMax = std::log10(maxQ*1.0001);

    for(int b=0;b<=nBins;b++)
        edges[b] = std::pow(10.0, logMin + (logMax-logMin)*b/nBins);

    TCanvas *c = new TCanvas(("cQDist_"+tag).c_str(),"",1000,750);
    c->SetBatch(kTRUE);
    c->SetLogx();
    c->SetLogy();
    c->SetGrid();

    TLegend *leg = new TLegend(0.50,0.50,0.98,0.92);

    TH1D *h[5];
    double globalMax = 1;

    for(int c4=0; c4<5; c4++){

        int   color = (c4<4) ? kChannelColor[c4] : kWeightedColor;
        const char *name = (c4<4) ? kChannelName[c4] : kWeightedName;

        h[c4] = new TH1D(Form("hQDist_%d_%s",c4,tag.c_str()),
                          (fullTitle+";q;windows").c_str(), nBins, edges.data());

        for(double q : vals[c4]) if(q>0) h[c4]->Fill(q);

        h[c4]->SetStats(0);
        h[c4]->SetLineColor(color);
        h[c4]->SetLineWidth(2);
        if(c4==4) h[c4]->SetLineStyle(2); // dash the weighted one so it's distinguishable from k1 (same shape, larger values)

        if(h[c4]->GetMaximum()>globalMax) globalMax=h[c4]->GetMaximum();

        leg->AddEntry(h[c4],name,"l");
    }

    h[0]->GetYaxis()->SetRangeUser(0.5,globalMax*3);
    h[0]->Draw("HIST");
    for(int c4=1; c4<5; c4++) h[c4]->Draw("HIST SAME");

    TLine *lCap = new TLine(kMaxQCut,0.5,kMaxQCut,globalMax*3);
    lCap->SetLineColor(kGray+2);
    lCap->SetLineStyle(2);
    lCap->SetLineWidth(2);
    lCap->Draw();
    leg->AddEntry(lCap,Form("scan cap (%.4f)",kMaxQCut),"l");

    if(lec.valid){
        for(int i=0;i<4;i++){
            if(lec.q(i)<=0) continue;
            TLine *lLec = new TLine(lec.q(i),0.5,lec.q(i),globalMax*3);
            lLec->SetLineColor(kChannelColor[i]);
            lLec->SetLineStyle(3);
            lLec->SetLineWidth(2);
            lLec->Draw();
            leg->AddEntry(lLec,Form("LEC %s (%.5f)",kChannelName[i],lec.q(i)),"l");
        }
        if(lec.weightedQ()>0){
            TLine *lLecW = new TLine(lec.weightedQ(),0.5,lec.weightedQ(),globalMax*3);
            lLecW->SetLineColor(kWeightedColor);
            lLecW->SetLineStyle(3);
            lLecW->SetLineWidth(2);
            lLecW->Draw();
            leg->AddEntry(lLecW,Form("LEC weighted (%.5f)",lec.weightedQ()),"l");
        }
    }

    leg->Draw();

    string out = dir+"/qcut_distribution_"+tag+".pdf";
    c->SaveAs(out.c_str());
    delete c;

    cout << "  wrote " << out << endl;
}

void plotAllQCutScans(const string &label, const string &title){

    vector<int> exposures = discoverTaggedSubdirs(label,'e');

    if(exposures.empty()){
        cout << "[" << label << "] no e<N> subfolders found under '" << label << "/' -- skipping" << endl;
        return;
    }

    for(int exposure : exposures){

        string expDir = label+"/e"+to_string(exposure);

        vector<int> radii = discoverTaggedSubdirs(expDir,'r');

        if(radii.empty()){
            cout << "[" << label << "/e" << exposure << "] no r<N> subfolders found -- skipping" << endl;
            continue;
        }

        for(int radius : radii){
            plotQCutScanForCombo(label,title,exposure,radius);
            plotQDistribution(label,title,exposure,radius);
        }
    }
}

int PlotQ1CutScan(){

    gROOT->SetBatch(kTRUE);

    plotAllQCutScans("goodquads",  "good quads");
    plotAllQCutScans("notbadquads","not-bad quads");

    return 0;
}
