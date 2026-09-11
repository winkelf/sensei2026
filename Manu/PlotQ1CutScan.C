// For each (category, exposure, radius), scans a q1 = k1/N cut from 0 up to
// min(kMaxQ1Cut, highest q1 value observed across windows) -- capped so a
// handful of near-fully-masked windows (e.g. N=1, k1=1, q1=1) don't swamp
// the scan -- keeping only windows with q1 < q1_cut at each step, and plots:
//   1) k1, k2, k3, k4 (summed over kept windows) vs. q1_cut
//   2) N (summed over kept windows) vs. q1_cut
//   3) k_i/N (the same two curves, divided) vs. q1_cut
// Each plot also overlays a single reference value computed once with the
// LEC ("Local Event Cut") mask bit added to the standard mask -- that
// selection isn't affected by the q1 cut, so it's just one number per
// channel, drawn as a dashed horizontal line.
//
// Requires MaskScanManye.C's standard output (<category>/e<exp>/r<radius>/)
// and, for the reference lines, MaskScanManyeLEC.C's output
// (<category>/lec/e<exp>/r<radius>/).

#include "PlotWindowStats.C"
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"

const int    kNCutSteps = 300;
const double kMaxQ1Cut  = 0.02;
const int   kChannelColor[4] = { kBlue+1, kRed+1, kGreen+2, kMagenta+1 };
const char *kChannelName[4]  = { "k_{1}", "k_{2}", "k_{3}", "k_{4+}" };

struct LecTotals{
    bool    valid = false;
    Long64_t sumN  = 0;
    Long64_t sumK[4] = {0,0,0,0};
};

// sums k1-k4 and N over every window in the LEC-mask output for this
// (category, exposure, radius) -- no q1 cut applied, computed once
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

// draws a dashed horizontal reference line spanning the pad's current
// x-range at y=value, in the given color, and adds it to the legend
TLine* referenceLine(double xmin, double xmax, double value, int color, TLegend *leg, const string &label){

    TLine *l = new TLine(xmin,value,xmax,value);
    l->SetLineColor(color);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();

    if(leg) leg->AddEntry(l,label.c_str(),"l");

    return l;
}

void plotQ1CutScanForCombo(const string &label, const string &title, int exposure, int radius){

    string dir = label+"/e"+to_string(exposure)+"/r"+to_string(radius);
    string tag = label+"_e"+to_string(exposure)+"_r"+to_string(radius);
    string fullTitle = title+", exposure="+to_string(exposure)+", r="+to_string(radius);

    cout << "[" << tag << "] q1-cut scan from '" << dir << "/' ..." << endl;

    TChain *ch = buildChain(dir,label);

    Long64_t nFiles = ch->GetListOfFiles() ? ch->GetListOfFiles()->GetEntries() : 0;

    if(nFiles==0){
        cout << "  no windowStatisticsManye_" << label << "_*.root files found -- skipping" << endl;
        delete ch;
        return;
    }

    int N,k1,k2,k3,k4;
    double r1;

    ch->SetBranchAddress("N", &N);
    ch->SetBranchAddress("k1",&k1);
    ch->SetBranchAddress("k2",&k2);
    ch->SetBranchAddress("k3",&k3);
    ch->SetBranchAddress("k4",&k4);
    ch->SetBranchAddress("r1",&r1);

    Long64_t n = ch->GetEntries();

    // pull the whole (small) per-window table into memory once
    vector<int>    vN(n), vK1(n), vK2(n), vK3(n), vK4(n);
    vector<double> vQ1(n);

    double maxQ1 = 0;

    for(Long64_t i=0;i<n;i++){
        ch->GetEntry(i);
        vN[i]=N; vK1[i]=k1; vK2[i]=k2; vK3[i]=k3; vK4[i]=k4; vQ1[i]=r1;
        if(r1>maxQ1) maxQ1=r1;
    }

    delete ch;

    // a handful of near-fully-masked windows (e.g. N=1, k1=1) can push the
    // true max q1 up to 1, which would swamp a scan meant to show the
    // informative low-q1 region -- cap the scan at kMaxQ1Cut instead
    double scanMax = std::min(maxQ1, kMaxQ1Cut);

    cout << "  " << n << " windows, max q1 = " << maxQ1
         << " (scanning up to " << scanMax << ")" << endl;

    // scan the cut from 0 to scanMax
    vector<double> cutVals(kNCutSteps+1);
    vector<double> sumK[4];
    vector<double> sumNv(kNCutSteps+1);

    for(int c=0;c<4;c++) sumK[c].assign(kNCutSteps+1,0);

    for(int s=0; s<=kNCutSteps; s++){

        double cut = scanMax * s / kNCutSteps;
        cutVals[s] = cut;

        Long64_t k[4] = {0,0,0,0};
        Long64_t Ns = 0;

        for(Long64_t i=0;i<n;i++){

            if(vQ1[i] >= cut) continue;

            k[0]+=vK1[i]; k[1]+=vK2[i]; k[2]+=vK3[i]; k[3]+=vK4[i];
            Ns+=vN[i];
        }

        for(int c=0;c<4;c++) sumK[c][s] = (double)k[c];
        sumNv[s] = (double)Ns;
    }

    LecTotals lec = loadLecTotals(label,exposure,radius);

    //////////////////////////////////////////////////////////
    // Plot 1: k1-k4 vs q1 cut
    //////////////////////////////////////////////////////////
    {
        TCanvas *c1 = new TCanvas(("ck1cut_"+tag).c_str(),"",1000,750);
        c1->SetBatch(kTRUE);
        c1->SetLogy();
        c1->SetGrid();

        TGraph *g[4];
        double globalMax = 1;

        for(int i=0;i<4;i++){
            g[i] = new TGraph(kNCutSteps+1, cutVals.data(), sumK[i].data());
            g[i]->SetLineColor(kChannelColor[i]);
            g[i]->SetLineWidth(2);
            for(double v : sumK[i]) if(v>globalMax) globalMax=v;
        }

        g[0]->SetTitle((fullTitle+";q_{1} cut;events kept (k)").c_str());
        g[0]->GetYaxis()->SetRangeUser(0.5, globalMax*3);
        g[0]->Draw("AL");
        for(int i=1;i<4;i++) g[i]->Draw("L SAME");

        TLegend *leg = new TLegend(0.72,0.15,0.98,0.45);
        for(int i=0;i<4;i++) leg->AddEntry(g[i],kChannelName[i],"l");

        if(lec.valid)
            for(int i=0;i<4;i++)
                referenceLine(0,scanMax,(double)lec.sumK[i],kChannelColor[i],leg,
                              string(kChannelName[i])+" (LEC)");

        leg->Draw();

        string out = dir+"/q1cut_k_scan_"+tag+".pdf";
        c1->SaveAs(out.c_str());
        delete c1;

        cout << "  wrote " << out << endl;
    }

    //////////////////////////////////////////////////////////
    // Plot 2: N vs q1 cut
    //////////////////////////////////////////////////////////
    {
        TCanvas *c2 = new TCanvas(("cNcut_"+tag).c_str(),"",1000,750);
        c2->SetBatch(kTRUE);
        c2->SetLogy();
        c2->SetGrid();

        TGraph *gN = new TGraph(kNCutSteps+1, cutVals.data(), sumNv.data());
        gN->SetLineColor(kAzure+2);
        gN->SetLineWidth(2);
        gN->SetTitle((fullTitle+";q_{1} cut;unmasked pixels kept (N)").c_str());

        double maxN = *std::max_element(sumNv.begin(),sumNv.end());
        gN->GetYaxis()->SetRangeUser(1, maxN*3);
        gN->Draw("AL");

        TLegend *leg = new TLegend(0.65,0.15,0.98,0.30);
        leg->AddEntry(gN,"N (q_{1} cut scan)","l");

        if(lec.valid)
            referenceLine(0,scanMax,(double)lec.sumN,kAzure+2,leg,"N (LEC)");

        leg->Draw();

        string out = dir+"/q1cut_N_scan_"+tag+".pdf";
        c2->SaveAs(out.c_str());
        delete c2;

        cout << "  wrote " << out << endl;
    }

    //////////////////////////////////////////////////////////
    // Plot 3: k_i / N (the ratio of the two plots above) vs q1 cut
    //////////////////////////////////////////////////////////
    {
        TCanvas *c3 = new TCanvas(("cRcut_"+tag).c_str(),"",1000,750);
        c3->SetBatch(kTRUE);
        c3->SetLogy();
        c3->SetGrid();

        // skip cut=0 (N=0, undefined ratio) when building each ratio graph
        vector<double> rc, rq[4];

        for(int s=0;s<=kNCutSteps;s++){
            if(sumNv[s]<=0) continue;
            rc.push_back(cutVals[s]);
            for(int i=0;i<4;i++)
                rq[i].push_back(sumK[i][s]/sumNv[s]);
        }

        TGraph *g[4];
        double globalMax = 1e-9;

        for(int i=0;i<4;i++){
            g[i] = new TGraph((int)rc.size(), rc.data(), rq[i].data());
            g[i]->SetLineColor(kChannelColor[i]);
            g[i]->SetLineWidth(2);
            for(double v : rq[i]) if(v>globalMax) globalMax=v;
        }

        g[0]->SetTitle((fullTitle+";q_{1} cut;q_{i} = k_{i}/N (kept windows)").c_str());
        g[0]->GetYaxis()->SetRangeUser(globalMax*1e-6, globalMax*3);
        g[0]->Draw("AL");
        for(int i=1;i<4;i++) g[i]->Draw("L SAME");

        TLegend *leg = new TLegend(0.72,0.65,0.98,0.95);
        for(int i=0;i<4;i++) leg->AddEntry(g[i],("q "+string(kChannelName[i])).c_str(),"l");

        if(lec.valid && lec.sumN>0)
            for(int i=0;i<4;i++)
                referenceLine(0,scanMax,(double)lec.sumK[i]/lec.sumN,kChannelColor[i],leg,
                              "q "+string(kChannelName[i])+" (LEC)");

        leg->Draw();

        string out = dir+"/q1cut_ratio_scan_"+tag+".pdf";
        c3->SaveAs(out.c_str());
        delete c3;

        cout << "  wrote " << out << endl;
    }
}

void plotAllQ1CutScans(const string &label, const string &title){

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

        for(int radius : radii)
            plotQ1CutScanForCombo(label,title,exposure,radius);
    }
}

int PlotQ1CutScan(){

    gROOT->SetBatch(kTRUE);

    plotAllQ1CutScans("goodquads",  "good quads");
    plotAllQ1CutScans("notbadquads","not-bad quads");

    return 0;
}
