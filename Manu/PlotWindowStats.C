#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

#include "TROOT.h"
#include "TSystemDirectory.h"
#include "TList.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPaveText.h"

using namespace std;

// finds "r<N>" subfolders of <label>/ (one per radius MaskScanManye.C
// produced) and returns the radii found, sorted ascending
vector<int> discoverRadii(const string &label){

    vector<int> radii;

    TSystemDirectory sysDir(label.c_str(),label.c_str());
    TList *files = sysDir.GetListOfFiles();

    if(!files) return radii;

    for(auto obj : *files){

        TSystemFile *f = (TSystemFile*)obj;

        if(!f->IsDirectory()) continue;

        string name = f->GetName();

        if(name.size()<2 || name[0]!='r') continue;

        string digits = name.substr(1);

        if(digits.empty() || !std::all_of(digits.begin(),digits.end(),::isdigit))
            continue;

        radii.push_back(std::stoi(digits));
    }

    delete files;

    std::sort(radii.begin(),radii.end());

    return radii;
}

// running totals across all windows in a category
struct Totals{
    Long64_t nWindows = 0;
    Long64_t sumN     = 0;
    Long64_t sumK1    = 0;
    Long64_t sumK2    = 0;
    Long64_t sumK3    = 0;
    Long64_t sumK4    = 0;
};

// chains the "windows" tree from every windowStatisticsManye_<label>_*.root
// file found in <dir> (does not touch the raw hitSumm/calPixTree data --
// this only re-reads the small per-window summaries MaskScanManye.C wrote)
TChain* buildChain(const string &dir, const string &label){

    TChain *ch = new TChain("windows");

    TSystemDirectory sysDir(dir.c_str(),dir.c_str());
    TList *files = sysDir.GetListOfFiles();

    if(!files){
        cout << "  no such directory: " << dir << endl;
        return ch;
    }

    string prefix = "windowStatisticsManye_"+label+"_";

    for(auto obj : *files){

        TSystemFile *f = (TSystemFile*)obj;

        if(f->IsDirectory()) continue;

        string name = f->GetName();

        if(name.size()<prefix.size()+5) continue;
        if(name.compare(0,prefix.size(),prefix)!=0) continue;
        if(name.substr(name.size()-5)!=".root") continue;

        ch->Add((dir+"/"+name).c_str());
    }

    delete files;

    return ch;
}

void styleAndDraw(TH1D *h, const char *xTitle){

    h->SetStats(0);
    h->SetLineColor(kAzure+2);
    h->GetXaxis()->SetTitle(xTitle);
    h->GetYaxis()->SetTitle("Windows");
    h->Draw("HIST");
    gPad->SetLogy();
}

void plotCategory(const string &label, const string &title, int radius){

    string dir = label+"/r"+to_string(radius);
    string tag = label+"_r"+to_string(radius);
    string fullTitle = title+", r="+to_string(radius);

    cout << "[" << tag << "] scanning '" << dir << "/' ..." << endl;

    TChain *ch = buildChain(dir,label);

    Long64_t nFiles = ch->GetListOfFiles() ? ch->GetListOfFiles()->GetEntries() : 0;

    if(nFiles==0){
        cout << "  no windowStatisticsManye_" << label << "_*.root files found -- skipping" << endl;
        delete ch;
        return;
    }

    cout << "  chained " << nFiles << " files" << endl;

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

    TH1D hN (("hN_" +tag).c_str(), ("N per window ("       +fullTitle+")").c_str(),100,0,3000);
    TH1D hK1(("hK1_"+tag).c_str(), ("k_{1} per window ("   +fullTitle+")").c_str(), 11,-0.5,10.5);
    TH1D hK2(("hK2_"+tag).c_str(), ("k_{2} per window ("   +fullTitle+")").c_str(), 11,-0.5,10.5);
    TH1D hK3(("hK3_"+tag).c_str(), ("k_{3} per window ("   +fullTitle+")").c_str(), 11,-0.5,10.5);
    TH1D hK4(("hK4_"+tag).c_str(), ("k_{4+} per window ("  +fullTitle+")").c_str(), 11,-0.5,10.5);
    TH1D hQ1(("hQ1_"+tag).c_str(), ("q_{1} = k_{1}/N ("    +fullTitle+")").c_str(), 50,0,0.10);
    TH1D hQ2(("hQ2_"+tag).c_str(), ("q_{2} = k_{2}/N ("    +fullTitle+")").c_str(), 50,0,0.05);
    TH1D hQ3(("hQ3_"+tag).c_str(), ("q_{3} = k_{3}/N ("    +fullTitle+")").c_str(), 50,0,0.02);
    TH1D hQ4(("hQ4_"+tag).c_str(), ("q_{4+} = k_{4+}/N ("  +fullTitle+")").c_str(), 50,0,0.02);

    Totals tot;

    Long64_t n = ch->GetEntries();

    for(Long64_t i=0;i<n;i++){

        ch->GetEntry(i);

        hN.Fill(N);
        hK1.Fill(k1);
        hK2.Fill(k2);
        hK3.Fill(k3);
        hK4.Fill(k4);
        hQ1.Fill(r1);
        hQ2.Fill(r2);
        hQ3.Fill(r3);
        hQ4.Fill(r4);

        tot.nWindows++;
        tot.sumN  += N;
        tot.sumK1 += k1;
        tot.sumK2 += k2;
        tot.sumK3 += k3;
        tot.sumK4 += k4;
    }

    cout << "  " << tot.nWindows << " windows" << endl;

    double meanN  = tot.nWindows>0 ? (double)tot.sumN /tot.nWindows : 0;
    double meanQ1 = tot.sumN>0     ? (double)tot.sumK1/tot.sumN     : 0;
    double meanQ2 = tot.sumN>0     ? (double)tot.sumK2/tot.sumN     : 0;
    double meanQ3 = tot.sumN>0     ? (double)tot.sumK3/tot.sumN     : 0;
    double meanQ4 = tot.sumN>0     ? (double)tot.sumK4/tot.sumN     : 0;

    auto makeStatsBox = [&](){
        TPaveText *pt = new TPaveText(0.02,0.02,0.98,0.98,"NDC");
        pt->SetTextAlign(12);
        pt->SetBorderSize(0);
        pt->SetFillColor(0);
        pt->AddText(Form("Category: %s",fullTitle.c_str()));
        pt->AddText(Form("Windows: %lld",(Long64_t)tot.nWindows));
        pt->AddText(Form("<N> = %.1f",meanN));
        pt->AddText(Form("Total k_{1} = %lld    <q_{1}> = %.2e",(Long64_t)tot.sumK1,meanQ1));
        pt->AddText(Form("Total k_{2} = %lld    <q_{2}> = %.2e",(Long64_t)tot.sumK2,meanQ2));
        pt->AddText(Form("Total k_{3} = %lld    <q_{3}> = %.2e",(Long64_t)tot.sumK3,meanQ3));
        pt->AddText(Form("Total k_{4+} = %lld    <q_{4+}> = %.2e",(Long64_t)tot.sumK4,meanQ4));
        return pt;
    };

    //////////////////////////////////////////////////////////
    // k statistics summary -- raw counts per window
    //////////////////////////////////////////////////////////

    TCanvas *ck = new TCanvas(("ck_"+tag).c_str(),"",1500,900);
    ck->SetBatch(kTRUE);
    ck->Divide(3,2);

    TH1D *kPlots[4] = { &hK1, &hK2, &hK3, &hK4 };

    for(int i=0;i<4;i++){
        ck->cd(i+1);
        styleAndDraw(kPlots[i],"k (events/window)");
    }

    ck->cd(5);
    styleAndDraw(&hN,"N (unmasked pixels/window)");

    ck->cd(6);
    makeStatsBox()->Draw();

    string kOut = dir+"/k_statistics_summary_"+tag+".pdf";
    ck->SaveAs(kOut.c_str());
    delete ck;

    //////////////////////////////////////////////////////////
    // q statistics summary -- ratios k/N
    //////////////////////////////////////////////////////////

    TCanvas *cq = new TCanvas(("cq_"+tag).c_str(),"",1500,900);
    cq->SetBatch(kTRUE);
    cq->Divide(3,2);

    cq->cd(1);
    styleAndDraw(&hN,"N (unmasked pixels/window)");

    TH1D *qPlots[4] = { &hQ1, &hQ2, &hQ3, &hQ4 };

    for(int i=0;i<4;i++){
        cq->cd(i+2);
        styleAndDraw(qPlots[i],"q = k/N");
    }

    cq->cd(6);
    makeStatsBox()->Draw();

    string qOut = dir+"/q_statistics_summary_"+tag+".pdf";
    cq->SaveAs(qOut.c_str());
    delete cq;

    cout << "  wrote " << kOut << " and " << qOut << endl;

    delete ch;
}

void plotAllRadii(const string &label, const string &title){

    vector<int> radii = discoverRadii(label);

    if(radii.empty()){
        cout << "[" << label << "] no r<N> subfolders found under '" << label << "/' -- skipping" << endl;
        return;
    }

    for(int radius : radii)
        plotCategory(label,title,radius);
}

int PlotWindowStats(){

    gROOT->SetBatch(kTRUE);

    plotAllRadii("goodquads",  "good quads");
    plotAllRadii("notbadquads","not-bad quads");

    return 0;
}
