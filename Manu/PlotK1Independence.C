// Independence tests between multi-electron events and the other per-window
// counts, using the real windows.
//
// A study (IndepSpec) is defined by
//   - the tag: which hits tag a window -- k2 (2e hits), or k3+ = k3 + k4
//     (hits with >= 3e) -- and whether the window needs at least one
//     (">= 1") or exactly one ("= 1") of them;
//   - the tested count: the per-window count whose distribution is compared
//     (k1, or k2).
//
// Observed: distribution of the tested count over windows that pass the tag.
// Expected under independence: keep each window's real counts and real N,
// and take the tagging hits as uniform over its image's unmasked pixels at
// that image's measured rate lambda = n_tag/Npix. With x = N*lambda the
// window passes the tag with probability
//     p = 1 - exp(-x)        (at least one)
//     p = x * exp(-x)        (exactly one)
// and contributes p to its bin (exact expectation, no simulation). Per-bin
// variance of the expected count is sum p(1-p) (sum of independent
// Bernoullis).
//
// Needs ImageRates.C's output (imagerates/e<exposure>/) for n1..n4 and Npix --
// window sums can't be used for these because windows overlap.
//
// Caveat: windows overlap (grid step = radius), so neighbouring windows
// share events and are not independent -- the chi2 p-value below treats
// them as independent and is therefore optimistic (too small errors).

#include "PlotWindowStats.C"
#include "PlotStyle.C"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TColor.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegendEntry.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"

#include <map>
#include <tuple>

// (input-file index, runID, LTA, OHDU): run IDs restart between data
// campaigns, so the same (runID, LTA, OHDU) can be two different images when
// a dataset mixes campaigns. The file index n is the same number used in
// windowStatisticsManye_<label>_<n>.root and imageRates_<n>.root (position in
// the config's file list), which makes the key unique.
typedef std::tuple<int,int,int,int> ImgKey;

struct ImgRate{ Long64_t Npix=0; Long64_t n[4]={0,0,0,0}; };

// trailing "_<n>.root" of a per-input-file output name, -1 if absent
int fileIndexOf(const string &path){
    size_t dot = path.rfind(".root");
    size_t us  = path.rfind('_',dot);
    if(dot==string::npos || us==string::npos) return -1;
    string num = path.substr(us+1,dot-us-1);
    if(num.empty() || !std::all_of(num.begin(),num.end(),::isdigit)) return -1;
    return std::stoi(num);
}

map<ImgKey,ImgRate> loadImageRates(int exposure){

    map<ImgKey,ImgRate> rates;

    string dir = "imagerates/e"+to_string(exposure);

    TSystemDirectory sysDir(dir.c_str(),dir.c_str());
    TList *files = sysDir.GetListOfFiles();
    if(!files) return rates;

    for(auto obj : *files){

        string name = obj->GetName();
        if(name.rfind("imageRates_",0)!=0) continue;
        int idx = fileIndexOf(name);
        if(idx<0) continue;

        TFile f((dir+"/"+name).c_str());
        TTree *t = f.IsZombie() ? nullptr : (TTree*)f.Get("images");
        if(!t) continue;

        int run, lta, ohdu;
        Long64_t Npix, nc[4];
        t->SetBranchAddress("runID",&run);
        t->SetBranchAddress("LTA",&lta);
        t->SetBranchAddress("OHDU",&ohdu);
        t->SetBranchAddress("Npix",&Npix);
        for(int c=0;c<4;c++) t->SetBranchAddress(Form("n%d",c+1),&nc[c]);

        for(Long64_t i=0;i<t->GetEntries();i++){
            t->GetEntry(i);
            auto &r = rates[std::make_tuple(idx,run,lta,ohdu)];
            r.Npix += Npix;
            for(int c=0;c<4;c++) r.n[c] += nc[c];
        }
    }

    delete files;

    return rates;
}

// plots show the tested count up to this value (histograms and chi2 still
// use the full range)
const int kIndepPlotMax = 20;

struct IndepSpec{
    vector<int> tagChannels;   // indices into (k1,k2,k3,k4): {1} = k2, {2,3} = k3+k4
    bool        exactlyOne;    // tag count == 1 (true) or >= 1 (false)
    int         testVar;       // index of the tested count: 0 = k1, 1 = k2
    string      tagLatex;      // "k_{2}", "k_{3+}"
    string      rateLatex;     // numerator of lambda, e.g. "n_{3}+n_{4+}"
    string      lambdaLatex;   // "#lambda_{2}", "#lambda_{3+}"
    string      fileTag;       // "k2ge1", "k3pe1"
    string      threeCasesPrefix; // output name stem of the three-case page
};

struct IndepResult{ double totO=0, totE=0, totV=0, chi2=0, pval=1; int ndf=0; };

// Statistics and pages for one study, from per-window inputs:
//   vX   = tested count (k1 or k2) of each window
//   vP   = probability that the window passes the tag under independence
//   vSel = whether the window really passes it
// With draw = false only the numbers are computed (closure toys). A
// non-empty footnote replaces the default one on the test page.
IndepResult independencePages(const vector<int> &vX, const vector<double> &vP, const vector<char> &vSel,
                              const IndepSpec &sp, const string &dir, const string &tag,
                              const string &title, int exposure, int radius,
                              bool draw=true, const string &footnote=""){

    const string varName  = "k"+to_string(sp.testVar+1);
    const string varL     = "k_{"+to_string(sp.testVar+1)+"}";
    const string selL     = sp.tagLatex+(sp.exactlyOne ? " = 1" : " #geq 1");
    const string selPlain = sp.fileTag;
    const string pL       = sp.exactlyOne ? "N"+sp.lambdaLatex+" e^{-N"+sp.lambdaLatex+"}"
                                          : "1 - e^{-N"+sp.lambdaLatex+"}";
    const string htag     = tag+"_"+sp.fileTag+"_"+varName;

    if(draw) applyPlotStyle();

    int maxX = 0;
    for(int x : vX) maxX = std::max(maxX,x);

    int nb = maxX+1;

    TH1D *hObs = new TH1D(("hObs_"+htag).c_str(),"",nb,-0.5,maxX+0.5);
    TH1D *hExp = new TH1D(("hExp_"+htag).c_str(),"",nb,-0.5,maxX+0.5);
    vector<double> var(nb+2,0.0);

    for(size_t i=0;i<vX.size();i++){
        if(vSel[i]) hObs->Fill(vX[i]);
        hExp->Fill(vX[i],vP[i]);
        var[hExp->FindBin(vX[i])] += vP[i]*(1.0-vP[i]);
    }

    for(int b=1;b<=nb;b++) hExp->SetBinError(b,std::sqrt(var[b]));

    // chi2: merge consecutive bins (left to right) until each group has
    // expected >= 5, last short group merged into the previous one
    vector<double> gO, gE, gV;
    double accO=0, accE=0, accV=0;
    for(int b=1;b<=nb;b++){
        accO += hObs->GetBinContent(b);
        accE += hExp->GetBinContent(b);
        accV += var[b];
        if(accE>=5){ gO.push_back(accO); gE.push_back(accE); gV.push_back(accV); accO=accE=accV=0; }
    }
    if(accE>0){
        if(gE.empty()){ gO.push_back(accO); gE.push_back(accE); gV.push_back(accV); }
        else{ gO.back()+=accO; gE.back()+=accE; gV.back()+=accV; }
    }

    double chi2 = 0;
    for(size_t g=0;g<gE.size();g++) if(gV[g]>0) chi2 += (gO[g]-gE[g])*(gO[g]-gE[g])/gV[g];
    int ndf = (int)gE.size();
    double pval = TMath::Prob(chi2,ndf);

    double totO = hObs->Integral(), totE = hExp->Integral();
    double totV = 0; for(int b=1;b<=nb;b++) totV += var[b];

    double meanO = hObs->GetMean(), meanE = hExp->GetMean();

    IndepResult res;
    res.totO = totO; res.totE = totE; res.totV = totV;
    res.chi2 = chi2; res.ndf = ndf;   res.pval = pval;

    if(!draw){
        delete hObs;
        delete hExp;
        return res;
    }

    cout << Form("  observed %s windows = %.0f, expected = %.1f +- %.1f; chi2/ndf = %.1f/%d (p = %.3g)",
                 selPlain.c_str(),totO,totE,std::sqrt(totV),chi2,ndf,pval) << endl;

    //////////////////////////////////////////////////////////
    // plots
    //////////////////////////////////////////////////////////

    applyPlotStyle();

    double xHi = std::min(maxX,kIndepPlotMax)+0.5;

    string subtitle = title+"  #upoint  exposure "+to_string(exposure)+"  #upoint  r = "+to_string(radius);

    TH1D *hAll = new TH1D(("hAll_"+htag).c_str(),"",nb,-0.5,maxX+0.5);
    for(int k : vX) hAll->Fill(k);

    // styled copies, so the same encoding is reused in every panel:
    // real data = ink points, independence = blue line + light-blue band,
    // all windows = orange line
    auto obsPoints = [&](TH1D *h){
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.9);
        h->SetMarkerColor(gStInk);
        h->SetLineColor(gStInk);
        h->SetLineWidth(1);
    };
    auto expBand = [&](TH1D *h){
        h->SetFillColor(gStBlueBand);
        h->SetFillStyle(1001);
        h->SetLineColor(gStBlueBand);
        h->SetMarkerSize(0);
    };
    auto expLine = [&](TH1D *h){
        h->SetLineColor(gStBlue);
        h->SetLineWidth(2);
        h->SetFillStyle(0);
    };
    auto allLine = [&](TH1D *h){
        h->SetLineColor(gStOrange);
        h->SetLineWidth(2);
        h->SetFillStyle(0);
    };

    //////////////////////////////////////////////////////////
    // 1) observed vs independence, with ratio
    //////////////////////////////////////////////////////////

    TCanvas *c = new TCanvas(("cind_"+htag).c_str(),"",1000,820);
    c->SetBatch(kTRUE);

    TPad *pTop = new TPad("pTop","",0.0,0.32,1.0,1.0);
    TPad *pBot = new TPad("pBot","",0.0,0.0,1.0,0.32);
    pTop->SetBottomMargin(0.02);
    pTop->SetTopMargin(0.15);
    pBot->SetTopMargin(0.03);
    pBot->SetBottomMargin(0.34);
    pTop->Draw();
    pBot->Draw();

    pTop->cd();
    pTop->SetLogy();
    pTop->SetGridy();

    TH1D *bandE = (TH1D*)hExp->Clone(("bandE_"+htag).c_str());
    TH1D *lineE = (TH1D*)hExp->Clone(("lineE_"+htag).c_str());
    TH1D *ptsO  = (TH1D*)hObs->Clone(("ptsO_"+htag).c_str());
    expBand(bandE);
    expLine(lineE);
    obsPoints(ptsO);

    double ymax = std::max(hExp->GetMaximum(),hObs->GetMaximum());
    bandE->GetYaxis()->SetRangeUser(0.3,ymax*8);
    bandE->GetYaxis()->SetTitle("Windows");
    bandE->GetXaxis()->SetLabelSize(0);
    bandE->GetXaxis()->SetRangeUser(-0.5,xHi);
    bandE->Draw("E2");
    lineE->Draw("HIST SAME");
    ptsO->Draw("E1 X0 SAME");
    gPad->RedrawAxis();

    padHeader(varL+" of windows with "+selL, subtitle, 0.085);

    TLegend *leg = styledLegend(0.54,0.58,0.95,0.80);
    leg->AddEntry(ptsO,("real windows, "+selL).c_str(),"pe");
    leg->AddEntry(lineE,("independence #sum ("+pL+")").c_str(),"l");
    leg->AddEntry(bandE,"#pm1#sigma, #sqrt{#sum p(1-p)}","f");
    leg->Draw();

    // test summary on one line under the title, outside the data area
    TLatex st;
    st.SetNDC();
    st.SetTextFont(43);
    st.SetTextSize(16);
    st.SetTextColor(gStInk2);
    st.SetTextAlign(11);
    st.DrawLatex(pTop->GetLeftMargin(),1-pTop->GetTopMargin()+0.025,
        Form("observed %.0f,  expected %.1f #pm %.1f     #LT%s#GT  %.2f observed,  %.2f expected     #chi^{2}/ndf = %.1f / %d,  p = %.2g",
             totO,totE,std::sqrt(totV),varL.c_str(),meanO,meanE,chi2,ndf,pval));

    TLatex note;
    note.SetNDC();
    note.SetTextFont(43);
    note.SetTextSize(14);
    note.SetTextColor(gStMuted);
    note.SetTextAlign(11);
    note.DrawLatex(pTop->GetLeftMargin()+0.015,0.80,
        (footnote.empty() ? sp.lambdaLatex+" = ("+sp.rateLatex+")/N_{pix} per image  #upoint  #chi^{2} bins merged to E #geq 5  #upoint  windows overlap, so p is optimistic"
                      : footnote).c_str());

    pBot->cd();
    pBot->SetGridy();

    TH1D *rBand = (TH1D*)hExp->Clone(("rBand_"+htag).c_str());
    TH1D *rPts  = (TH1D*)hObs->Clone(("rPts_"+htag).c_str());
    for(int b=1;b<=nb;b++){
        double E = hExp->GetBinContent(b);
        double O = hObs->GetBinContent(b);
        if(E>=1){
            rBand->SetBinContent(b,1);
            rBand->SetBinError(b,std::sqrt(var[b])/E);
            rPts->SetBinContent(b,O/E);
            rPts->SetBinError(b,std::sqrt(O)/E);
        }
        else{
            // less than one window expected: ratio not meaningful, keep off the plot
            rBand->SetBinContent(b,-10); rBand->SetBinError(b,0);
            rPts->SetBinContent(b,-10);  rPts->SetBinError(b,0);
        }
    }
    expBand(rBand);
    obsPoints(rPts);

    rBand->GetYaxis()->SetRangeUser(0,3.2);
    rBand->GetYaxis()->SetNdivisions(503);
    rBand->GetYaxis()->SetTitle("real / indep.");
    rBand->GetXaxis()->SetTitle((varL+" per window").c_str());
    rBand->GetXaxis()->SetRangeUser(-0.5,xHi);
    rBand->GetXaxis()->SetTitleOffset(1.0);
    rBand->Draw("E2");

    TLine *one = new TLine(-0.5,1,xHi,1);
    one->SetLineColor(gStBlue);
    one->SetLineWidth(2);
    one->Draw();

    rPts->Draw("E1 X0 SAME");
    gPad->RedrawAxis();

    string out = dir+"/"+varName+"_independence_"+sp.fileTag+"_"+tag+".pdf";
    c->SaveAs(out.c_str());
    delete c;
    tightPdf(out);

    cout << "  wrote " << out << endl;

    //////////////////////////////////////////////////////////
    // 2) the three cases (all windows, real, independence), and P(tag) vs
    //    the tested count
    //////////////////////////////////////////////////////////

    TCanvas *c3 = new TCanvas(("c3_"+htag).c_str(),"",1400,520);
    c3->SetBatch(kTRUE);

    TPad *pHead = new TPad("pHead","",0.0,0.91,1.0,1.0);
    pHead->Draw();
    pHead->cd();
    {
        TLatex t;
        t.SetNDC();
        t.SetTextFont(43);
        t.SetTextSize(20);
        t.SetTextColor(gStInk);
        t.SetTextAlign(12);
        t.DrawLatex(0.035,0.5,(varL+" per window,  windows with "+selL+"  #upoint  "+subtitle).c_str());
    }
    c3->cd();

    TPad *pp[3];
    for(int i=0;i<3;i++){
        pp[i] = new TPad(Form("pp%d",i),"",i/3.0,0.0,(i+1)/3.0,0.91);
        pp[i]->SetTopMargin(0.09);
        pp[i]->SetLeftMargin(0.18);
        pp[i]->SetRightMargin(0.04);
        pp[i]->SetBottomMargin(0.14);
        pp[i]->Draw();
    }

    // panel 1: counts
    pp[0]->cd();
    gPad->SetLogy();
    gPad->SetGridy();
    {
        TH1D *a  = (TH1D*)hAll->Clone(("aAll_"+htag).c_str());
        TH1D *eb = (TH1D*)hExp->Clone(("aExpB_"+htag).c_str());
        TH1D *el = (TH1D*)hExp->Clone(("aExpL_"+htag).c_str());
        TH1D *o  = (TH1D*)hObs->Clone(("aObs_"+htag).c_str());
        allLine(a); expBand(eb); expLine(el); obsPoints(o);

        a->GetYaxis()->SetRangeUser(0.3,a->GetMaximum()*8);
        a->GetYaxis()->SetTitle("Windows");
        a->GetXaxis()->SetTitle((varL+" per window").c_str());
        a->GetXaxis()->SetRangeUser(-0.5,xHi);
        a->Draw("HIST");
        eb->Draw("E2 SAME");
        el->Draw("HIST SAME");
        a->Draw("HIST SAME");
        o->Draw("E1 X0 SAME");
        gPad->RedrawAxis();

        padHeader("Counts","");

        TLegend *l = styledLegend(0.42,0.70,0.95,0.89);
        TLegendEntry *la = l->AddEntry(a,Form("all windows (%.0f)",hAll->Integral()),"l");
        la->SetLineColor(gStOrange); la->SetLineWidth(2);
        l->AddEntry(o,Form("real, %s (%.0f)",selL.c_str(),totO),"pe");
        TLegendEntry *le = l->AddEntry(el,Form("independence (%.1f)",totE),"l");
        le->SetLineColor(gStBlue); le->SetLineWidth(2);
        l->Draw();
    }

    // panel 2: unit-area shapes
    pp[1]->cd();
    gPad->SetLogy();
    gPad->SetGridy();
    {
        TH1D *a  = (TH1D*)hAll->Clone(("nAll_"+htag).c_str());
        TH1D *eb = (TH1D*)hExp->Clone(("nExpB_"+htag).c_str());
        TH1D *el = (TH1D*)hExp->Clone(("nExpL_"+htag).c_str());
        TH1D *o  = (TH1D*)hObs->Clone(("nObs_"+htag).c_str());
        a->Scale(1.0/a->Integral());
        eb->Scale(1.0/eb->Integral());
        el->Scale(1.0/el->Integral());
        o->Scale(1.0/o->Integral());
        allLine(a); expBand(eb); expLine(el); obsPoints(o);

        double m = std::max({a->GetMaximum(),o->GetMaximum(),el->GetMaximum()});
        a->GetYaxis()->SetRangeUser(0.3/hAll->Integral(),m*8);
        a->GetYaxis()->SetTitle("Fraction of windows");
        a->GetXaxis()->SetTitle((varL+" per window").c_str());
        a->GetXaxis()->SetRangeUser(-0.5,xHi);
        a->Draw("HIST");
        eb->Draw("E2 SAME");
        el->Draw("HIST SAME");
        a->Draw("HIST SAME");
        o->Draw("E1 X0 SAME");
        gPad->RedrawAxis();

        padHeader("Normalized to unit area","");

        TLegend *l = styledLegend(0.42,0.70,0.95,0.89);
        TLegendEntry *la = l->AddEntry(a,Form("all windows  #LT%s#GT = %.2f",varL.c_str(),hAll->GetMean()),"l");
        la->SetLineColor(gStOrange); la->SetLineWidth(2);
        l->AddEntry(o,Form("real, %s  #LT%s#GT = %.2f",selL.c_str(),varL.c_str(),meanO),"pe");
        TLegendEntry *le = l->AddEntry(el,Form("independence  #LT%s#GT = %.2f",varL.c_str(),meanE),"l");
        le->SetLineColor(gStBlue); le->SetLineWidth(2);
        l->Draw();
    }

    // panel 3: P(tag | tested count) -- real = fraction of windows in the bin
    // that pass the tag (Clopper-Pearson error); independence = mean p over
    // the same windows (band = sqrt(sum p(1-p))/n). Only bins that contain
    // windows are drawn.
    pp[2]->cd();
    gPad->SetLogy();
    gPad->SetGridy();
    {
        vector<double> px, pex, poy, poeLo, poeHi, pey, pee;
        double pLo = 1, pHi = 0;
        for(int b=1;b<=nb;b++){
            double nAllB = hAll->GetBinContent(b);
            if(nAllB<=0) continue;
            if(hAll->GetBinCenter(b)>xHi) break;
            double f = hObs->GetBinContent(b)/nAllB;
            px.push_back(hAll->GetBinCenter(b));
            pex.push_back(0.5);
            int nTot = (int)nAllB, nPass = (int)hObs->GetBinContent(b);
            poy.push_back(f);
            poeLo.push_back(f-TEfficiency::ClopperPearson(nTot,nPass,0.683,false));
            poeHi.push_back(TEfficiency::ClopperPearson(nTot,nPass,0.683,true)-f);
            pey.push_back(hExp->GetBinContent(b)/nAllB);
            pee.push_back(std::sqrt(var[b])/nAllB);
            for(double v : {poy.back(),pey.back()}){
                if(v>0 && v<pLo) pLo=v;
                if(v>pHi) pHi=v;
            }
        }

        TH1F *fr = gPad->DrawFrame(-0.5,pLo*0.3,xHi,pHi*25);
        fr->GetXaxis()->SetTitle((varL+" per window").c_str());
        fr->GetYaxis()->SetTitle(("P("+selL+")").c_str());

        vector<double> zero(px.size(),0.0);

        TGraphErrors *gB = new TGraphErrors((int)px.size(),px.data(),pey.data(),pex.data(),pee.data());
        gB->SetFillColor(gStBlueBand);
        gB->SetLineColor(gStBlueBand);
        gB->Draw("2");

        TGraph *gL = new TGraph((int)px.size(),px.data(),pey.data());
        gL->SetLineColor(gStBlue);
        gL->SetLineWidth(2);
        gL->Draw("L");

        TGraphAsymmErrors *gO = new TGraphAsymmErrors((int)px.size(),px.data(),poy.data(),zero.data(),zero.data(),poeLo.data(),poeHi.data());
        gO->SetMarkerStyle(20);
        gO->SetMarkerSize(0.9);
        gO->SetMarkerColor(gStInk);
        gO->SetLineColor(gStInk);
        gO->Draw("P");
        gPad->RedrawAxis();

        padHeader("P("+selL+") vs "+varL,"");

        TLegend *l = styledLegend(0.21,0.74,0.75,0.89);
        l->AddEntry(gO,("real: fraction with "+selL+" (68% C-P)").c_str(),"pe");
        TLegendEntry *le = l->AddEntry(gL,("independence #LT"+pL+"#GT").c_str(),"l");
        le->SetLineColor(gStBlue); le->SetLineWidth(2);
        l->Draw();
    }

    string out3 = dir+"/"+sp.threeCasesPrefix+"_"+tag+".pdf";
    c3->SaveAs(out3.c_str());
    delete c3;
    tightPdf(out3);

    cout << "  wrote " << out3 << endl;

    return res;
}

void plotIndependence(const string &label, const string &title, int exposure, int radius,
                      const IndepSpec &sp){

    string dir = label+"/e"+to_string(exposure)+"/r"+to_string(radius);
    string tag = label+"_e"+to_string(exposure)+"_r"+to_string(radius);

    const string varName  = "k"+to_string(sp.testVar+1);
    const string selPlain = sp.fileTag;

    cout << "[" << tag << "] " << varName << " independence test (" << selPlain << ") ..." << endl;

    // before any histogram is created, so every object starts from the
    // shared style (the first call in a session otherwise kept ROOT's
    // default relative font sizes in the ratio panel)
    applyPlotStyle();

    map<ImgKey,ImgRate> rates = loadImageRates(exposure);

    if(rates.empty()){
        cout << "  no imagerates/e" << exposure << "/ output -- run ImageRates.C first" << endl;
        return;
    }

    TChain *ch = buildChain(dir,label);
    if(!ch->GetListOfFiles() || ch->GetListOfFiles()->GetEntries()==0){
        cout << "  no window files -- skipping" << endl;
        delete ch;
        return;
    }

    int run, lta, ohdu, N, k[4];
    ch->SetBranchAddress("runID",&run);
    ch->SetBranchAddress("LTA",&lta);
    ch->SetBranchAddress("OHDU",&ohdu);
    ch->SetBranchAddress("N",&N);
    for(int c=0;c<4;c++) ch->SetBranchAddress(Form("k%d",c+1),&k[c]);

    Long64_t n = ch->GetEntries();

    vector<int> vX; vX.reserve(n);
    vector<double> vP; vP.reserve(n);
    vector<char> vSel; vSel.reserve(n);

    Long64_t nMissing = 0;
    int maxX = 0;

    int curTree = -1, curIdx = -1;

    for(Long64_t i=0;i<n;i++){
        ch->GetEntry(i);

        if(ch->GetTreeNumber()!=curTree){
            curTree = ch->GetTreeNumber();
            curIdx  = fileIndexOf(ch->GetFile()->GetName());
        }

        auto it = rates.find(std::make_tuple(curIdx,run,lta,ohdu));
        if(it==rates.end() || it->second.Npix<=0){ nMissing++; continue; }

        Long64_t nTag = 0;
        int kTag = 0;
        for(int c : sp.tagChannels){ nTag += it->second.n[c]; kTag += k[c]; }

        double x = N*(double)nTag/it->second.Npix;

        vX.push_back(k[sp.testVar]);
        vP.push_back(sp.exactlyOne ? x*std::exp(-x) : 1.0-std::exp(-x));
        vSel.push_back(sp.exactlyOne ? kTag==1 : kTag>=1);

        if(k[sp.testVar]>maxX) maxX=k[sp.testVar];
    }

    delete ch;

    if(nMissing>0)
        cout << "  WARNING: " << nMissing << " windows had no matching image rate and were skipped" << endl;

    independencePages(vX,vP,vSel,sp,dir,tag,title,exposure,radius);
}

// k1 of windows with at least one 2e hit (the original study)
IndepSpec specK2ge1(){
    return { {1}, false, 0, "k_{2}", "n_{2}", "#lambda_{2}", "k2ge1", "k1_threecases" };
}

// k1 or k2 of windows with exactly one hit of >= 3e
IndepSpec specK3pe1(int testVar){
    string v = "k"+to_string(testVar+1);
    return { {2,3}, true, testVar, "k_{3+}", "n_{3}+n_{4+}", "#lambda_{3+}", "k3pe1", v+"_threecases_k3pe1" };
}

void plotIndependenceForRadius(const string &label, const string &title, int radius, const IndepSpec &sp){

    for(int exposure : discoverTaggedSubdirs(label,'e')){
        vector<int> radii = discoverTaggedSubdirs(label+"/e"+to_string(exposure),'r');
        if(std::find(radii.begin(),radii.end(),radius)==radii.end()) continue;
        plotIndependence(label,title,exposure,radius,sp);
    }
}

int PlotK1Independence(){

    gROOT->SetBatch(kTRUE);

    vector<IndepSpec> specs = { specK2ge1(), specK3pe1(0), specK3pe1(1) };

    for(auto &sp : specs){
        plotIndependenceForRadius("goodquads",  "good quads",   30, sp);
        plotIndependenceForRadius("notbadquads","not-bad quads",30, sp);
    }

    return 0;
}
