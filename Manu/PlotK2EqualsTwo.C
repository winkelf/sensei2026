// Selection studies on the per-window counts:
//   - windows with exactly v tagging hits (k2 == 1, k2 == 2, or
//     k3+ == 1 with k3+ = k3 + k4): distributions of N and of every count
//     the selection does not fix, with their rates q = k/N;
//   - all windows (no cut).
//
// Restricted to a single radius at a time (r30), other radii can be added
// the same way.
//
// #include's PlotWindowStats.C to reuse buildChain()/discoverTaggedSubdirs().

#include <functional>

#include "PlotWindowStats.C"
#include "PlotStyle.C"

const char *kSelKTitle[4] = {"k_{1}","k_{2}","k_{3}","k_{4+}"};
const char *kSelQTitle[4] = {"q_{1} = k_{1}/N","q_{2} = k_{2}/N","q_{3} = k_{3}/N","q_{4+} = k_{4+}/N"};
const char *kSelKHead[4]  = {"1e hits per window","2e hits per window","3e hits per window","#geq4e hits per window"};
const char *kSelQHead[4]  = {"1e rate per live pixel","2e rate per live pixel","3e rate per live pixel","#geq4e rate per live pixel"};

// windows whose tag count (sum of k over tagCh) equals tagVal
void plotTagSelection(const string &label, const string &title, int exposure, int radius,
                      const vector<int> &tagCh, int tagVal, const string &tagLatex,
                      const string &fileStem){

    string dir = label+"/e"+to_string(exposure)+"/r"+to_string(radius);
    string tag = label+"_e"+to_string(exposure)+"_r"+to_string(radius);
    string htag = fileStem+"_"+tag;

    cout << "[" << tag << "] " << fileStem << " selection from '" << dir << "/' ..." << endl;

    TChain *ch = buildChain(dir,label);

    Long64_t nFiles = ch->GetListOfFiles() ? ch->GetListOfFiles()->GetEntries() : 0;

    if(nFiles==0){
        cout << "  no windowStatisticsManye_" << label << "_*.root files found -- skipping" << endl;
        delete ch;
        return;
    }

    // counts not fixed by the selection
    vector<int> freeCh;
    for(int c=0;c<4;c++) if(std::find(tagCh.begin(),tagCh.end(),c)==tagCh.end()) freeCh.push_back(c);
    const int nF = freeCh.size();

    int N, k[4];
    double r[4];
    ch->SetBranchAddress("N",&N);
    for(int c=0;c<4;c++){
        ch->SetBranchAddress(Form("k%d",c+1),&k[c]);
        ch->SetBranchAddress(Form("r%d",c+1),&r[c]);
    }

    Long64_t n = ch->GetEntries();

    // first pass: collect the selected windows' values so histogram ranges
    // can be set to the full observed range instead of a fixed guess
    vector<int> selN;
    vector<vector<int>>    selK(nF);
    vector<vector<double>> selQ(nF);
    vector<Long64_t> tagSplit(tagCh.size(),0);   // selected windows with the tag hit(s) in each tag channel

    for(Long64_t i=0;i<n;i++){

        ch->GetEntry(i);

        int kTag = 0;
        for(int c : tagCh) kTag += k[c];
        if(kTag!=tagVal) continue;

        selN.push_back(N);
        for(int j=0;j<nF;j++){
            selK[j].push_back(k[freeCh[j]]);
            selQ[j].push_back(r[freeCh[j]]);
        }
        for(size_t j=0;j<tagCh.size();j++) if(k[tagCh[j]]>0) tagSplit[j]++;
    }

    delete ch;

    Long64_t nSel = (Long64_t)selN.size();

    cout << "  " << nSel << " / " << n << " windows selected" << endl;

    if(nSel==0){
        cout << "  no windows selected -- skipping plot" << endl;
        return;
    }

    int maxN = 0;
    Long64_t sumN = 0;
    for(int v : selN){ maxN = std::max(maxN,v); sumN += v; }

    vector<int> maxK(nF,0);
    vector<double> maxQ(nF,0);
    vector<Long64_t> sumK(nF,0);
    for(int j=0;j<nF;j++){
        for(Long64_t i=0;i<nSel;i++){
            maxK[j] = std::max(maxK[j],selK[j][i]);
            maxQ[j] = std::max(maxQ[j],selQ[j][i]);
            sumK[j] += selK[j][i];
        }
        // guard against a degenerate (zero-width) axis when every selected
        // window has the same (often zero) value for a ratio
        if(maxQ[j]<=0) maxQ[j] = 1e-6;
    }

    // integer-valued counts: one bin per integer, centered on that integer
    // (bin [k-0.5,k+0.5] holds the count k), over the full observed range
    TH1D *hN = new TH1D(("hN_sel_"+htag).c_str(),"",60,0,maxN+1);
    vector<TH1D*> hK(nF), hQ(nF);
    for(int j=0;j<nF;j++){
        hK[j] = new TH1D(Form("hK%d_sel_%s",freeCh[j]+1,htag.c_str()),"",maxK[j]+1,-0.5,maxK[j]+0.5);
        hQ[j] = new TH1D(Form("hQ%d_sel_%s",freeCh[j]+1,htag.c_str()),"",50,0,maxQ[j]*1.001);
    }
    for(Long64_t i=0;i<nSel;i++){
        hN->Fill(selN[i]);
        for(int j=0;j<nF;j++){ hK[j]->Fill(selK[j][i]); hQ[j]->Fill(selQ[j][i]); }
    }

    double meanN = (double)sumN/nSel;

    applyPlotStyle();

    const int ncol = 1+nF;
    const int CW = 400*ncol, CH = 780;
    TCanvas *c = new TCanvas(("csel_"+htag).c_str(),"",CW,CH);
    c->SetBatch(kTRUE);

    string fixedNote;
    if(tagCh.size()==1)
        fixedNote = tagLatex+" is fixed by the selection, so it is not shown";
    else{
        fixedNote = tagLatex+" = ";
        for(size_t j=0;j<tagCh.size();j++) fixedNote += string(j ? " + " : "")+kSelKTitle[tagCh[j]];
        fixedNote += " is fixed by the selection (windows with ";
        for(size_t j=0;j<tagCh.size();j++)
            fixedNote += Form("%s%s > 0: %lld",j ? ", " : "",kSelKTitle[tagCh[j]],tagSplit[j]);
        fixedNote += ")";
    }

    string context = title+"  #upoint  exposure "+to_string(exposure)+"  #upoint  r = "+to_string(radius);
    double used = canvasHeader(c,CH,Form("Windows with %s = %d",tagLatex.c_str(),tagVal),context,
        Form("%lld of %lld windows (%.2g%%)     #LTN#GT = %.0f live pixels     %s",
             (Long64_t)nSel,(Long64_t)n,100.0*nSel/n,meanN,fixedNote.c_str()));

    vector<TPad*> pads = gridPads(ncol,2,1-used-0.01,"psel_"+htag+"_");

    pads[0]->cd(); drawDistribution(hN,nColor(),nFill(),"N"); padHeader("Live pixels per window","");
    for(int j=0;j<nF;j++){
        int cc = freeCh[j];
        pads[1+j]->cd();
        drawDistribution(hK[j],chanColor(cc),chanFill(cc),kSelKTitle[cc]);
        padHeader(kSelKHead[cc],"");
        pads[ncol+j]->cd();
        drawDistribution(hQ[j],chanColor(cc),chanFill(cc),kSelQTitle[cc]);
        padHeader(kSelQHead[cc],"");
    }

    // summary table in the last cell
    pads[2*ncol-1]->cd();
    {
        TLatex t;
        t.SetNDC();
        t.SetTextFont(43);
        t.SetTextAlign(11);
        t.SetTextSize(18);
        t.SetTextColor(gStInk);
        t.DrawLatex(0.08,0.80,"Totals over the selected windows");
        t.SetTextSize(17);
        t.SetTextColor(gStInk2);
        t.DrawLatex(0.08,0.66,"channel");
        t.DrawLatex(0.42,0.66,"total k");
        t.DrawLatex(0.66,0.66,"#sum k / #sum N");
        for(int j=0;j<nF;j++){
            int cc = freeCh[j];
            double y = 0.53-0.12*j;
            TLine *sw = new TLine(0.08,y+0.012,0.14,y+0.012);
            sw->SetNDC(); sw->SetLineColor(chanColor(cc)); sw->SetLineWidth(3); sw->Draw();
            t.SetTextColor(gStInk);
            t.DrawLatex(0.17,y,kSelKTitle[cc]);
            t.DrawLatex(0.42,y,Form("%lld",sumK[j]));
            t.DrawLatex(0.66,y,Form("%.2e",sumN>0 ? (double)sumK[j]/sumN : 0.0));
        }
    }

    string out = dir+"/"+fileStem+"_selection_"+tag+".pdf";
    c->SaveAs(out.c_str());
    delete c;
    tightPdf(out);
    cout << "  wrote " << out << endl;
}

// windows with exactly k2Val 2e hits
void plotK2EqSelection(const string &label, const string &title, int exposure, int radius, int k2Val=2){
    plotTagSelection(label,title,exposure,radius,{1},k2Val,"k_{2}","k2eq"+to_string(k2Val));
}

// windows with exactly v hits of >= 3e
void plotK3pEqSelection(const string &label, const string &title, int exposure, int radius, int v=1){
    plotTagSelection(label,title,exposure,radius,{2,3},v,"k_{3+}","k3pe"+to_string(v));
}

void forEachExposure(const string &label, int radius, const std::function<void(int)> &f){

    vector<int> exposures = discoverTaggedSubdirs(label,'e');

    if(exposures.empty()){
        cout << "[" << label << "] no e<N> subfolders found under '" << label << "/' -- skipping" << endl;
        return;
    }

    for(int exposure : exposures){
        vector<int> radii = discoverTaggedSubdirs(label+"/e"+to_string(exposure),'r');
        if(std::find(radii.begin(),radii.end(),radius)==radii.end()){
            cout << "[" << label << "/e" << exposure << "] no r" << radius << " subfolder -- skipping" << endl;
            continue;
        }
        f(exposure);
    }
}

void plotK2SelectionForRadius(const string &label, const string &title, int radius, int k2Val=2){
    forEachExposure(label,radius,[&](int e){ plotK2EqSelection(label,title,e,radius,k2Val); });
}

void plotK3pSelectionForRadius(const string &label, const string &title, int radius, int v=1){
    forEachExposure(label,radius,[&](int e){ plotK3pEqSelection(label,title,e,radius,v); });
}

// no selection: every stored window (windows with N = 0 are never written
// by MaskScanManye.C). Same layout as the k2 selections, plus k2 and q2
// since nothing fixes them here. q is shown on a log axis: it spans several
// decades (nearly fully masked windows reach q = 1) and most windows have
// q = 0, which is counted separately.
void plotNoCutDistributions(const string &label, const string &title, int exposure, int radius){

    string dir = label+"/e"+to_string(exposure)+"/r"+to_string(radius);
    string tag = label+"_e"+to_string(exposure)+"_r"+to_string(radius);

    cout << "[" << tag << "] no-cut distributions from '" << dir << "/' ..." << endl;

    TChain *ch = buildChain(dir,label);
    if(!ch->GetListOfFiles() || ch->GetListOfFiles()->GetEntries()==0){
        cout << "  no windowStatisticsManye_" << label << "_*.root files found -- skipping" << endl;
        delete ch;
        return;
    }

    int N, k[4];
    double r[4];
    ch->SetBranchAddress("N",&N);
    for(int c=0;c<4;c++){
        ch->SetBranchAddress(Form("k%d",c+1),&k[c]);
        ch->SetBranchAddress(Form("r%d",c+1),&r[c]);
    }

    Long64_t n = ch->GetEntries();
    vector<int> vN(n), vK[4];
    vector<double> vQ[4];
    for(int c=0;c<4;c++){ vK[c].resize(n); vQ[c].resize(n); }

    int maxN=0, maxK[4]={0,0,0,0};
    double qLo[4]={0,0,0,0}, qHi[4]={0,0,0,0};
    Long64_t sumN=0, sumK[4]={0,0,0,0}, nZero[4]={0,0,0,0}, nWith[4]={0,0,0,0};

    for(Long64_t i=0;i<n;i++){
        ch->GetEntry(i);
        vN[i]=N; sumN+=N; maxN=std::max(maxN,N);
        for(int c=0;c<4;c++){
            vK[c][i]=k[c]; vQ[c][i]=r[c];
            sumK[c]+=k[c];
            maxK[c]=std::max(maxK[c],k[c]);
            if(k[c]>0) nWith[c]++;
            if(r[c]<=0){ nZero[c]++; continue; }
            if(qLo[c]==0 || r[c]<qLo[c]) qLo[c]=r[c];
            qHi[c]=std::max(qHi[c],r[c]);
        }
    }
    delete ch;

    cout << "  " << n << " windows" << endl;
    if(n==0) return;

    applyPlotStyle();

    TH1D *hN = new TH1D(("hN_nocut_"+tag).c_str(),"",60,0,maxN+1);
    TH1D *hK[4], *hQ[4];
    for(int c=0;c<4;c++){
        hK[c] = new TH1D(Form("hK%d_nocut_%s",c+1,tag.c_str()),"",maxK[c]+1,-0.5,maxK[c]+0.5);
        if(qHi[c]<=0){ qLo[c]=1e-6; qHi[c]=1e-5; }
        if(qLo[c]>=qHi[c]) qLo[c]=qHi[c]/10;
        const int nb=40;
        vector<double> e(nb+1);
        double l0=std::log10(qLo[c]*0.95), l1=std::log10(qHi[c]*1.05);
        for(int b=0;b<=nb;b++) e[b]=std::pow(10.0,l0+(l1-l0)*b/nb);
        hQ[c] = new TH1D(Form("hQ%d_nocut_%s",c+1,tag.c_str()),"",nb,e.data());
    }
    for(Long64_t i=0;i<n;i++){
        hN->Fill(vN[i]);
        for(int c=0;c<4;c++){
            hK[c]->Fill(vK[c][i]);
            if(vQ[c][i]>0) hQ[c]->Fill(vQ[c][i]);
        }
    }

    const int CW = 1900, CH = 780;
    TCanvas *cv = new TCanvas(("cnocut_"+tag).c_str(),"",CW,CH);
    cv->SetBatch(kTRUE);

    string context = title+"  #upoint  exposure "+to_string(exposure)+"  #upoint  r = "+to_string(radius);
    double used = canvasHeader(cv,CH,"All windows (no cut)",context,
        Form("%lld windows     #LTN#GT = %.0f live pixels     fully masked windows (N = 0) are not stored",
             (Long64_t)n,(double)sumN/n));

    vector<TPad*> pads = gridPads(5,2,1-used-0.01,"pnc_"+tag+"_");

    const char *kT[4]  = {"k_{1}","k_{2}","k_{3}","k_{4+}"};
    const char *kH[4]  = {"1e hits per window","2e hits per window","3e hits per window","#geq4e hits per window"};
    const char *qT[4]  = {"q_{1} = k_{1}/N","q_{2} = k_{2}/N","q_{3} = k_{3}/N","q_{4+} = k_{4+}/N"};
    const char *qH[4]  = {"1e rate","2e rate","3e rate","#geq4e rate"};

    pads[0]->cd(); drawDistribution(hN,nColor(),nFill(),"N"); padHeader("Live pixels per window","");
    for(int c=0;c<4;c++){
        pads[1+c]->cd();
        drawDistribution(hK[c],chanColor(c),chanFill(c),kT[c]);
        padHeader(kH[c],"");

        pads[6+c]->cd();
        gPad->SetLogx();
        drawDistribution(hQ[c],chanColor(c),chanFill(c),qT[c]);
        padHeader(qH[c],Form("+ %lld at q = 0",nZero[c]));
    }

    pads[5]->cd();
    {
        TLatex t;
        t.SetNDC();
        t.SetTextFont(43);
        t.SetTextAlign(11);
        t.SetTextSize(18);
        t.SetTextColor(gStInk);
        t.DrawLatex(0.06,0.82,"Totals over all windows");
        t.SetTextSize(16);
        t.SetTextColor(gStInk2);
        t.DrawLatex(0.06,0.68,"channel");
        t.DrawLatex(0.34,0.68,"total k");
        t.DrawLatex(0.56,0.68,"#sum k / #sum N");
        t.DrawLatex(0.06,0.10,"windows are overlapping circles, so a hit");
        t.DrawLatex(0.06,0.03,"is counted in every window that contains it");
        for(int c=0;c<4;c++){
            double y = 0.56-0.105*c;
            TLine *sw = new TLine(0.06,y+0.012,0.12,y+0.012);
            sw->SetNDC(); sw->SetLineColor(chanColor(c)); sw->SetLineWidth(3); sw->Draw();
            t.SetTextColor(gStInk);
            t.DrawLatex(0.15,y,kT[c]);
            t.DrawLatex(0.34,y,Form("%lld",sumK[c]));
            t.DrawLatex(0.56,y,Form("%.2e",sumN>0 ? (double)sumK[c]/sumN : 0.0));
        }
    }

    string out = dir+"/nocut_distributions_"+tag+".pdf";
    cv->SaveAs(out.c_str());
    delete cv;
    tightPdf(out);

    cout << "  wrote " << out << endl;
}

void plotNoCutForRadius(const string &label, const string &title, int radius){

    for(int exposure : discoverTaggedSubdirs(label,'e')){
        vector<int> radii = discoverTaggedSubdirs(label+"/e"+to_string(exposure),'r');
        if(std::find(radii.begin(),radii.end(),radius)==radii.end()) continue;
        plotNoCutDistributions(label,title,exposure,radius);
    }
}

int PlotK2EqualsTwo(){

    gROOT->SetBatch(kTRUE);

    for(auto cat : vector<pair<string,string>>{{"goodquads","good quads"},{"notbadquads","not-bad quads"}}){
        plotNoCutForRadius(cat.first,cat.second,30);
        plotK2SelectionForRadius(cat.first,cat.second,30,1);
        plotK2SelectionForRadius(cat.first,cat.second,30,2);
        plotK3pSelectionForRadius(cat.first,cat.second,30,1);
    }

    return 0;
}
