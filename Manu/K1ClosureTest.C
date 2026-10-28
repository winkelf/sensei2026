// Closure test of the k1 / (k2 >= 1) independence method.
//
// Only images with no real 2e hit are used, so k1 cannot be correlated with
// 2e hits there. In each such image, fake "2e" hits are placed uniformly on
// the live pixels (same mask bits as MaskScanManye.C) at the dataset's real
// 2e rate for that quad selection, and a window is "selected" if a fake hit's
// position falls inside its circle -- the same geometry as a real hit, so the
// overlap between neighbouring windows is reproduced. The selection is then
// compared to the independence expectation with exactly the code of the real
// test (independencePages in PlotK1Independence.C):
//     each window contributes p = 1 - exp(-N*lambda_fake) to its k1 bin,
//     lambda_fake = (fake hits in the image) / (live pixels in the image).
// If the method is right, real and expected agree here. Repeating with
// kClosureToys independent sets of fake hits also measures how well the
// chi2 p-values are calibrated when windows overlap.
//
// Run from a dataset folder (like the other macros); needs that folder's
// config.yaml (raw files), window trees and imagerates/.
//
// Output per (quad selection, exposure), in <quads>/e<exp>/r<radius>/:
//   k1_independence_k2fake_<tag>.pdf   all toys combined, same page as the real test
//   k1_threecases_k2fake_<tag>.pdf     all windows / fake selection / independence
//   k1_closure_calibration_<tag>.pdf   per-toy p-values and pulls

#include "PlotK1Independence.C"
#include "TRandom3.h"
#include <yaml-cpp/yaml.h>

const int kClosureToys     = 20;
const int kClosureMaskBits = 0x56FC;     // same as MaskScanManye.C
const int kCNX = 3200, kCNY = 520;

struct CWin{ int x, y, N, k1; };

struct CImg{
    Long64_t Npix = 0;
    bool inCat[2] = {false,false};
    vector<CWin> wins[2];        // this image's windows in goodquads / notbadquads
};

IndepSpec specClosure(){
    return { {1}, false, 0, "k_{2}^{fake}", "n_{2}^{fake}", "#lambda^{fake}", "k2fake", "k1_threecases_k2fake" };
}

unsigned int closureSeed(int toy, int cat, const ImgKey &k){
    return 1u + 1000003u*toy + 7919u*cat + 104729u*std::get<0>(k)
              + 131u*std::get<1>(k) + 17u*std::get<2>(k) + 3u*std::get<3>(k);
}

void calibrationPage(const vector<IndepResult> &toys, const IndepResult &all, double lambda,
                     int nImg, Long64_t nWinPerToy, const string &dir, const string &tag,
                     const string &title, int exposure, int radius){

    applyPlotStyle();

    TH1D *hp = new TH1D(("hpval_"+tag).c_str(),"",10,0,1);
    TH1D *hz = new TH1D(("hpull_"+tag).c_str(),"",16,-4,4);
    double sz=0, szz=0, schi=0;
    int nLow=0;
    for(auto &r : toys){
        double z = r.totV>0 ? (r.totO-r.totE)/std::sqrt(r.totV) : 0;
        hp->Fill(std::min(r.pval,0.9999));
        hz->Fill(std::max(-3.999,std::min(3.999,z)));
        sz += z; szz += z*z; schi += r.ndf>0 ? r.chi2/r.ndf : 0;
        if(r.pval<0.05) nLow++;
    }
    int n = toys.size();
    double zMean = sz/n, zRms = std::sqrt(std::max(0.0,szz/n-zMean*zMean));

    const int CW = 1500, CH = 520;
    TCanvas *c = new TCanvas(("ccal_"+tag).c_str(),"",CW,CH);
    c->SetBatch(kTRUE);
    string context = title+"  #upoint  exposure "+to_string(exposure)+"  #upoint  r = "+to_string(radius);
    double used = canvasHeader(c,CH,"Closure test: are the p-values and errors calibrated?",context,
        Form("%d toys of fake 2e hits on the %d images with no real 2e hit     each toy: %lld windows, about %.0f selected",
             n,nImg,nWinPerToy,all.totO/n));
    vector<TPad*> pads = gridPads(3,1,1-used-0.01,"pcal_"+tag+"_");

    pads[0]->cd();
    drawDistribution(hp,gStBlue,TColor::GetColorTransparent(gStBlue,0.3),"#chi^{2} p-value of one toy","Toys",false);
    hp->SetMinimum(0);
    hp->SetMaximum(std::max(hp->GetMaximum(),n/10.0)*1.6);
    TLine *u = new TLine(0,n/10.0,1,n/10.0);
    u->SetLineColor(gStInk2); u->SetLineStyle(2); u->SetLineWidth(2); u->Draw();
    padHeader("p-value per toy","dashed: uniform");

    pads[1]->cd();
    drawDistribution(hz,gStBlue,TColor::GetColorTransparent(gStBlue,0.3),"(real - expected) / #sigma, total selected windows","Toys",false);
    hz->SetMinimum(0);
    hz->SetMaximum(std::max(hz->GetMaximum(),n*0.5*0.4)*1.6);
    {
        vector<double> gx, gy;
        for(int i=0;i<=160;i++){ double x=-4+8.0*i/160; gx.push_back(x); gy.push_back(n*0.5*std::exp(-0.5*x*x)/std::sqrt(2*M_PI)); }
        TGraph *g = new TGraph(gx.size(),gx.data(),gy.data());
        g->SetLineColor(gStInk2); g->SetLineStyle(2); g->SetLineWidth(2); g->Draw("L");
    }
    padHeader("Pull of the total per toy","dashed: N(0,1)");

    pads[2]->cd();
    {
        TLatex t;
        t.SetNDC();
        t.SetTextFont(43);
        t.SetTextAlign(11);
        t.SetTextSize(18);
        t.SetTextColor(gStInk);
        t.DrawLatex(0.06,0.84,"Summary");
        t.SetTextSize(16);
        auto line = [&](double y, const string &a, const string &b){
            t.SetTextColor(gStInk2); t.DrawLatex(0.06,y,a.c_str());
            t.SetTextColor(gStInk);  t.DrawLatex(0.62,y,b.c_str());
        };
        line(0.72,"fake 2e rate #lambda",Form("%.2e / live pixel",lambda));
        line(0.63,"pull mean (expect 0)",Form("%.2f",zMean));
        line(0.54,"pull rms (expect 1)",Form("%.2f",zRms));
        line(0.45,"mean #chi^{2}/ndf per toy (expect 1)",Form("%.2f",schi/n));
        line(0.36,"toys with p < 0.05 (expect 5%)",Form("%d of %d",nLow,n));
        line(0.24,"all toys combined",Form("%.0f vs %.1f #pm %.1f",all.totO,all.totE,std::sqrt(all.totV)));
        line(0.15,"  #chi^{2}/ndf, p",Form("%.1f / %d,  %.2g",all.chi2,all.ndf,all.pval));
    }

    string out = dir+"/k1_closure_calibration_"+tag+".pdf";
    c->SaveAs(out.c_str());
    delete c;
    tightPdf(out);
    cout << "  wrote " << out << endl;
}

void closureForExposure(int exposure, int radius, int nToys){

    const char *cats[2]   = {"goodquads","notbadquads"};
    const char *titles[2] = {"good quads","not-bad quads"};

    cout << "=== closure test, exposure " << exposure << " ===" << endl;

    map<ImgKey,ImgRate> rates = loadImageRates(exposure);
    if(rates.empty()){
        cout << "  no imagerates/e" << exposure << "/ -- run ImageRates.C first" << endl;
        return;
    }

    // windows of the images with no 2e hit, and the real 2e rate of each
    // quad selection (over all its images)
    map<ImgKey,CImg> imgs;
    double lambdaRef[2] = {0,0};

    for(int c=0;c<2;c++){

        string dir = string(cats[c])+"/e"+to_string(exposure)+"/r"+to_string(radius);
        TChain *ch = buildChain(dir,cats[c]);
        if(!ch->GetListOfFiles() || ch->GetListOfFiles()->GetEntries()==0){ delete ch; continue; }

        int run, lta, ohdu, x, y, N, k1, k2;
        ch->SetBranchAddress("runID",&run);
        ch->SetBranchAddress("LTA",&lta);
        ch->SetBranchAddress("OHDU",&ohdu);
        ch->SetBranchAddress("x",&x);
        ch->SetBranchAddress("y",&y);
        ch->SetBranchAddress("N",&N);
        ch->SetBranchAddress("k1",&k1);
        ch->SetBranchAddress("k2",&k2);

        set<ImgKey> catImgs;
        Long64_t nBadK2 = 0;
        int curTree = -1, curIdx = -1;

        for(Long64_t i=0;i<ch->GetEntries();i++){
            ch->GetEntry(i);
            if(ch->GetTreeNumber()!=curTree){
                curTree = ch->GetTreeNumber();
                curIdx  = fileIndexOf(ch->GetFile()->GetName());
            }
            ImgKey key = std::make_tuple(curIdx,run,lta,ohdu);
            auto it = rates.find(key);
            if(it==rates.end()) continue;
            catImgs.insert(key);
            if(it->second.n[1]!=0) continue;
            if(k2!=0) nBadK2++;
            CImg &im = imgs[key];
            im.Npix = it->second.Npix;
            im.inCat[c] = true;
            im.wins[c].push_back({x,y,N,k1});
        }
        delete ch;

        double n2=0, np=0;
        for(auto &key : catImgs){ n2 += rates[key].n[1]; np += rates[key].Npix; }
        lambdaRef[c] = np>0 ? n2/np : 0;

        if(nBadK2>0) cout << "  WARNING: " << nBadK2 << " windows with k2 > 0 in images with n2 = 0" << endl;
    }

    if(imgs.empty()){ cout << "  no images without 2e hits" << endl; return; }

    // live-pixel masks of those images, rebuilt from the raw files
    YAML::Node cfg = YAML::LoadFile("config.yaml");
    vector<string> files = cfg["exposures"][to_string(exposure)].as<vector<string>>();

    map<int,vector<ImgKey>> byFile;
    for(auto &kv : imgs) byFile[std::get<0>(kv.first)].push_back(kv.first);

    map<ImgKey,vector<bool>> masks;
    for(auto &kv : byFile){
        int idx = kv.first;
        if(idx<1 || idx>(int)files.size()){ cout << "  WARNING: no input file for index " << idx << endl; continue; }
        TFile *f = TFile::Open(files[idx-1].c_str());
        if(!f || f->IsZombie()){ cout << "  WARNING: cannot open " << files[idx-1] << endl; continue; }
        TTree *cal = (TTree*)f->Get("calPixTree");
        int cx, cy, cm, cl, co, cr;
        cal->SetBranchStatus("*",0);
        for(const char *b : {"x","y","mask","LTANAME","ohdu","RUNID"}) cal->SetBranchStatus(b,1);
        cal->SetBranchAddress("x",&cx); cal->SetBranchAddress("y",&cy); cal->SetBranchAddress("mask",&cm);
        cal->SetBranchAddress("LTANAME",&cl); cal->SetBranchAddress("ohdu",&co); cal->SetBranchAddress("RUNID",&cr);
        for(auto &key : kv.second) masks[key].assign(kCNX*kCNY,false);
        for(Long64_t i=0;i<cal->GetEntries();i++){
            cal->GetEntry(i);
            if(cx<0 || cx>=kCNX || cy<0 || cy>=kCNY) continue;
            if((cm & kClosureMaskBits)==0) continue;
            auto it = masks.find(std::make_tuple(idx,cr,cl,co));
            if(it!=masks.end()) it->second[cx*kCNY+cy] = true;
        }
        f->Close();
        delete f;
    }

    // check: N recomputed from the rebuilt mask must equal the stored N
    {
        Long64_t nChk=0, nDiff=0;
        for(auto &kv : imgs){
            auto mit = masks.find(kv.first);
            if(mit==masks.end()) continue;
            const vector<bool> &m = mit->second;
            for(int c=0;c<2;c++) for(auto &w : kv.second.wins[c]){
                int live = 0;
                for(int dx=-radius;dx<=radius;dx++) for(int dy=-radius;dy<=radius;dy++){
                    if(dx*dx+dy*dy>radius*radius) continue;
                    if(!m[(w.x+dx)*kCNY+(w.y+dy)]) live++;
                }
                nChk++;
                if(live!=w.N) nDiff++;
            }
        }
        cout << "  mask check: " << nChk-nDiff << " of " << nChk << " windows have the stored N" << endl;
    }

    IndepSpec sp = specClosure();

    for(int c=0;c<2;c++){

        string dir = string(cats[c])+"/e"+to_string(exposure)+"/r"+to_string(radius);
        string tag = string(cats[c])+"_e"+to_string(exposure)+"_r"+to_string(radius);

        int nImg = 0;
        Long64_t nWin = 0;
        for(auto &kv : imgs) if(kv.second.inCat[c] && masks.count(kv.first)){ nImg++; nWin += kv.second.wins[c].size(); }
        if(nImg==0) continue;

        double lam = lambdaRef[c];
        cout << "[" << tag << "] " << nImg << " images with no 2e hit, " << nWin
             << " windows, fake rate " << lam << " per live pixel, " << nToys << " toys" << endl;

        vector<int> aX; vector<double> aP; vector<char> aS;
        vector<IndepResult> toyRes;

        for(int t=0;t<nToys;t++){

            vector<int> vX; vector<double> vP; vector<char> vS;
            vX.reserve(nWin); vP.reserve(nWin); vS.reserve(nWin);

            for(auto &kv : imgs){
                const CImg &im = kv.second;
                if(!im.inCat[c] || im.Npix<=0) continue;
                auto mit = masks.find(kv.first);
                if(mit==masks.end()) continue;
                const vector<bool> &m = mit->second;

                TRandom3 rng(closureSeed(t,c,kv.first));
                int nFake = rng.Poisson(lam*im.Npix);
                vector<double> hx(nFake), hy(nFake);
                for(int h=0;h<nFake;h++){
                    int px, py;
                    do{ px = rng.Integer(kCNX); py = rng.Integer(kCNY); } while(m[px*kCNY+py]);
                    hx[h] = px+rng.Uniform(-0.499,0.499);
                    hy[h] = py+rng.Uniform(-0.499,0.499);
                }

                double lt = (double)nFake/im.Npix;
                for(auto &w : im.wins[c]){
                    bool sel = false;
                    for(int h=0;h<nFake && !sel;h++){
                        double dx = hx[h]-w.x, dy = hy[h]-w.y;
                        if(dx*dx+dy*dy<=radius*radius) sel = true;
                    }
                    vX.push_back(w.k1);
                    vP.push_back(1.0-std::exp(-w.N*lt));
                    vS.push_back(sel);
                }
            }

            toyRes.push_back(independencePages(vX,vP,vS,sp,dir,tag,titles[c],exposure,radius,false));
            aX.insert(aX.end(),vX.begin(),vX.end());
            aP.insert(aP.end(),vP.begin(),vP.end());
            aS.insert(aS.end(),vS.begin(),vS.end());
        }

        string foot = Form("%d toys combined  #upoint  fake 2e hits at #lambda = %.2e per live pixel on the %d images with no real 2e hit  #upoint  windows overlap",
                           nToys,lam,nImg);
        IndepResult all = independencePages(aX,aP,aS,sp,dir,tag,titles[c],exposure,radius,true,foot);

        calibrationPage(toyRes,all,lam,nImg,nWin,dir,tag,titles[c],exposure,radius);
    }
}

int K1ClosureTest(){

    gROOT->SetBatch(kTRUE);

    const int radius = 30;
    set<int> exposures;
    for(const char *cat : {"goodquads","notbadquads"})
        for(int e : discoverTaggedSubdirs(cat,'e')){
            vector<int> radii = discoverTaggedSubdirs(string(cat)+"/e"+to_string(e),'r');
            if(std::find(radii.begin(),radii.end(),radius)!=radii.end()) exposures.insert(e);
        }

    for(int e : exposures) closureForExposure(e,radius,kClosureToys);

    return 0;
}
