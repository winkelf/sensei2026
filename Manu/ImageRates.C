// Per-image totals needed for the k1/k2 independence test: for every
// detector image (runID, LTA, OHDU) in the not-bad quads (a superset of the
// good quads, so one pass serves both categories), stores
//   Npix   -- unmasked pixels over the whole image (same maskBits as the
//             window scan)
//   n1..n4 -- hits with (int)e == 1,2,3,>=4 whose barycenter lands on an
//             unmasked pixel (same veto as MaskScanManye.C)
// Radius-independent. Window sums can't be used for this because the
// windows overlap (grid step = radius).
//
// Output: imagerates/e<exposure>/imageRates_<n>.root, tree "images", one
// file per input file so an interrupted run resumes where it stopped.

#include "MaskScanManye.C"

int imageRates(const char *fileName, int nFile, const char *exposure,
                const QuadMap &quadMap, int maskBits=0x56FC){

    string outDir  = string("imagerates/e")+exposure;
    string outName = outDir+"/imageRates_"+to_string(nFile)+".root";

    if(!gSystem->AccessPathName(outName.c_str())){
        cout << "Already exists, skipping: " << outName << endl;
        return 0;
    }

    TFile *file = TFile::Open(fileName);
    if(!file || file->IsZombie()){
        cout << "Cannot open " << fileName << endl;
        return 1;
    }

    TTree *hitTree    = (TTree*)file->Get("hitSumm");
    TTree *calPixTree = (TTree*)file->Get("calPixTree");

    int runID, lta, ohdu;
    float e, xBary, yBary;

    hitTree->SetBranchStatus("*",0);
    for(const char *b : {"runID","LTANAME","ohdu","e","xBary","yBary"})
        hitTree->SetBranchStatus(b,1);

    hitTree->SetBranchAddress("runID",&runID);
    hitTree->SetBranchAddress("LTANAME",&lta);
    hitTree->SetBranchAddress("ohdu",&ohdu);
    hitTree->SetBranchAddress("e",&e);
    hitTree->SetBranchAddress("xBary",&xBary);
    hitTree->SetBranchAddress("yBary",&yBary);

    int cal_x, cal_y, cal_mask, cal_lta, cal_ohdu, cal_runID;

    calPixTree->SetBranchStatus("*",0);
    for(const char *b : {"x","y","mask","LTANAME","ohdu","RUNID"})
        calPixTree->SetBranchStatus(b,1);

    calPixTree->SetBranchAddress("x",&cal_x);
    calPixTree->SetBranchAddress("y",&cal_y);
    calPixTree->SetBranchAddress("mask",&cal_mask);
    calPixTree->SetBranchAddress("LTANAME",&cal_lta);
    calPixTree->SetBranchAddress("ohdu",&cal_ohdu);
    calPixTree->SetBranchAddress("RUNID",&cal_runID);

    // images, discovered from hitSumm exactly as in MaskScanManye.C
    set<DetKey> detectors;
    Long64_t nEvents = hitTree->GetEntries();
    for(Long64_t i=0;i<nEvents;i++){
        hitTree->GetEntry(i);
        if(!quadAllowed(quadMap,lta,ohdu)) continue;
        detectors.insert(make_tuple(runID,lta,ohdu));
    }

    map<DetKey,vector<bool>> maskedPixels;
    Long64_t nPixels = calPixTree->GetEntries();
    for(Long64_t i=0;i<nPixels;i++){
        calPixTree->GetEntry(i);
        if(cal_x<0 || cal_x>=NX || cal_y<0 || cal_y>=NY) continue;
        if((cal_mask & maskBits)==0) continue;
        DetKey key = make_tuple(cal_runID,cal_lta,cal_ohdu);
        if(!detectors.count(key)) continue;
        getMask(maskedPixels,key)[pixIdx(cal_x,cal_y)] = true;
    }

    struct Counts{ Long64_t Npix=0; Long64_t n[4]={0,0,0,0}; };
    map<DetKey,Counts> counts;

    for(auto &det : detectors){
        auto &mask = getMask(maskedPixels,det);
        Long64_t masked = 0;
        for(bool b : mask) if(b) masked++;
        counts[det].Npix = (Long64_t)NX*NY - masked;
    }

    for(Long64_t i=0;i<nEvents;i++){
        hitTree->GetEntry(i);
        DetKey key = make_tuple(runID,lta,ohdu);
        auto it = counts.find(key);
        if(it==counts.end()) continue;

        int ix = (int)round(xBary);
        int iy = (int)round(yBary);
        if(ix<0 || ix>=NX || iy<0 || iy>=NY) continue;
        if(getMask(maskedPixels,key)[pixIdx(ix,iy)]) continue;

        int ne = (int)e;
        if(ne<1) continue;
        it->second.n[std::min(ne,4)-1]++;
    }

    file->Close();
    delete file;

    gSystem->mkdir(outDir.c_str(),kTRUE);
    TFile out(outName.c_str(),"RECREATE");
    TTree tree("images","per-image unmasked pixels and hit counts");

    int oRun, oLTA, oOHDU;
    Long64_t oNpix, oN1, oN2, oN3, oN4;
    tree.Branch("runID",&oRun,"runID/I");
    tree.Branch("LTA",&oLTA,"LTA/I");
    tree.Branch("OHDU",&oOHDU,"OHDU/I");
    tree.Branch("Npix",&oNpix,"Npix/L");
    tree.Branch("n1",&oN1,"n1/L");
    tree.Branch("n2",&oN2,"n2/L");
    tree.Branch("n3",&oN3,"n3/L");
    tree.Branch("n4",&oN4,"n4/L");

    for(auto &kv : counts){
        oRun  = get<0>(kv.first);
        oLTA  = get<1>(kv.first);
        oOHDU = get<2>(kv.first);
        oNpix = kv.second.Npix;
        oN1 = kv.second.n[0]; oN2 = kv.second.n[1];
        oN3 = kv.second.n[2]; oN4 = kv.second.n[3];
        tree.Fill();
    }

    tree.Write();
    out.Close();

    cout << "Wrote " << outName << " (" << counts.size() << " images)" << endl;

    return 0;
}

int ImageRates(){

    YAML::Node configFile = YAML::LoadFile("config.yaml");
    QuadMap notBadQuads = loadQuadMap(configFile["notbad_quads"]);

    YAML::Node exposureNode = configFile["exposures"];

    for(auto it = exposureNode.begin(); it!=exposureNode.end(); ++it){

        string exposure = it->first.as<string>();
        vector<string> files = it->second.as<vector<string>>();

        cout << "=== exposure=" << exposure << " ===" << endl;

        int nFile = 0;
        for(auto &f : files) imageRates(f.c_str(),++nFile,exposure.c_str(),notBadQuads);
    }

    cout << "Done. Run PlotK1Independence.C next." << endl;

    return 0;
}
