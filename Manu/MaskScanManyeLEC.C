// Same category x exposure x radius x file scan as MaskScanManye.C, but
// with the LEC ("Local Event Cut") pixel-mask bit (0x2000 = 8192) added on
// top of the standard bad-pixel mask (0x56FC). This gives a single
// reference (k1-k4, N) per (category, exposure, radius) to compare against
// the q1-cut scan done on the standard-mask output -- see PlotQ1CutScan.C.
//
// Output goes to <category>/lec/e<exposure>/r<radius>/, separate from the
// standard-mask output at <category>/e<exposure>/r<radius>/, so this can be
// run at any time without touching (or needing) the existing results.

#include "MaskScanManye.C"

int MaskScanManyeLEC(){

    const int maskBits = 0x56FC | 8192;
    const char *pathTag = "lec/";

    YAML::Node configFile = YAML::LoadFile("config.yaml");
    vector<int> radii = configFile["radii"].as<vector<int>>();

    map<string,vector<string>> exposureFiles;

    YAML::Node exposureNode = configFile["exposures"];

    for(auto it = exposureNode.begin(); it!=exposureNode.end(); ++it){

        string exposure = it->first.as<string>();
        exposureFiles[exposure] = it->second.as<vector<string>>();
    }

    QuadMap goodQuads   = loadQuadMap(configFile["good_quads"]);
    QuadMap notBadQuads = loadQuadMap(configFile["notbad_quads"]);

    struct Category{
        string label;
        QuadMap quads;
    };

    vector<Category> categories = {
        { "goodquads",   goodQuads   },
        { "notbadquads", notBadQuads }
    };

    for(auto &cat : categories){

        for(auto &exp : exposureFiles){

            const string &exposure          = exp.first;
            const vector<string> &inputFiles = exp.second;

            for(int radius : radii){

                cout << "=== [LEC] " << cat.label << ", exposure=" << exposure
                     << ", radius=" << radius << " ===" << endl;

                int nImages = 0;

                for (const auto& rf : inputFiles) {
                  nImages+=1;
                  int qs = windowStudy(rf.c_str(),nImages,radius,exposure.c_str(),
                                        cat.quads,cat.label.c_str(),maskBits,pathTag);
                }
            }
        }
    }

    cout << "Done. Run PlotQ1CutScan.C to produce the comparison plots." << endl;

    return 0;
}
