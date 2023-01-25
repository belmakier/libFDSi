#ifndef FDSI_DEGAI_HH
#define FDSI_DEGAI_HH

#include <string.h> 

#include "Definitions.hh"
#include "DEGAiHit.hh"
#include "Gamma.hh"

namespace FDSi {
  class DEGAiConf {
  public:
    int nClovers; //number of clovers
    int CloverID[MAX_CLOVERS * 4]; //don't know how big this should be: crystal ID's
    int nSubDets[MAX_CLOVERS * 4]; //don't know how big this should be: crystal ID's
    int Cr[MAX_CLOVERS][4]; //crate ID    - these are all indexed by 
    int Sl[MAX_CLOVERS][4]; //slot IDs    - 0-3: crystals
    int Ch[MAX_CLOVERS][4]; //channel IDs
    int pos[MAX_CLOVERS+1]; //these are indexed by ID not read-in order
    int present[MAX_CLOVERS+1]; //these are indexed by ID not read-in order

    int CloverIDMap[FD_MAX_CRATES][FD_MAX_SLOTS_PER_CRATE][FD_MAX_CHANNELS_PER_BOARD];
    int CrystalMap[FD_MAX_CRATES][FD_MAX_SLOTS_PER_CRATE][FD_MAX_CHANNELS_PER_BOARD];

    double offset[MAX_CLOVERS * 4][4];
    double gain[MAX_CLOVERS * 4][4];
    double toff[MAX_CLOVERS * 4][4];
    
    double ct_cal[MAX_CLOVERS][4][4];    //when addback occurs
    bool ct_cal_loaded = false;

    double theta[MAX_CLOVERS][5];
    double phi[MAX_CLOVERS][5];

    bool AddBack = true;
    double EnThresh = 30; //keV
    double RawEnThresh = 80; //ADC units
    double RawUpThresh = 60000; //ADC units
    double AddBackTDiff = 500; //in ns
    double ClusterAngle = 30.0; //in degrees
    
    std::string conf_name;

    DEGAiConf() : conf_name("uninitialized") {};
    DEGAiConf(std::string na, std::string conffile);

    int ReadCal(std::string calfile);
    int ReadCTCal(std::string calfile);
    int ReadAngleMap(std::string mapfile);

    float AbsEff;
    float EfficiencyCal[7];

    float AbsEff_i[MAX_CLOVERS+1];
    float EffCal_i[MAX_CLOVERS+1][7]; //per detector
    
    float Efficiency(const double &e) const;
    float Efficiency(const double &e, const int &ID) const;

    int ReadEfficiencyCal(std::string calfile);
    int ReadEffIDCal(std::string calfile);
      
    void Print();
  };

  class DEGAi {
  public:
    DEGAiConf conf;
    
    int nHits;
    int nValid;
    int nGammas;
    int nPileUp;
    int nNonPrompt;
    
    DEGAiHit hits[MAX_HITS];
    Gamma gammas[MAX_GAMMAS];

    DEGAi() { };
    DEGAi(std::string na, std::string conffile) : conf(na, conffile) {};

    int ReadAngleMap(std::string mapfile) { return conf.ReadAngleMap(mapfile); }

    void AddHit(const PIXIE::Measurement &meas, int id, int cryst);
    void MakeGammas();
    void PrintConf() { conf.Print(); }
  };
}

#endif
