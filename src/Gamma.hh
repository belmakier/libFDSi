#ifndef FDSI_GAMMA_HH
#define FDSI_GAMMA_HH

#include <vector>

#include "DEGAiHit.hh"
#include "Definitions.hh"

namespace FDSi {
  class Gamma {
  public:
    double Energy;
    double MaxEn;
    int MaxEnInd;
    int CloverID;
    int CrystalID;

    unsigned long long int timestamp;
    unsigned long long int ctimestamp;
    
    int HitInds[MAX_HITS];
    int nHits;
    
    double Theta;
    double Phi;

  public:
    Gamma() {};
    Gamma(const GammaHit &hit, int ind, int id, int xtlid);
    void Set(const GammaHit &hit, int ind, int id, int xtlid);
  };
}

#endif
