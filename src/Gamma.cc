#include "Gamma.hh"
#include "DEGAiHit.hh"

namespace FDSi {
  Gamma::Gamma(const GammaHit &hit, int ind, int id, int xtlid) {
    Energy = hit.Energy;
    MaxEn = hit.Energy;
    HitInds[0] = ind;
    MaxEnInd = ind;
    CloverID = id;
    CrystalID = xtlid;
    timestamp = hit.timestamp;
    ctimestamp = hit.timestamp;
    nHits = 1;
  }
  void Gamma::Set(const GammaHit &hit, int ind, int id, int xtlid) {
    Energy = hit.Energy;
    MaxEn = hit.Energy;
    HitInds[0] = ind;
    MaxEnInd = ind;
    CloverID = id;
    CrystalID = xtlid;
    timestamp = hit.timestamp;
    ctimestamp = hit.timestamp;
    nHits = 1;
  }
}
