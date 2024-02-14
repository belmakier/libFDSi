#include "DEGAiHit.hh"
#include "DEGAi.hh"
#include "LaBr.hh"

#include "libpixie/reader.hh"

namespace FDSi {
  GammaHit::GammaHit(const PIXIE::Measurement &meas, int id, int cryst) { Set(meas, id, cryst); }

  void GammaHit::Set(const PIXIE::Measurement &meas, int id, int cryst) {
    Valid = true;
    CloverID = id;
    CrystalID = cryst;
    RawEnergy = meas.eventEnergy;
    timestamp = meas.eventTime;
    PileUp = meas.finishCode;
    CFDForce = meas.CFDForce;
    OutOfRange = meas.outOfRange;

    if (PileUp == 1) {
      Valid = false;
    }
    if (OutOfRange == 1) {
      Valid = false;
    }
    if (CFDForce == 1) {
      Valid = false;
    }
  }

  DEGAiHit::DEGAiHit(const PIXIE::Measurement &meas, int id, int cryst, const DEGAiConf &conf, int measInd) { Set(meas, id, cryst, conf, measInd); }

  void DEGAiHit::Set(const PIXIE::Measurement &meas, int id, int cryst, const DEGAiConf &conf, int measInd) {
    GammaHit::Set(meas, id, cryst);
    timestamp = timestamp + conf.toff[id][cryst]*3276.8;

    Theta = conf.theta[CloverID][CrystalID]*3.1415926525/180.0;
    Phi = conf.phi[CloverID][CrystalID]*3.1415926525/180.0;

    Energy = conf.gain[id][cryst]*(RawEnergy + PIXIE::Reader::Dither()) + conf.offset[id][cryst];

    measIndx = measInd;

    if (Energy < conf.EnThresh) {
      Valid = false;
    }
    if (RawEnergy < conf.RawEnThresh || RawEnergy > conf.RawUpThresh) {
       Valid = false;
    }
  }

  LaBrHit::LaBrHit(const PIXIE::Measurement &meas, int id, int cryst, const LaBrConf &conf) { Set(meas, id, cryst, conf); }

  void LaBrHit::Set(const PIXIE::Measurement &meas, int id, int cryst, const LaBrConf &conf) {
    GammaHit::Set(meas, id, cryst);
    timestamp = timestamp + conf.toff[id][cryst]*3276.8;

    Theta = conf.theta[CloverID][CrystalID]*3.1415926525/180.0;
    Phi = conf.phi[CloverID][CrystalID]*3.1415926525/180.0;

    float dither_en = RawEnergy + PIXIE::Reader::Dither();
    Energy = conf.quad[id][cryst]*dither_en*dither_en + conf.gain[id][cryst]*(dither_en) + conf.offset[id][cryst];

    if (Energy < conf.EnThresh) {
      Valid = false;
    }
    if (RawEnergy < conf.RawEnThresh || RawEnergy > conf.RawUpThresh) {
       Valid = false;
    }
  }
}
