#ifndef FDSI_DEGAIHIT_HH
#define FDSI_DEGAIHIT_HH

#include "libpixie/measurement.hh"

namespace FDSi {
  class DEGAiConf;
  class LaBrConf;
  
  class GammaHit {
  public:
    int CloverID;
    int CrystalID;
    
    double RawEnergy;
    double Energy;

    unsigned long long int timestamp;

    int PileUp;
    int CFDForce;
    int OutOfRange;

    int measIndx;

    bool Valid;
    
    double Theta;
    double Phi;
    GammaHit() {}; 
    GammaHit(const PIXIE::Measurement &meas, int id, int cryst);
    void Set(const PIXIE::Measurement &meas, int id, int cryst); 
  };

  class DEGAiHit : public GammaHit {
  public:
    DEGAiHit() {};
    DEGAiHit(const PIXIE::Measurement &meas, int id, int cryst, const DEGAiConf &conf, int measInd);
    void Set(const PIXIE::Measurement &meas, int id, int cryst, const DEGAiConf &conf, int measInd);
    void Print() { 
      std::cout << "Clover ID: " << CloverID << std::endl;
      std::cout << "Crystal ID: " << CrystalID << std::endl;
      std::cout << "Raw Energy: " << RawEnergy << std::endl;
      std::cout << "Timestamp: " << timestamp << std::endl;
      std::cout << "Pileup: " << PileUp << std::endl;
      std::cout << "Out of Range: " << OutOfRange << std::endl;
    } 
  };

  class LaBrHit : public GammaHit {
  public:
    LaBrHit() {};
    LaBrHit(const PIXIE::Measurement &meas, int id, int cryst, const LaBrConf &conf);
    void Set(const PIXIE::Measurement &meas, int id, int cryst, const LaBrConf &conf);
    void Print() { 
      std::cout << "Port ID: " << CloverID << std::endl;
      std::cout << "Detector ID: " << CrystalID << std::endl;
      std::cout << "Raw Energy: " << RawEnergy << std::endl;
      std::cout << "Timestamp: " << timestamp << std::endl;
      std::cout << "Pileup: " << PileUp << std::endl;
      std::cout << "Out of Range: " << OutOfRange << std::endl;
    } 
  };
  
}

#endif
