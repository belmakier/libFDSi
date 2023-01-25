#ifndef FDSI_IMPLANT_HH
#define FDSI_IMPLANT_HH

#include <string>

#include "libpixie/measurement.hh"

#include "Definitions.hh"

namespace FDSi {
  class ImplantChannel {
    public:
      virtual void SetMeas(PIXIE::Measurement &meas, int indx) {}
  };      

  class pspmtCal {
    public:
      bool loaded = false;
      float *ens[2304];
      float *posxs[2304];
      float *posys[2304];
      
      float *cal_xposx;
      float *cal_xposy;
  
      float *cal_yposx;
      float *cal_yposy;

      int lens[2304];
    public:
      pspmtCal() : loaded(false) { cal_xposx=new float[128*48*48]; cal_xposy=new float[128*48*48]; cal_yposx=new float[128*48*48]; cal_yposy=new float[128*48*48];} ;
      void read_calgraphs(std::string file);
      void make_map();
      std::pair<float,float> interp_en(float val, int xind, int yind);
      int find_indx(float x, float y, int ien, int axis);

      float interp(float val, float *x, float *y);
      std::pair<float,float> interp2(float x, float y, float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4);
      std::pair<float,float> cal(float en, float x, float y); 
  };

  class Implant : ImplantChannel {
    public:
      double TOF;
      double dE;
      unsigned long long time;

      double energy[3][5];

      double xpos[3];
      double ypos[3];
      double xcal[3];
      double ycal[3];

      double thresh[5];
      int upperthresh;

      bool present = false;
      bool valid = false;
      int nHits = 0;
      int cutID = 0;

      int fired[5];

      uint32_t tracelen;
      uint16_t *trace; 
    public:
      Implant();
      void SetMeas(PIXIE::Measurement &meas, int indx);
      int firedAnodes();
      int firedDynode();
      int EnergySum(int type);
      int SetPos(int type, pspmtCal &cal);
      void Calibrate(int type, pspmtCal &cal);
      void Reset();
  };

  class Beta : ImplantChannel {
    public:
      unsigned long long time = 0;
      bool present = false;
      bool valid = false;
      int nHits = 0;
      double tdiff = -999; //correlation time
      int nCuts;
      int firstCut;
      int cutIDs[MAX_BETACUTS];
      double tdiffs[MAX_BETACUTS];

      double energy[3][5];
      int fired[5];

      double xpos[3];
      double ypos[3];
      double xcal[3];
      double ycal[3];

      double thresh[5];
      int upperthresh;

      uint16_t *trace[5];
      int tracelen;
    
      bool promptGamma;
      int gammaIndex;
      double gammaEnergy;

    public:
      Beta();
      void SetMeas(PIXIE::Measurement &meas, int indx);
      int firedDynode();
      int firedAnodes();
      int EnergySum(int type);
      void Reset();
      int SetPos(int type, pspmtCal &cal);
      void Calibrate(int type, pspmtCal &cal);
  };

  class IonTrigger : ImplantChannel { 
    public:
      double thresh = 0;
      double energy = 0;
      unsigned long long time = 0;
      bool valid = 0;
      int fired = 0;
    public:
      IonTrigger() : energy(0), time(0), valid(false), fired(0) {}
      void SetMeas(PIXIE::Measurement &meas, int indx);
      void Reset();
  };

  class Pin : ImplantChannel {
    public:
      double energy[2];
      unsigned long long time[2];
      int fired[2];
      bool valid = false;
      double ecal;

      double gain;
      double offset;

      Pin();
      void SetMeas(PIXIE::Measurement &meas, int indx);
      void Reset();
      void SetCal(double g, double o);
      void Calibrate();

  };

  class Scint : ImplantChannel {
    public:
      double energy[2];
      unsigned long long time[2];
      int fired[2];
      double thresh[2];

      unsigned long long avtime = 0;
      bool valid = false;
      Scint();
      void SetMeas(PIXIE::Measurement &meas, int indx);
      void Reset();
  };

  class PPAC  : ImplantChannel {
    public:
      unsigned long long time[4];
      double energy[4];
      int fired[4];

      unsigned long long avtime;
      bool valid = false;
      PPAC(); 
      void SetMeas(PIXIE::Measurement &meas, int indx);
      void Reset();
      void validate();
  };

  class ImplantEvent {
    public:
      //maps (chan,slot,crate)-> index
      int indexMap[FD_MAX_CRATES][FD_MAX_SLOTS_PER_CRATE][FD_MAX_CHANNELS_PER_BOARD];
      //maps (ID)->(object,index)
      ImplantChannel **implantMap[FD_MAX_CRATES][FD_MAX_SLOTS_PER_CRATE][FD_MAX_CHANNELS_PER_BOARD];

      const static int nStored{200};
      int impCtr = 0;
      int betaCtr = 0;
      Implant imps[nStored];
      Beta betas[nStored];

      Beta *beta = &betas[0];    //current beta, implant
      Implant *imp = &imps[0];

      IonTrigger *fit;     
      IonTrigger *rit;     
      Pin *pin0;
      Pin *pin1;
      Scint *cross_scint;
      Scint *img_scint;
      PPAC *ppac;

    public:
      ImplantEvent();
      void ReadConf(std::string conffile);
      void Set(PIXIE::Measurement &meas) {
        ImplantChannel **ch = implantMap[meas.crateID][meas.slotID][meas.channelNumber];
        int indx = indexMap[meas.crateID][meas.slotID][meas.channelNumber];
        if (ch != NULL && indx >= 0) {
          (*ch)->SetMeas(meas, indx);
        }
      }
      void SetBetaThresh(int indx, float val);
      void SetImplantThresh(int indx, float val);
      void Reset();
  };

}

#endif
