#include <cmath>

#include "libpixie/reader.hh"

#include "Implant.hh"

namespace FDSi {
  ImplantEvent::ImplantEvent() {
    for (int cr=0; cr<FD_MAX_CRATES; ++cr) {
      for (int sl=0; sl<FD_MAX_SLOTS_PER_CRATE; ++sl) {
        for (int ch=0; ch<FD_MAX_CHANNELS_PER_BOARD; ++ch) {
          indexMap[cr][sl][ch] = -1;
          implantMap[cr][sl][ch] = NULL;
        }
      }
    }
    fit = new IonTrigger();
    rit = new IonTrigger();
    pin0 = new Pin();
    pin1 = new Pin();
    cross_scint = new Scint();
    img_scint = new Scint();
    ppac = new PPAC();
  }

  void ImplantEvent::ReadConf(std::string conffile)  {
    std::ifstream file(conffile.c_str());
    if (!file.is_open()) {
      std::cout << "Warning! " << conffile << " not open!" << std::endl;
    }
    std::string line; 

    int n = 0;
    while (std::getline(file, line)) {
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      std::stringstream ss(line);

      int crID, slID, chID;
      std::string name;
      int indx;
      ss >> crID >>  slID >>  chID >> name >> indx;

      indexMap[crID][slID][chID] = indx;
      if (name == "beta") {
        implantMap[crID][slID][chID] = (ImplantChannel**)&beta;
      }
      else if (name == "imp") { 
        implantMap[crID][slID][chID] = (ImplantChannel**)&imp;
      }
      else if (name == "fit") {
        implantMap[crID][slID][chID] = (ImplantChannel**)&fit;
      }
      else if (name == "rit") {
        implantMap[crID][slID][chID] = (ImplantChannel**)&rit;        
      }
      else if (name == "pin0") {
        implantMap[crID][slID][chID] = (ImplantChannel**)&pin0;
      }
      else if (name == "pin1") {
        implantMap[crID][slID][chID] = (ImplantChannel**)&pin1;
      }
      else if (name == "cross_scint") {
        implantMap[crID][slID][chID] = (ImplantChannel**)&cross_scint;
      }
      else if (name == "img_scint") {
        implantMap[crID][slID][chID] = (ImplantChannel**)&img_scint;
      }
      else if (name == "ppac") {
        implantMap[crID][slID][chID] = (ImplantChannel**)&ppac;
      }
    }
  }

  void ImplantEvent::SetBetaThresh(int indx, float val) {
    for (int i=0; i<nStored; ++i) {
      betas[i].thresh[indx] = val;
    }
  }

  void ImplantEvent::SetImplantThresh(int indx, float val) {
    for (int i=0; i<nStored; ++i) {
      imps[i].thresh[indx] = val;
    }
  }

  void ImplantEvent::Reset() {
    beta = &betas[betaCtr];
    imp = &imps[impCtr];
    beta->Reset();
    imp->Reset();
    fit->Reset();
    rit->Reset();
    pin0->Reset();
    pin1->Reset();
    cross_scint->Reset();
    img_scint->Reset();
  }

  Implant::Implant()
    : TOF(0), dE(0), time(0), present(false), valid(false), nHits(0), tracelen(0), trace(NULL), upperthresh(60000) {
      for (int i=0; i<5; ++i) {
        energy[0][i] = 0;
        energy[1][i] = 0;
        energy[2][i] = 0;
        fired[i] = 0;
      }
    }

  void Implant::SetMeas(PIXIE::Measurement &meas, int indx) {
    if (meas.finishCode) { return; }
    if (meas.outOfRange) { return; }
    if (meas.CFDForce) { return; }
    if (meas.eventEnergy < thresh[indx]) { return; }
    if (meas.eventEnergy > upperthresh) { return; }

    fired[indx] += 1;

    if (fired[indx] == 1 || meas.eventEnergy > energy[0][indx]) {
      if (indx == 0) {
        time = meas.eventTime;
      }

      energy[0][indx] = meas.eventEnergy;
      //todo: handle if traces are not enabled
      energy[1][indx] = meas.trace_meas[0].datum/16;
      energy[2][indx] = meas.trace_meas[5].datum;
      trace = &meas.trace[0];
      tracelen = meas.traceLength;
    }
  }

  int Implant::firedDynode() {
    bool valid = true;
    if (fired[0] < 1) { valid = false; }
    return valid;
  }

  int Implant::firedAnodes() {
    bool valid = true;
    for (int i=1; i<5; ++i) { 
      if (fired[i] < 1) { valid = false; }
    }
    return valid; 
  }

  int Implant::EnergySum(int type) {
    double sum = 0;
    for (int i=1; i<5; ++i) {
      sum += energy[type][i];
    }
    return sum;
  }

  int Implant::SetPos(int type, pspmtCal &cal) {
    if (!firedAnodes()) { return -1; }
    double sum = EnergySum(type);
    ypos[type] = (energy[type][1] + energy[type][2])/sum;
    xpos[type] = (energy[type][2] + energy[type][3])/sum;
    Calibrate(type, cal);
    return 0;
  }

  void Implant::Calibrate(int type, pspmtCal &cal) {
    std::pair<float,float> calpos = cal.cal(energy[type][0], xpos[type], ypos[type]);
    xcal[type] = calpos.first;
    ycal[type] = calpos.second;
  }

  void Implant::Reset() {
    valid = false;
    nHits = 0; 
    for (int i=0; i<5; ++i) {
      fired[i] = 0;
    }
  }

  Beta::Beta() : nCuts(0), time(0), present(false), valid(false), nHits(0), tdiff(-999), tracelen(0), upperthresh(60000) {
    for (int i=0; i<5; ++i) {
      energy[0][i] = 0;
      energy[1][i] = 0;
      energy[2][i] = 0;
      fired[i] = 0;
      trace[i] = NULL;
    }
  }

  void Beta::SetMeas(PIXIE::Measurement &meas, int indx) {
    if (meas.finishCode) { return; }
    if (meas.outOfRange) { return; }
    if (meas.CFDForce) { return; }
    if (meas.eventEnergy < thresh[indx]) { return; }
    if (meas.eventEnergy > upperthresh) { return; }

    fired[indx] += 1;

    if (fired[indx] == 1 || meas.eventEnergy > energy[0][indx]) {
      if (indx == 0) {
        time = meas.eventTime;
      }

      energy[0][indx] = meas.eventEnergy;
      //todo: handle if traces are not enabled
      energy[1][indx] = meas.trace_meas[0].datum/16;
      energy[2][indx] = meas.trace_meas[5].datum;

      trace[indx] = &meas.trace[0];
      tracelen = meas.traceLength;
    }
  }

  int Beta::firedDynode() {
    bool valid = true;
    if (fired[0] < 1) { valid = false; }
    return valid;
  }

  int Beta::firedAnodes() {
    bool valid = true;
    for (int i=1; i<5; ++i) { 
      if (fired[i] < 1) { valid = false; }
    }
    return valid; 
  }

  int Beta::EnergySum(int type) {
    double sum = 0;
    for (int i=1; i<5; ++i) {
      sum += energy[type][i];
    }
    return sum;
  }

  int Beta::SetPos(int type, pspmtCal &cal) {
    if (!firedAnodes()) { return -1; }
    double sum = EnergySum(type);
    ypos[type] = (energy[type][1] + energy[type][2])/sum;
    xpos[type] = (energy[type][2] + energy[type][3])/sum;
    Calibrate(type, cal);
    return 0;
  }

  void Beta::Calibrate(int type, pspmtCal &cal) {
    std::pair<float,float> calpos = cal.cal(energy[type][0], xpos[type], ypos[type]);
    xcal[type] = calpos.first;
    ycal[type] = calpos.second;
  }

  void Beta::Reset() {
    nHits = 0; 
    nCuts = 0;
    for (int i=0; i<5; ++i) {
      fired[i] = 0;
    }
    valid = false;
    promptGamma = false;
  }

  void IonTrigger::SetMeas(PIXIE::Measurement &meas, int indx) {
    if (meas.finishCode) { return; }
    if (meas.outOfRange) { return; }
    if (meas.CFDForce) { return; }
    if (fired == 0 || meas.eventEnergy > energy) {
      ++fired;
      time = meas.eventTime;
      energy = meas.eventEnergy;
    }
  }
  void IonTrigger::Reset() {
    fired = 0;
    energy = 0;
    valid = false;
  }

  Pin::Pin()
    : valid(false), ecal(0), gain(1), offset(0) {
      for (int i=0; i<2; ++i) {
        energy[i] = 0;
        time[i] = 0;
        fired[i] = 0;
      }
    }

  void Pin::SetMeas(PIXIE::Measurement &meas, int indx) {
    if (meas.finishCode) { return; }
    if (meas.outOfRange) { return; }
    if (meas.CFDForce) { return; }
    if (fired[indx] == 0 || meas.eventEnergy > energy[indx]) {
      ++fired[indx];
      time[indx] = meas.eventTime;
      energy[indx] = meas.eventEnergy;
    }
  }

  void Pin::Reset() {
    valid = false;
    fired[0] = 0;
    fired[1] = 0;
  }

  void Pin::SetCal(double g, double o) {
    gain = g;
    offset = o;
  }

  void Pin::Calibrate() {
    ecal = (energy[0] + PIXIE::Reader::Dither())*gain + offset;
  }

  Scint::Scint()
    : avtime(0), valid(false) {
      for (int i=0; i<2; ++i) {
        energy[i] = 0;
        time[i] = 0;
        fired[i] = 0;
        thresh[i] = 0;
      }
    }

  void Scint::SetMeas(PIXIE::Measurement &meas, int indx) {
    if (meas.finishCode) { return; }
    if (meas.outOfRange) { return; }
    if (meas.CFDForce) { return; }
    if (meas.eventEnergy <= thresh[indx]) { return; }
    if (fired[indx] == 0 || meas.eventEnergy > energy[indx]) {
      ++fired[indx];
      time[indx] = meas.eventTime;
      energy[indx] = meas.eventEnergy;
    }
  }

  void Scint::Reset() {
    valid = false;
    fired[0] = 0;
    fired[1] = 0;
  }

  PPAC::PPAC()
    : valid(false), avtime(0) {
      for (int i=0; i<4; ++i) {
        fired[i] = 0;
        time[i] = 0;
        energy[i] = 0;
      }
    }

  void PPAC::SetMeas(PIXIE::Measurement &meas, int indx) {
    if (meas.finishCode) { return; }
    if (meas.outOfRange) { return; }
    if (meas.CFDForce) { return; }
    if (fired[indx] == 0 || meas.eventEnergy > energy[indx]) {
      ++fired[indx];
      time[indx] = meas.eventTime;
      energy[indx] = meas.eventEnergy;
    }
  }

  void PPAC::validate()  { 
    valid = (fired[0] == 1 && fired[1] == 1 & fired[2] == 1 && fired[3] == 1);
    avtime = (time[0] + time[1] + time[2] + time[3])/4;
  }

  void PPAC::Reset() {
    valid = false;
    for (int i=0; i<4; ++i) {
      fired[i] = 0;
    }
  }

  void pspmtCal::read_calgraphs(std::string file) {
    FILE *cfile= fopen(file.c_str(), "ra");
    std::stringstream ss;
    char cline[2048];

    int ind = 0;
    while(std::fgets(cline, sizeof cline, cfile)!=NULL) {
      std::string line(cline);
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      ss.clear();
      ss.str(line);

      int id, npoints;
      ss >> id;
      ss >> npoints;
    
      lens[id] = npoints;
      ens[id] = new float[npoints];
      posxs[id] = new float[npoints];
      posys[id] = new float[npoints];

      int ct = 0;
      while (ct < npoints) {
        if (std::fgets(cline, sizeof cline, cfile) == NULL ) { return; }

        std::string line(cline);
        if (line.size() == 0) { continue; }
        if (line[0] == '#') { continue; }
        if (line[0] == ';') { continue; }

        ss.clear();
        ss.str(line);

        double e, x, y;
        ss >> e;
        ss >> x;
        ss >> y;

        ens[id][ct] = e;
        posxs[id][ct] = x;
        posys[id][ct] = y;

        ++ct;
      }

      ++ind;
    }  
    loaded = true;
  }

std::pair<float,float> pspmtCal::interp_en(float val, int i, int j) {
  float* en = &ens[i*48 + j][0];
  float* x = &posxs[i*48 + j][0];
  float* y = &posys[i*48 + j][0];
  int len = lens[i*48+j];

  if (val >= en[0]) { return {x[0], y[0]}; }
  if (val <= en[len-1]) { return {x[len-1], y[len-1]}; }

  for (int i=0; i<len-1; ++i) {
    if (val < en[i] && val >= en[i+1]) {
      float mx = (x[i+1]-x[i])/(en[i+1]-en[i]);
      float cx = x[i] - en[i]*mx;

      float my = (y[i+1]-y[i])/(en[i+1]-en[i]);
      float cy = y[i] - en[i]*my;
      return {mx*val + cx, my*val + cy};
    }
  }
  return {0,0};
}

  void pspmtCal::make_map() {
    for (int ien=0; ien<128; ++ien) {
      float energy = ien*512;
      for (int i=0; i<48; ++i) {
        for (int j=0; j<48; ++j) {
          std::pair<float,float> coord = interp_en(energy, i, j);
          cal_xposx[ien*48*48+i*48+j] = coord.first;
          cal_xposy[ien*48*48+i*48+j] = coord.second;

          cal_yposx[ien*48*48+j*48+i] = coord.first;
          cal_yposy[ien*48*48+j*48+i] = coord.second;
        }
      }
    }
  }

  int pspmtCal::find_indx(float x, float y, int ien, int axis) {

  float *x_cal;
  float *y_cal;

  int indx = -999;

  if (axis == 0) {
    x_cal = &cal_xposx[ien*48*48];
    y_cal = &cal_xposy[ien*48*48];
  }
  else if (axis == 1) {
    x_cal = &cal_yposy[ien*48*48];
    y_cal = &cal_yposx[ien*48*48];
  }
  else {
    std::cout << "Invalid axis passed to find_indx" << std::endl;
    return indx;
  }

  float low = interp(y, &y_cal[0], &x_cal[0]);
  
  if (x < low) {
    indx = -1;
    return indx;
  }
  
  float high = interp(y, &y_cal[47*48], &x_cal[47*48]);
  if (x >= high) {
    indx = 48;
    return indx;
  }
  for (int i = 0; i<47; ++i) {
    high = interp(y, &y_cal[(i+1)*48], &x_cal[(i+1)*48]);
    if (x >= low && x < high) {
      indx=i;
      break;
    }
    low = high;
  }
  return indx;
}

float pspmtCal::interp(float val, float* x, float* y) {
  int len=48;
  if (val <= x[0]) { return y[0]; }
  if (val >= x[len-1]) { return y[len-1]; }

  for (int i=0; i<len-1; ++i) {
    if (val > x[i] && val <= x[i+1]) {
      float m = (y[i+1]-y[i])/(x[i+1]-x[i]);
      float c = y[i] - x[i]*m;
      return m*val + c;
    }
  }
  return 0;
}

std::pair<float,float> pspmtCal::interp2(float x, float y, float x1, float y1, float x2, float y2, float x3, float y3, float x4, float y4) {
  float a1 = x1;
  float a2 = x2 - x1;
  float a3 = x4 - x1;
  float a4 = x1 - x2 + x3 - x4;

  float b1 = y1;
  float b2 = y2 - y1;
  float b3 = y4 - y1;
  float b4 = y1 - y2 + y3 - y4;

  float m,l;
  if (abs(a3) > abs(a2)) {
    float maa = a4*b2 - a2*b4;
    float mbb = a4*b1 - a1*b4 + a2*b3 - a3*b2 + x*b4 - y*a4;
    float mcc = a3*b1 - a1*b3 + x*b3 - y*a3;
  
    float ml = 0;
    if (maa != 0) {
      float mdet = std::sqrt(mbb*mbb - 4*maa*mcc);
      ml = (-mbb+mdet)/(2*maa);
    }
    else {
      ml = -mcc/mbb;
    }
    m = (x-a1-a2*ml)/(a3+a4*ml);
    l = ml;
  }
  else {
    float laa = a4*b3 - a3*b4;
    float lbb = a4*b1 - a1*b4 + a2*b3 - a3*b2 + x*b4 - y*a4;
    float lcc = a2*b1 - a1*b2 + x*b2 - y*a2;

    float lm = 0;
    if (laa != 0) {
      float ldet = std::sqrt(lbb*lbb - 4*laa*lcc);
      lm = (-lbb+ldet)/(2*laa);
      l = (x-a1-a3*lm)/(a2+a4*lm);
    }
    else {
      lm = -lcc/lbb;
      l = (x-a1-a3*lm)/(a2+a4*lm);
    }
    m = lm;
  }
  
  return {l,m};
}
  
  std::pair<float, float> pspmtCal::cal(float en, float x, float y) {
    if (loaded == false) { return {x,y}; }
        
        int iz = (int)(en/512.0);

        int indx = find_indx(x,y,iz,0);
        int indy = find_indx(y,x,iz,1);

        float x1,x2,x3,x4,y1,y2,y3,y4;
        x1 = cal_xposx[iz*48*48+indx*48+indy];
        x2 = cal_xposx[iz*48*48+(indx+1)*48+indy];
        x3 = cal_xposx[iz*48*48+(indx+1)*48+(indy+1)];
        x4 = cal_xposx[iz*48*48+(indx)*48+(indy+1)];

        y1 = cal_xposy[iz*48*48+indx*48+indy];
        y2 = cal_xposy[iz*48*48+(indx+1)*48+indy];
        y3 = cal_xposy[iz*48*48+(indx+1)*48+(indy+1)];
        y4 = cal_xposy[iz*48*48+(indx)*48+(indy+1)];

        std::pair<float,float> newpos = interp2(x,y,x1,y1,x2,y2,x3,y3,x4,y4);

        float newx,newy;
        newx = 0.2 + ((float)indx + newpos.first)/47.0 * (0.8 - 0.2);
        newy = 0.2 + ((float)indy + newpos.second)/47.0 * (0.8 - 0.2);
        return {newx, newy};
}
}
