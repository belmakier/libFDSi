#include <fstream>
#include <cmath>

#include "DEGAi.hh"

namespace FDSi {
  DEGAiConf::DEGAiConf(std::string na, std::string conffile) {
    conf_name = na; 
    for (int i=0; i<FD_MAX_CRATES; ++i) {
      for (int j=0; j<FD_MAX_SLOTS_PER_CRATE; ++j) {
        for (int k=0; k<FD_MAX_CHANNELS_PER_BOARD; ++k) {
          CloverIDMap[i][j][k] = -1;
          CrystalMap[i][j][k] = -1;
        }
      }
    }

    for (int i=0; i<MAX_CLOVERS+1; ++i) {
      present[i] = 0;
      pos[i] = 0; //default furthest in
    }
    
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

      ss >> CloverID[n];
      ss >> nSubDets[n];
      ss >> pos[CloverID[n]];
      for (int i=0; i<nSubDets[n]; ++i) {
        ss >> Cr[n][i] >> Sl[n][i] >> Ch[n][i];
        CloverIDMap[Cr[n][i]][Sl[n][i]][Ch[n][i]] = CloverID[n];
        CrystalMap[Cr[n][i]][Sl[n][i]][Ch[n][i]] = i;
      } 

      present[CloverID[n]] = 1;

      n += 1;
    }
    nClovers = n;

    for (int j=0; j<nClovers; ++j) {
      int id = CloverID[j];
      for (int i=0; i<4; ++i) {
        gain[id][i] = 1.0;
        offset[id][i] = 0.0;
        toff[id][i] = 0.0;        
      }
      for (int i=0; i<4; ++i) {
        for (int i2=0; i2<4; ++i2) {
          ct_cal[id][i][i2] = 0.0;
        }
      }
    }
  }

  int DEGAiConf::ReadCal(std::string calfile) {
    std::ifstream file(calfile.c_str());
    if (!file.is_open()) {
      std::cout << "Warning! " << calfile << " not open" << std::endl;
      return -1;
    }
    std::string line;

    int n = 0;
    while (std::getline(file, line)) {
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      std::stringstream ss(line);

      int id, i;
      ss >> id >> i;
      ss >> offset[id][i] >> gain[id][i] >> toff[id][i];
      ++n;
    }
    return n;
  }

  int DEGAiConf::ReadCTCal(std::string calfile) {
    std::ifstream file(calfile.c_str());
    if (!file.is_open()) {
      std::cout << "Warning! " << calfile << " not open" << std::endl;
      return -1;
    }
    file.precision(7);
    std::string line;

    ct_cal_loaded = true;
    int n = 0;
    while (std::getline(file, line)) {
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      std::stringstream ss(line);
      ss.precision(7);
      
      int id, i, i2;
      double gain;
      ss >> id >> i >> i2 >> gain;
 
      ct_cal[id][i][i2] = gain;

      ++n;
    }
    return n;
  }
  
  int DEGAiConf::ReadAngleMap(std::string mapfile) {
    std::ifstream file(mapfile.c_str());
    if (!file.is_open()) {
      std::cout << "Warning! " << mapfile << " not open" << std::endl;
      return -1;
    }
    std::string line;

    int n = 0;
    while (std::getline(file, line)) {
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      std::stringstream ss(line);

      int id, i;
      double phin, phim, phif, thetan, thetam, thetaf;
      ss >> id >> i >> thetan >> phin >> thetam >> phim >> thetaf >> phif;

      if (present[id] != 1) { continue; }
      if (pos[id] == 0) { //near
        theta[id][i] = thetan;
        phi[id][i] = phin;
        ++n;
      }
      else if (pos[id] == 1) { //mid
        theta[id][i] = thetam;
        phi[id][i] = phim;
        ++n;
      }
      else if (pos[id] == 2) { //far
        theta[id][i] = thetam;
        phi[id][i] = phim;
        ++n;
      }
      else { //wherever you are
        std::cout << "Error: Invalid Clover " << id << " position " << pos[id] << " - should be: " << std::endl;
        std::cout << "    0: nearest" << std::endl;
        std::cout << "    1: middle" << std::endl;
        std::cout << "    2: furthest" << std::endl;
      }
    }
    return n;
  }

  int DEGAiConf::ReadEfficiencyCal(std::string fn) {
    std::ifstream file(fn.c_str());
    if (!file.is_open()) {
      std::cout << "Warning! " << fn << " not open" << std::endl;
      return -1;
    }    
    std::string line;
    getline(file, line);
    getline(file, line);
    std::stringstream ss(line);
    float abseff, a, aerr, b, berr, c, cerr, d, derr, e, eerr, f, ferr, g, gerr;
    ss >> abseff >> a >> aerr >> b >> berr >> c >> cerr >> d >> derr >> e >> eerr >> f >> ferr >> g >> gerr;
    AbsEff = abseff;
    EfficiencyCal[0] = a;
    EfficiencyCal[1] = b;
    EfficiencyCal[2] = c;
    EfficiencyCal[3] = d;
    EfficiencyCal[4] = e;
    EfficiencyCal[5] = f;
    EfficiencyCal[6] = g;
    return 1;
  }

  int DEGAiConf::ReadEffIDCal(std::string fn) {
    std::ifstream file(fn.c_str());
    if (!file.is_open()) {
      std::cout << "Warning! " << fn << " not open" << std::endl;
      return -1;
    }
    std::string line;
    getline(file, line);
    int n = 0;
    while (std::getline(file, line)) {
      if (line.size() == 0) { continue; }
      if (line[0] == '#') { continue; }
      if (line[0] == ';') { continue; }

      ++n;
      std::stringstream ss(line);
      int ID; 
      float abseff, relnorm, a, aerr, b, berr, c, cerr, d, derr, e, eerr, f, ferr, g, gerr;
      ss >> ID >> abseff >> relnorm >> a >> aerr >> b >> berr >> c >> cerr >> d >> derr >> e >> eerr >> f >> ferr >> g >> gerr;
      AbsEff_i[ID] = abseff;
      EffCal_i[ID][0] = a;
      EffCal_i[ID][1] = b;
      EffCal_i[ID][2] = c;
      EffCal_i[ID][3] = d;
      EffCal_i[ID][4] = e;
      EffCal_i[ID][5] = f;
      EffCal_i[ID][6] = g;
    }
    return n;
  }

  float DEGAiConf::Efficiency(const double &e) const {
    double x = std::log(e/100.0);
    double y = std::log(e/1000.0);
    float eff = std::pow((EfficiencyCal[0] + EfficiencyCal[1]*x + EfficiencyCal[2]*x*x), -EfficiencyCal[6]) +
      std::pow((EfficiencyCal[3] + EfficiencyCal[4]*y + EfficiencyCal[5]*y*y), -EfficiencyCal[6]);
    eff = std::pow(eff, -1.0/EfficiencyCal[6]);
    eff = AbsEff*std::exp(eff);
    return eff;
  }

  float DEGAiConf::Efficiency(const double &e, const int &ID) const {
    double x = std::log(e/100.0);
    double y = std::log(e/1000.0);
    float eff = std::pow((EffCal_i[ID][0] + EffCal_i[ID][1]*x + EffCal_i[ID][2]*x*x), -EffCal_i[ID][6]) +
      std::pow((EffCal_i[ID][3] + EffCal_i[ID][4]*y + EffCal_i[ID][5]*y*y), -EffCal_i[ID][6]);
    eff = std::pow(eff, -1.0/EffCal_i[ID][6]);
    eff = AbsEff_i[ID]*std::exp(eff);
    return eff;
  }

  void DEGAiConf::Print() {
    std::cout << conf_name << " configuration: " << std::endl;
    for (int i=0; i<nClovers; ++i) {
      if (present[CloverID[i]] == 0) { continue; }
      printf("%i    %i   ", CloverID[i], pos[CloverID[i]]);
      for (int j=0; j<nSubDets[i]; ++j) {
        printf("%i.%i.%i   ", Cr[i][j],Sl[i][j],Ch[i][j]);
      }
      printf("\n");
    }
    std::cout << conf_name << " calibration: " << std::endl;
    for (int i=0; i<nClovers; ++i) {
      if (present[CloverID[i]] == 0) { continue; }
      for (int j=0; j<nSubDets[i]; ++j) {
        printf("%i  %i     %5.4f    %5.4f   %5.4f\n", CloverID[i], j, offset[CloverID[i]][j], gain[CloverID[i]][j], toff[CloverID[i]][j]);
      }
    }
    if (ct_cal_loaded) {
    std::cout << conf_name << " CT correction: " << std::endl;
    for (int i=0; i<nClovers; ++i) {
      if (present[CloverID[i]] == 0) { continue; }
      std::cout << "Clover " << CloverID[i] << std::endl;
      printf("( 1.0000   %6.5f   %6.5f   %6.5f )\n", ct_cal[CloverID[i]][0][1], ct_cal[CloverID[i]][0][2], ct_cal[CloverID[i]][0][3]);
      printf("( %6.5f   1.0000   %6.5f   %6.5f )\n", ct_cal[CloverID[i]][1][0], ct_cal[CloverID[i]][1][2], ct_cal[CloverID[i]][1][3]);
      printf("( %6.5f   %6.5f   1.0000   %6.5f )\n", ct_cal[CloverID[i]][2][0], ct_cal[CloverID[i]][2][1], ct_cal[CloverID[i]][2][3]);
      printf("( %6.5f   %6.5f   %6.5f  1.0000 )\n", ct_cal[CloverID[i]][3][0], ct_cal[CloverID[i]][3][1], ct_cal[CloverID[i]][3][2]);
    }
    }
    std::cout << conf_name << " angles: " << std::endl;
    for (int i=0; i<nClovers; ++i) {
      if (present[CloverID[i]] == 0) { continue; }
      for (int j=0; j<nSubDets[i]; ++j) {
        printf("%i  %i     %5.4f    %5.4f\n", CloverID[i], j, theta[CloverID[i]][j], phi[CloverID[i]][j]);
      }
    }
    
  }

  void DEGAi::AddHit(const PIXIE::Measurement &meas, int id, int cryst) {
    if (nHits >= MAX_HITS) {
      std::cout << "Error! More than " << MAX_HITS << " hits" << std::endl;
      return;
    }
    if (cryst <= 3) {
      hits[nHits].Set(meas, id, cryst, conf);
      ++nHits;
    }
  }

  void DEGAi::MakeGammas() {
    nGammas = 0;
    nValid = 0;
    nPileUp = 0;
    nNonPrompt = 0;
    Gamma *gam;

    int nAddbacks = 0;
    std::vector<int> addbacks(nHits);
    for (int i=0; i<nHits; ++i) {
      long long int hitit = hits[i].timestamp; 
      if (hits[i].PileUp == 1 ) { nPileUp += 1; continue; }
      if (hits[i].Valid == false) { continue; }

      ++nValid;

      for (int j=0; j<nHits; ++j) {
        long long int hitjt = hits[j].timestamp; 
        if (hits[j].Valid == false) { continue; }

        if (i==j) { continue; }
        if (addbacks[i] == 1 || addbacks[j] == 1) { continue; }
        double tdiff = (double)((long long)(hitit - hitjt))/3276.8;
        if (abs(tdiff) > conf.AddBackTDiff) { continue; }
          
          double costhet = std::sin(hits[i].Theta)*std::sin(hits[j].Theta)*std::cos(hits[i].Phi - hits[j].Phi) + std::cos(hits[i].Theta)*std::cos(hits[j].Theta);
          if (costhet >= std::cos(conf.ClusterAngle*3.14159265/180.0) && costhet <= 1.0) {
            ++nAddbacks;

            addbacks[i] = 1;
            addbacks[j] = 1;

            gammas[nGammas].Set(hits[i], i, hits[i].CloverID, hits[i].CrystalID);
            gam = &gammas[nGammas];
            gam->Energy += hits[j].Energy;                         // add hit j
            gam->HitInds[gam->nHits] = j;
            ++gam->nHits;
            if (hits[j].Energy > gam->MaxEn) {
              gam->MaxEn = hits[j].Energy;
              gam->MaxEnInd = j;
              gam->CloverID = hits[j].CloverID;
              gam->CrystalID = hits[j].CrystalID;
              gam->timestamp = hits[j].timestamp;
              gam->ctimestamp = hits[j].timestamp;
            }

            if (hits[i].CloverID == hits[j].CloverID) {
               gam->Energy += hits[j].Energy * conf.ct_cal[hits[i].CloverID][hits[i].CrystalID][hits[j].CrystalID];
               gam->Energy += hits[i].Energy * conf.ct_cal[hits[i].CloverID][hits[j].CrystalID][hits[i].CrystalID];
               /*if (hits[i].CloverID == 7) {
                 if (gam->Energy > 1345 && gam->Energy < 1360) { 
               std::cout << hits[i].CloverID << "   " << gam->Energy << std::endl;
               std::cout << hits[i].Energy << "  " << hits[j].Energy << std::endl;
               std::cout << hits[i].CrystalID << "   " << hits[j].CrystalID << std::endl;
               std::cout << hits[i].Theta*180.0/3.141592 << "   " << hits[j].Theta*180.0/3.141592 << std::endl;
               }
               }*/
               //std::cout << gam->Energy << std::endl;
            }

            ++nGammas;
          }
        }
      }
    

    for (int i=0; i<nHits; ++i) {
      if (hits[i].Valid == false) { continue; }

      if (addbacks[i] == 1) { continue; }
      
      gammas[nGammas].Set(hits[i], i, hits[i].CloverID, hits[i].CrystalID);
      ++nGammas;
    }
  }

}
