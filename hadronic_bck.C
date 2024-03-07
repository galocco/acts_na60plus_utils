#include <cmath>
#include <iostream>
#include <string>
#include <bitset>
#include <TF1.h>
#include <TRandom.h>

TF1 *fdNdYPi, *fdNdYK, *fdNdYP, *fdNdPtPi, *fdNdPtK, *fdNdPtP;
Double_t fNChPi, fNChK, fNChP;

FILE *fpcsv;

void InitBgGenerationPart(double, double, double, double, double,
                          double, double, double, double, double,
                          double, double, double, double,
                          double, double,
                          double, double, double, double, double);
void GenBgEvent(double, double, double);

const Double_t kMassP = 0.938;
const Double_t kMassK = 0.4937;
const Double_t kMassPi = 0.1396;
const Double_t kMassMu = 0.1057;
const Double_t kMassE = 0.0005;

void hadronic_bck(int nev, double beamSigma = 0, bool fullTargetSystem = false)
{

  long int particleNumber = 4503599644147712;
  // (1) 40 GeV pi K PRC66 (2002)054902
  // (2) 40 GeV p PRC83 (2011) 014901
  double y0BG = 2.22;  // gaussian y mean - 40 GeV
  double sigyBG = 1.2; // .. sigma
  double TBG = 0.17;   // inv.slope of thermal pt distribution

  double y0BGPi = 0.666;       // (1) mean of the two gaussians displaced symmetrically wrt y=0
  double y0BGKplus = 0.694;    // (1)
  double y0BGKminus = 0.569;   // (1)
  double y0BGP = 0.907;        // (2) checked approx
  double sigyBGPi = 0.872;     // (1) sigma of the two gaussians displaced symmetrically wrt y=0
  double sigyBGKplus = 0.725;  // (1)
  double sigyBGKminus = 0.635; // (1)
  double sigyBGP = 0.798;      // (2) checked approx
  double NBGPi = 74.;          //(1) normalization paramters of the sum of two gaussian
  double NBGKplus = 16.2;      // (1)
  double NBGKminus = 6.03;     // (1)
  double NBGP = 37.5;          // (2) checked approx
  double Piratio = 0.91;
  double TBGpi = 0.169;  //(1) from a fit in 0.2<M_t-m<0.7 GeV
  double TBGK = 0.229;   // (1) average of K+ (232 MeV) and K- (226 MeV)
  double TBGP = 0.25;    // ?
  double yminBG = 0;     // min y to generate
  double ymaxBG = 6;     // max y to generate
  double ptminBG = 0.01; // min pt to generate
  double ptmaxBG = 3;    // max pt to generate

  InitBgGenerationPart(NBGPi, NBGKplus, NBGKminus, NBGP, Piratio, y0BG, y0BGPi, y0BGKplus, y0BGKminus, y0BGP, sigyBGPi, sigyBGKplus, sigyBGKminus, sigyBGP, yminBG, ymaxBG, TBGpi, TBGK, TBGP, ptminBG, ptmaxBG);

  printf("pion   multiplicity in %f<y<%f = %f\n", yminBG, ymaxBG, fNChPi);
  printf("kaon   multiplicity in %f<y<%f = %f\n", yminBG, ymaxBG, fNChK);
  printf("proton multiplicity in %f<y<%f = %f\n", yminBG, ymaxBG, fNChP);

  double vX = 0, vY = 0, vZ = 0;

  char csvname[30];

  gRandom->SetSeed(1801);
  for (int i = 0; i < nev; i++)
  {

    int nzeros = 9;
    int eventNumber = i;
    while (eventNumber != 0) {
        eventNumber /= 10;
        nzeros--;
    }

    std::string eventFile = "event" + std::string(nzeros, '0') + std::to_string(eventNumber) + "-particles.csv";
    sprintf(csvname, eventFile.c_str(), i);
    fpcsv = fopen(csvname, "w");

    // fprintf(fp,"Event no. %d\n",i);
    printf("\nEvent no. %d\n", i);
    if (beamSigma > 0)
    {
      vX = gRandom->Gaus(0, beamSigma);
      vY = gRandom->Gaus(0, beamSigma);
    }
    if (fullTargetSystem)
    {
      int nTargets = 5;
      double targetThickness = 1.5; // mm
      double rnd = gRandom->Rndm();
      double distBetweenTargets = 12; // mm
      for (int i = 5; i > 0; i++)
      {
        if (rnd > 4. / 5.)
        {
          vZ = -(gRandom->Rndm() * targetThickness + i * (distBetweenTargets+targetThickness));
          break;
        }
      }
    }
    GenBgEvent(vX, vY, vZ);

    fclose(fpcsv);
  }
}

void InitBgGenerationPart(double NPi, double NKplus, double NKminus, double NP, double Piratio,
                          double y0, double y0Pi, double y0Kplus, double y0Kminus, double y0P,
                          double sigyPi, double sigyKplus, double sigyKminus, double sigyP,
                          double ymin, double ymax,
                          double Tpi, double TK, double TP, double ptmin, double ptmax)
{
  // initialize bg generation routines
  fNChPi = NPi * (1 + Piratio) * sigyPi * TMath::Sqrt(TMath::Pi() / 2.) * (TMath::Erf((ymax - y0 - y0Pi) / sqrt(2.) / sigyPi) + TMath::Erf((ymax - y0 + y0Pi) / sqrt(2.) / sigyPi) - TMath::Erf((ymin - y0 - y0Pi) / sqrt(2.) / sigyPi) - TMath::Erf((ymin - y0 + y0Pi) / sqrt(2.) / sigyPi));
  fNChK = NKplus * sigyKplus * TMath::Sqrt(TMath::Pi() / 2.) * (TMath::Erf((ymax - y0 - y0Kplus) / sqrt(2.) / sigyKplus) + TMath::Erf((ymax - y0 + y0Kplus) / sqrt(2.) / sigyKplus) - TMath::Erf((ymin - y0 - y0Kplus) / sqrt(2.) / sigyKplus) - TMath::Erf((ymin - y0 + y0Kplus) / sqrt(2.) / sigyKplus)) +
          NKminus * sigyKminus * TMath::Sqrt(TMath::Pi() / 2.) * (TMath::Erf((ymax - y0 - y0Kminus) / sqrt(2.) / sigyKminus) + TMath::Erf((ymax - y0 + y0Kminus) / sqrt(2.) / sigyKminus) - TMath::Erf((ymin - y0 - y0Kminus) / sqrt(2.) / sigyKminus) - TMath::Erf((ymin - y0 + y0Kminus) / sqrt(2.) / sigyKminus));
  fNChP = NP * sigyP * TMath::Sqrt(TMath::Pi() / 2.) * (TMath::Erf((ymax - y0 - y0P) / sqrt(2.) / sigyP) + TMath::Erf((ymax - y0 + y0P) / sqrt(2.) / sigyP) - TMath::Erf((ymin - y0 - y0P) / sqrt(2.) / sigyP) - TMath::Erf((ymin - y0 + y0P) / sqrt(2.) / sigyP));
  //
  fdNdYPi = new TF1("dndy", "exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )", ymin, ymax);
  fdNdYPi->SetParameters(y0, y0Pi, sigyPi);
  fdNdYK = new TF1("dndy", "exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2) )", ymin, ymax);
  fdNdYK->SetParameters(y0, y0Kplus, sigyKplus);
  fdNdYP = new TF1("dndy", "37.45*(exp( -0.5*pow( (x-[0]-[1])/[2],2) )+ exp( -0.5*pow( (x-[0]+[1])/[2],2)) )", ymin, ymax);
  fdNdYP->SetParameters(y0, y0P, sigyP);

  fdNdPtPi = new TF1("dndptPi", "x*exp(-sqrt(x*x+1.949e-02)/[0])", ptmin, ptmax); // pion
  fdNdPtK = new TF1("dndptK", "x*exp(-sqrt(x*x+0.493*0.493)/[0])", ptmin, ptmax); // kaon
  fdNdPtP = new TF1("dndptP", "x*exp(-sqrt(x*x+0.938*0.938)/[0])", ptmin, ptmax); // proton
  fdNdPtPi->SetParameter(0, Tpi);
  fdNdPtK->SetParameter(0, TK);
  fdNdPtP->SetParameter(0, TP);
}

void GenBgEvent(double x, double y, double z)
{

  if (fNChPi < 0 && fNChK < 0 && fNChP < 0)
    return;
  // generate bg. events from simple thermalPt-gaussian Y parameterization
  if (!fdNdYPi || !fdNdYK || !fdNdYP || !fdNdPtPi || !fdNdPtK || !fdNdPtP)
  {
    printf("Background generation was not initialized\n");
    exit(1);
  }

  int ntrTot = 0;
  // pions
  double ntrPi = gRandom->Poisson(fNChPi);
  printf("fNChPi=%f ntrPi=%f\n", fNChPi, ntrPi);
  printf("vX=%f vY=%f vZ=%f \n", x, y, z);
  // 2|0|14|0|0
  int particleNumber = 1;
  fprintf(fpcsv, "particle_id,particle_type,process,vx,vy,vz,vt,px,py,pz,m,q\n");
  for (int itr = 0; itr < ntrPi; itr++)
  {
    double yrap = fdNdYPi->GetRandom();
    double pt = fdNdPtPi->GetRandom();
    double phi = gRandom->Rndm() * TMath::Pi() * 2;
    float charge = gRandom->Rndm() > 0.52 ? 1 : -1;
    int ptypecsv = 211;
    double masscsv = 0.13957039;
    double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + kMassPi * kMassPi) * TMath::SinH(yrap)};

    //               in pdg pr x  y  z t  px py pz m  
    /* from BARCODE
    constexpr Value vertexPrimary() const { return level(0); }
    /// Return the secondary vertex identifier.
    constexpr Value vertexSecondary() const { return level(1); }
    /// Return the particle identifier.
    constexpr Value particle() const { return level(2); }
    /// Return the generation identifier.
    constexpr Value generation() const { return level(3); }
    /// Return the sub-particle identifier.
    constexpr Value subParticle() const { return level(4); 
    */
    //auto  = 

    std::string binaryString = "1" + std::string(12, '0') + std::bitset<16>(particleNumber++).to_string() + std::string(24, '0');
    unsigned long long barcode = std::stoull(binaryString, nullptr, 2);
    fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, ptypecsv, 0, x, y, z, 0., pxyz[0], pxyz[1], pxyz[2], masscsv, charge);
  }
  // kaons
  double ntrK = gRandom->Poisson(fNChK);
  printf("fNChK=%f ntrK=%f\n", fNChK, ntrK);
  for (int itr = 0; itr < ntrK; itr++)
  {
    double yrap = fdNdYK->GetRandom();
    double pt = fdNdPtK->GetRandom();
    double phi = gRandom->Rndm() * TMath::Pi() * 2;
    float charge = gRandom->Rndm() > 0.3 ? 1 : -1;
    int ptypecsv = 321;
    double masscsv = 0.493677;
    double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + kMassK * kMassK) * TMath::SinH(yrap)};

    std::string binaryString = "1" + std::string(12, '0') + std::bitset<16>(particleNumber++).to_string() + std::string(24, '0');
    unsigned long long barcode = std::stoull(binaryString, nullptr, 2);
    fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, ptypecsv, 0, x, y, z, 0., pxyz[0], pxyz[1], pxyz[2], masscsv, charge);
  }
  // protons
  double ntrP = gRandom->Poisson(fNChP);
  printf("fNChP=%f ntrP=%f\n", fNChP, ntrP);
  for (int itr = 0; itr < ntrP; itr++)
  {
    double yrap = fdNdYP->GetRandom();
    double pt = fdNdPtP->GetRandom();
    double phi = gRandom->Rndm() * TMath::Pi() * 2;
    double charge = 1;
    int ptypecsv = 2212;
    double masscsv = 0.938272119;
    double pxyz[3] = {pt * TMath::Cos(phi), pt * TMath::Sin(phi), TMath::Sqrt(pt * pt + kMassP * kMassP) * TMath::SinH(yrap)};
    std::string binaryString = "1" + std::string(12, '0') + std::bitset<16>(particleNumber++).to_string() + std::string(24, '0');
    unsigned long long barcode = std::stoull(binaryString, nullptr, 2);
    fprintf(fpcsv, "%llu,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", barcode, ptypecsv, 0, x, y, z, 0., pxyz[0], pxyz[1], pxyz[2], masscsv, charge);
  }
}
