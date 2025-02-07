#include "TApplication.h"
#include "TROOT.h"
#include "sPhenixStyle.C"
#include "sPhenixStyle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TRandom2.h"
#include "TMath.h"
#include <utility>  // for std::pair
#include <cstdio>
#include <iostream>
#include "TGraph.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TColor.h>
#include "stdlib.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <TMinuit.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>
#include <algorithm>
#include "dlUtility.h"
#include <TDatime.h>
static const float radius_EM = 93.5;
static const float minz_EM = -130.23;
static const float maxz_EM = 130.23;

static const float radius_IH = 127.503;
static const float minz_IH = -170.299;
static const float maxz_IH = 170.299;

static const float radius_OH = 225.87;
static const float minz_OH = -301.683;
static const float maxz_OH = 301.683;

float get_emcal_mineta_zcorrected(float zvertex) {
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zertex);
  float minlimit = emcal_mineta;
  if (ihcal_mineta > minlimit) minlimit = ihcal_mineta;
  if (ohcal_mineta > minlimit) minlimit = ohcal_mineta;
  float maxlimit = emcal_maxeta;
  if (ihcal_maxeta < maxlimit) maxlimit = ihcal_maxeta;
  if (ohcal_maxeta < maxlimit) maxlimit = ohcal_maxeta;
  minlimit += jet_radius;
  maxlimit -= jet_radius;
  return jet_eta < minlimit || jet_eta > maxlimit;
}

void get_scaledowns(int runnumber, int scaledowns[])
{

  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","","");

  if (db)
    {
      printf("Server info: %s\n", db->ServerInfo());
    }
  else
    {
      printf("bad\n");
    }


  TSQLRow *row;
  TSQLResult *res;
  TString cmd = "";
  char sql[1000];

  //cout << runnumber << endl;
  for (int is = 0; is < 64; is++)
    {
      sprintf(sql, "select scaledown%02d from gl1_scaledown where runnumber = %d;", is, runnumber);
      //printf("%s \n" , sql);                                                                      
      res = db->Query(sql);

      int nrows = res->GetRowCount();

      int nfields = res->GetFieldCount();
      for (int i = 0; i < nrows; i++) {
        row = res->Next();
        for (int j = 0; j < nfields; j++) {
          scaledowns[is] = stoi(row->GetField(j));

          if(is == 10 || is==17 || is==18 || is==19) cout  << is << ":" << scaledowns[is] << " ";
        }
        delete row;
      }


      delete res;
    }
  //cout << endl;
  delete db;
}

int get_scaledown17(int runnumber)
{

  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","","");

  if (db)
    {
      printf("Server info: %s\n", db->ServerInfo());
    }
  else
    {
      printf("bad\n");
    }

  TSQLRow *row;
  TSQLResult *res;
  char sql[1000];

  sprintf(sql, "select scaledown17 from gl1_scaledown where runnumber = %d;", runnumber);
  res = db->Query(sql);
  row = res->Next();
  int sd17 = stoi(row->GetField(0));
      
  delete row;
  delete res;
  delete db;
  return sd17;
}

int get_scaledown10(int runnumber)
{

  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","","");

  if (db)
    {
      printf("Server info: %s\n", db->ServerInfo());
    }
  else
    {
      printf("bad\n");
    }

  TSQLRow *row;
  TSQLResult *res;
  char sql[1000];

  sprintf(sql, "select scaledown10 from gl1_scaledown where runnumber = %d;", runnumber);
  res = db->Query(sql);
  row = res->Next();
  int sd10 = stoi(row->GetField(0));
      
  delete row;
  delete res;
  delete db;
  return sd10;
}

long unsigned int get_nmb(int runnumber)
{

  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","","");

  if (db)
    {
      printf("Server info: %s\n", db->ServerInfo());
    }
  else
    {
      printf("bad\n");
    }

  TSQLRow *row;
  TSQLResult *res;
  char sql[1000];

  sprintf(sql, "select scaled from gl1_scalers where runnumber = %d and index = 10;", runnumber);
  res = db->Query(sql);
  row = res->Next();
  long unsigned int nmb = stoi(row->GetField(0));
      
  delete row;
  delete res;
  delete db;
  return nmb;
}

int build_chi2hists(string filebase, int runnumber)
{
  TCanvas* c = new TCanvas("","",1000,1000);
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0);
  
  
  string inputfilename = filebase+".root";
  TFile *file = TFile::Open(inputfilename.c_str()); // Replace with your ROOT file name
    if (!file || file->IsZombie()) {
      //std:://cerr << "Error opening file!" << std::endl;
      return 1;
    }

    // Get the tree
    TTree *jet_tree = (TTree*)file->Get("jet_tree");
    if (!jet_tree) {
      //std:://cerr << "Error getting jet_tree!" << std::endl;
      return 1;
    }

    // Declare variables for branches
    Float_t ecc, theta, frcoh, frcem, eta, phi, jet_ET, dphi, subjet_ET, maxTowE, subTowE, maxTowDiff, maxETowChi2, zvtx;
    int isdijet, nBadChi2, maxETowChi2Det, maxETowIsZS;
    Float_t jetcompE[3][512], jetcompEta[3][512], jetcompPhi[3][512];
    Float_t maxTowChi2[3], l2pcEta;
    float dPhi2pc[1000], dEta2pc[1000];
    float jet_eta[10], jet_phi[10], jet_et[10], alljetfrcem[10], alljetfrcoh[10], emLayerJetPhi[10], emLayerJetEta[10], emLayerJetET[10], ohLayerJetPhi[10], ohLayerJetEta[10], ohLayerJetET[10], dPhiLayer[10];
    int jet_n, n2pc;//, nLayerEm, nLayerOh;
    long long unsigned int triggervec;
    unsigned int bbfqavec, elmbgvec;
    // Set branch addresses
    jet_tree->SetBranchAddress("ecc", &ecc);
    jet_tree->SetBranchAddress("theta", &theta);
    jet_tree->SetBranchAddress("frcoh", &frcoh);
    jet_tree->SetBranchAddress("frcem", &frcem);
    jet_tree->SetBranchAddress("eta", &eta);
    jet_tree->SetBranchAddress("phi", &phi);
    jet_tree->SetBranchAddress("jet_ET", &jet_ET);
    jet_tree->SetBranchAddress("dphi", &dphi);
    jet_tree->SetBranchAddress("subjet_ET", &subjet_ET);
    //jet_tree->SetBranchAddress("jetcompE", jetcompE);
    //jet_tree->SetBranchAddress("jetcompEta", jetcompEta);
    //jet_tree->SetBranchAddress("jetcompPhi", jetcompPhi);
    jet_tree->SetBranchAddress("isdijet",&isdijet);
    jet_tree->SetBranchAddress("maxTowChi2", maxTowChi2);
    jet_tree->SetBranchAddress("maxTowE", &maxTowE);
    jet_tree->SetBranchAddress("subTowE", &subTowE);
    jet_tree->SetBranchAddress("maxTowDiff", &maxTowDiff);
    jet_tree->SetBranchAddress("nBadChi2", &nBadChi2);
    jet_tree->SetBranchAddress("maxETowChi2", &maxETowChi2);
    jet_tree->SetBranchAddress("maxETowIsZS", &maxETowIsZS);
    jet_tree->SetBranchAddress("maxETowChi2Det",&maxETowChi2Det);
    jet_tree->SetBranchAddress("zvtx",&zvtx);
    jet_tree->SetBranchAddress("bbfqavec",&bbfqavec);
    jet_tree->SetBranchAddress("elmbgvec",&elmbgvec);
    jet_tree->SetBranchAddress("jet_et",jet_et);
    jet_tree->SetBranchAddress("jet_eta",jet_eta);
    jet_tree->SetBranchAddress("jet_phi",jet_phi);
    jet_tree->SetBranchAddress("jet_n",&jet_n);
    jet_tree->SetBranchAddress("alljetfrcem",alljetfrcem);
    jet_tree->SetBranchAddress("alljetfrcoh",alljetfrcoh);
    jet_tree->SetBranchAddress("n2pc",&n2pc);
    //jet_tree->SetBranchAddress("nLayerEm",&nLayerEm);
    //jet_tree->SetBranchAddress("nLayerOh",&nLayerOh);
    jet_tree->SetBranchAddress("dPhiLayer",dPhiLayer);
    jet_tree->SetBranchAddress("dPhi2pc",dPhi2pc);
    jet_tree->SetBranchAddress("dEta2pc",dEta2pc);
    jet_tree->SetBranchAddress("l2pcEta",&l2pcEta);
    /*
    jet_tree->SetBranchAddress("emLayerJetPhi",emLayerJetPhi);
    jet_tree->SetBranchAddress("emLayerJetEta",emLayerJetEta);
    jet_tree->SetBranchAddress("emLayerJetET",emLayerJetET);
    jet_tree->SetBranchAddress("ohLayerJetPhi",ohLayerJetPhi);
    jet_tree->SetBranchAddress("ohLayerJetEta",ohLayerJetEta);
    jet_tree->SetBranchAddress("ohLayerJetET",ohLayerJetET);
    */

    
    TH2F* asdich_all = new TH2F("asdich_all","",25,0,1,100,0,100);
    TH2F* asdich_fail = new TH2F("asdich_fail","",25,0,1,100,0,100);
    TH2F* asdich_dphi = new TH2F("asdich_dphi","",25,0,1,100,0,100);

    TH2F* h2_n2pc[2];
    TH2F* h2_dPhiLayer[2];
    TH2F* h2_n2pcMinus[2];
    TH2F* h2_dPhiLayerMinus[2];
    TH2F* h2_n2pcPlus[2];
    TH2F* h2_dPhiLayerPlus[2];
    for(int i=0; i<2; ++i)
      {
	h2_n2pc[i] = new TH2F(("dancheck_n2pc"+to_string(i)).c_str(),"",121,-1.2,1.2,257,-M_PI,M_PI);
	h2_n2pcMinus[i] = new TH2F(("dancheck_n2pcMinus"+to_string(i)).c_str(),"",121,-1.2,1.2,257,-M_PI,M_PI);
	h2_n2pcPlus[i] = new TH2F(("dancheck_n2pcPlus"+to_string(i)).c_str(),"",121,-1.2,1.2,257,-M_PI,M_PI);
	h2_dPhiLayer[i] = new TH2F(("dancheck_dPhiLayer"+to_string(i)).c_str(),"",80,-0.4,0.4,120,-0.1,1.1);
	h2_dPhiLayerMinus[i] = new TH2F(("dancheck_dPhiLayerMinus"+to_string(i)).c_str(),"",80,-0.4,0.4,120,-0.1,1.1);
	h2_dPhiLayerPlus[i] = new TH2F(("dancheck_dPhiLayerPlus"+to_string(i)).c_str(),"",80,-0.4,0.4,120,-0.1,1.1);
      }


    //std:://cerr << "set branches" << endl;
    // Create 2D histograms for all combinations of variables
    const int numHistograms = 36; // Number of histograms in each array
    const int numTypes = 92;
    // Arrays to hold the histograms
    //cout << "got branches" << endl;
    TH2F* h2_maxETowChi2_nBadChi2[numHistograms];
    TH2F* h2_maxETowChi2_maxTowDiff[numHistograms];
    TH2F* h2_maxETowChi2_chi2[numHistograms];
    TH2F* h2_maxETowChi2_frcoh[numHistograms];
    TH2F* h2_maxETowChi2_frcem[numHistograms];
    TH2F* h2_maxETowChi2_eta[numHistograms];
    TH2F* h2_maxETowChi2_phi[numHistograms];
    TH2F* h2_maxETowChi2_jet_ET[numHistograms];
    TH2F* h2_maxETowChi2_dphi[numHistograms];
    TH2F* h2_maxETowChi2_subjet_ET[numHistograms];
    TH2F* h2_maxETowChi2_ecc[numHistograms];
    TH2F* h2_maxETowChi2_maxTowE[numHistograms];
    TH2F* h2_maxETowChi2_subTowE[numHistograms];

    TH2F* h2_nBadChi2_maxTowDiff[numHistograms];
    TH2F* h2_nBadChi2_chi2[numHistograms];
    TH2F* h2_nBadChi2_frcoh[numHistograms];
    TH2F* h2_nBadChi2_frcem[numHistograms];
    TH2F* h2_nBadChi2_eta[numHistograms];
    TH2F* h2_nBadChi2_phi[numHistograms];
    TH2F* h2_nBadChi2_jet_ET[numHistograms];
    TH2F* h2_nBadChi2_dphi[numHistograms];
    TH2F* h2_nBadChi2_subjet_ET[numHistograms];
    TH2F* h2_nBadChi2_ecc[numHistograms];
    TH2F* h2_nBadChi2_maxTowE[numHistograms];
    TH2F* h2_nBadChi2_subTowE[numHistograms];

    TH2F* h2_maxTowDiff_chi2[numHistograms];
    TH2F* h2_maxTowDiff_frcoh[numHistograms];
    TH2F* h2_maxTowDiff_frcem[numHistograms];
    TH2F* h2_maxTowDiff_eta[numHistograms];
    TH2F* h2_maxTowDiff_phi[numHistograms];
    TH2F* h2_maxTowDiff_jet_ET[numHistograms];
    TH2F* h2_maxTowDiff_dphi[numHistograms];
    TH2F* h2_maxTowDiff_subjet_ET[numHistograms];
    TH2F* h2_maxTowDiff_ecc[numHistograms];
    TH2F* h2_maxTowDiff_maxTowE[numHistograms];
    TH2F* h2_maxTowDiff_subTowE[numHistograms];

    TH2F* h2_subTowE_chi2[numHistograms];
    TH2F* h2_subTowE_frcoh[numHistograms];
    TH2F* h2_subTowE_frcem[numHistograms];
    TH2F* h2_subTowE_eta[numHistograms];
    TH2F* h2_subTowE_phi[numHistograms];
    TH2F* h2_subTowE_jet_ET[numHistograms];
    TH2F* h2_subTowE_dphi[numHistograms];
    TH2F* h2_subTowE_subjet_ET[numHistograms];
    TH2F* h2_subTowE_ecc[numHistograms];
    TH2F* h2_subTowE_maxTowE[numHistograms];

    TH2F* h2_maxTowE_chi2[numHistograms];
    TH2F* h2_maxTowE_frcoh[numHistograms];
    TH2F* h2_maxTowE_frcem[numHistograms];
    TH2F* h2_maxTowE_eta[numHistograms];
    TH2F* h2_maxTowE_phi[numHistograms];
    TH2F* h2_maxTowE_jet_ET[numHistograms];
    TH2F* h2_maxTowE_dphi[numHistograms];
    TH2F* h2_maxTowE_subjet_ET[numHistograms];
    TH2F* h2_maxTowE_ecc[numHistograms];

    TH2F* h2_ecc_chi2[numHistograms];
    TH2F* h2_ecc_frcoh[numHistograms];
    TH2F* h2_ecc_frcem[numHistograms];
    TH2F* h2_ecc_eta[numHistograms];
    TH2F* h2_ecc_phi[numHistograms];
    TH2F* h2_ecc_jet_ET[numHistograms];
    TH2F* h2_ecc_dphi[numHistograms];
    TH2F* h2_ecc_subjet_ET[numHistograms];

    TH2F* h2_chi2_frcoh[numHistograms];
    TH2F* h2_chi2_frcem[numHistograms];
    TH2F* h2_chi2_eta[numHistograms];
    TH2F* h2_chi2_phi[numHistograms];
    TH2F* h2_chi2_jet_ET[numHistograms];
    TH2F* h2_chi2_dphi[numHistograms];
    TH2F* h2_chi2_subjet_ET[numHistograms];

    TH2F* h2_frcoh_frcem[numHistograms];
    TH2F* h2_frcoh_eta[numHistograms];
    TH2F* h2_frcoh_phi[numHistograms];
    TH2F* h2_frcoh_jet_ET[numHistograms];
    TH2F* h2_frcoh_dphi[numHistograms];
    TH2F* h2_frcoh_subjet_ET[numHistograms];

    TH2F* h2_frcem_eta[numHistograms];
    TH2F* h2_frcem_phi[numHistograms];
    TH2F* h2_frcem_jet_ET[numHistograms];
    TH2F* h2_frcem_dphi[numHistograms];
    TH2F* h2_frcem_subjet_ET[numHistograms];

    TH2F* h2_eta_phi[numHistograms];
    TH2F* h2_eta_jet_ET[numHistograms];
    TH2F* h2_eta_dphi[numHistograms];
    TH2F* h2_eta_subjet_ET[numHistograms];

    TH2F* h2_phi_jet_ET[numHistograms];
    TH2F* h2_phi_dphi[numHistograms];
    TH2F* h2_phi_subjet_ET[numHistograms];

    TH2F* h2_jet_ET_dphi[numHistograms];
    TH2F* h2_jet_ET_subjet_ET[numHistograms];

    TH2F* h2_subjet_ET_dphi[numHistograms];

    TH2F* h2_AJ_dphi[numHistograms];
    //cout << "made all hist pointers" << endl;
    TH2F** allhists[numTypes] = {
      h2_maxETowChi2_nBadChi2, h2_maxETowChi2_maxTowDiff, h2_maxETowChi2_chi2, h2_maxETowChi2_frcoh, h2_maxETowChi2_frcem, h2_maxETowChi2_eta, h2_maxETowChi2_phi, h2_maxETowChi2_jet_ET, h2_maxETowChi2_dphi, h2_maxETowChi2_subjet_ET, h2_maxETowChi2_ecc, h2_maxETowChi2_maxTowE, h2_maxETowChi2_subTowE,

      h2_nBadChi2_maxTowDiff, h2_nBadChi2_chi2, h2_nBadChi2_frcoh, h2_nBadChi2_frcem, h2_nBadChi2_eta, h2_nBadChi2_phi, h2_nBadChi2_jet_ET, h2_nBadChi2_dphi, h2_nBadChi2_subjet_ET, h2_nBadChi2_ecc, h2_nBadChi2_maxTowE, h2_nBadChi2_subTowE,

      h2_maxTowDiff_chi2, h2_maxTowDiff_frcoh, h2_maxTowDiff_frcem, h2_maxTowDiff_eta, h2_maxTowDiff_phi, h2_maxTowDiff_jet_ET, h2_maxTowDiff_dphi, h2_maxTowDiff_subjet_ET, h2_maxTowDiff_ecc, h2_maxTowDiff_maxTowE, h2_maxTowDiff_subTowE,
      h2_subTowE_chi2, h2_subTowE_frcoh, h2_subTowE_frcem, h2_subTowE_eta, h2_subTowE_phi, h2_subTowE_jet_ET, h2_subTowE_dphi, h2_subTowE_subjet_ET, h2_subTowE_ecc, h2_subTowE_maxTowE,
      h2_maxTowE_chi2, h2_maxTowE_frcoh, h2_maxTowE_frcem, h2_maxTowE_eta, h2_maxTowE_phi, h2_maxTowE_jet_ET, h2_maxTowE_dphi, h2_maxTowE_subjet_ET, h2_maxTowE_ecc,
      h2_ecc_chi2, h2_ecc_frcoh, h2_ecc_frcem, h2_ecc_eta, h2_ecc_phi, h2_ecc_jet_ET, h2_ecc_dphi, h2_ecc_subjet_ET,
      h2_chi2_frcoh, h2_chi2_frcem, h2_chi2_eta, h2_chi2_phi, h2_chi2_jet_ET, h2_chi2_dphi, h2_chi2_subjet_ET,
      h2_frcoh_frcem, h2_frcoh_eta, h2_frcoh_phi, h2_frcoh_jet_ET, h2_frcoh_dphi ,h2_frcoh_subjet_ET,
      h2_frcem_eta, h2_frcem_phi, h2_frcem_jet_ET, h2_frcem_dphi, h2_frcem_subjet_ET,
      h2_eta_phi, h2_eta_jet_ET, h2_eta_dphi, h2_eta_subjet_ET,
      h2_phi_jet_ET, h2_phi_dphi, h2_phi_subjet_ET,
      h2_jet_ET_dphi, h2_jet_ET_subjet_ET,
      h2_subjet_ET_dphi,
      h2_AJ_dphi
    };
    //cout << "put hists in array" << endl;
    // Create the histograms using a loop
    const char* names[] = {
      "maxETowChi2_nBadChi2", "maxETowChi2_maxTowDiff", "maxETowChi2_chi2", "maxETowChi2_frcoh", "maxETowChi2_frcem", "maxETowChi2_eta", "maxETowChi2_phi", "maxETowChi2_jet_ET", "maxETowChi2_dphi", "maxETowChi2_subjet_ET", "maxETowChi2_ecc", "maxETowChi2_maxTowE", "maxETowChi2_subTowE",
      "nBadChi2_maxTowDiff", "nBadChi2_chi2", "nBadChi2_frcoh", "nBadChi2_frcem", "nBadChi2_eta", "nBadChi2_phi", "nBadChi2_jet_ET", "nBadChi2_dphi", "nBadChi2_subjet_ET", "nBadChi2_ecc", "nBadChi2_maxTowE", "nBadChi2_subTowE",
      "maxTowDiff_chi2", "maxTowDiff_frcoh", "maxTowDiff_frcem", "maxTowDiff_eta", "maxTowDiff_phi", "maxTowDiff_jet_ET", "maxTowDiff_dphi", "maxTowDiff_subjet_ET", "maxTowDiff_ecc", "maxTowDiff_maxTowE", "maxTowDiff_subTowE",
      "subTowE_chi2", "subTowE_frcoh", "subTowE_frcem", "subTowE_eta", "subTowE_phi", "subTowE_jet_ET", "subTowE_dphi", "subTowE_subjet_ET", "subTowE_ecc", "subTowE_maxTowE",
      "maxTowE_chi2", "maxTowE_frcoh", "maxTowE_frcem", "maxTowE_eta", "maxTowE_phi", "maxTowE_jet_ET", "maxTowE_dphi", "maxTowE_subjet_ET", "maxTowE_ecc",
      "ecc_chi2", "ecc_frcoh", "ecc_frcem", "ecc_eta", "ecc_phi", "ecc_jet_ET", "ecc_dphi", "ecc_subjet_ET",
      "chi2_frcoh", "chi2_frcem", "chi2_eta", "chi2_phi", "chi2_jet_ET", "chi2_dphi", "chi2_subjet_ET",
      "frcoh_frcem", "frcoh_eta", "frcoh_phi", "frcoh_jet_ET", "frcoh_dphi", "frcoh_subjet_ET",
      "frcem_eta", "frcem_phi", "frcem_jet_ET", "frcem_dphi", "frcem_subjet_ET",
      "eta_phi", "eta_jet_ET", "eta_dphi", "eta_subjet_ET",
      "phi_jet_ET", "phi_dphi", "phi_subjet_ET",
      "jet_ET_dphi", "jet_ET_subjet_ET",
      "subjet_ET_dphi",
      "AJ_dphi"
    };
    
    const double xbins[] = {
      1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
      100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
      100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
      100, 100, 100, 100, 100, 100, 100, 100, 100,
      1, 1, 1, 1, 1, 1, 1, 1,
      1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6,
      1.25, 1.25, 1.25, 1.25, 1.25, 1.25,
      1.25, 1.25, 1.25, 1.25, 1.25,
      0.7, 0.7, 0.7, 0.7,
      3.14, 3.14, 3.14,
      100, 100,
      100,
      1
    };

    const double ybins[] = {
      64, 100, 1e6, 1.25, 1.25, 0.7, 3.14, 100, 3.15, 100, 1, 100, 100,
      100, 1e6, 1.25, 1.25, 0.7, 3.14, 100, 3.15, 100, 1, 100, 100,
      1e6, 1.25, 1.25, 0.7, 3.14, 100, 3.15, 100, 1, 100, 100,
      1e6, 1.25, 1.25, 0.7, 3.14, 100, 3.15, 100, 1, 100,
      1e6, 1.25, 1.25, 0.7, 3.14, 100, 3.15, 100, 1,
      1e6, 1.25, 1.25, 0.7, 3.14, 100, 3.15, 100,
      1.25, 1.25, 0.7, 3.14, 100, 3.15, 100,
      1.25, 0.7, 3.14, 100, 3.15, 100,
      0.7, 3.14, 100, 3.15, 100,
      3.14, 100, 3.15, 100,
      100, 3.15, 100,
      3.15, 100,
      3.15,
      3.15
    };
    //std:://cerr << "prep to make hists" << std::endl;
    //std:://cout << "got all arrays for making hists" << endl;
    // Loop to create histograms
    for (int i = 0; i < numTypes; ++i) {
      for (int j = 0; j < numHistograms; ++j) {
	std::string name = "h2_" + std::string(names[i]) + "_" + std::to_string(j);
	std::string title = ""; //std::string(titles[i]) + " " + std::to_string(j);

        // Set ranges based on index
        double xMin = xbins[i]==0.7?-0.7:0;
	if(xbins[i]==3.14) xMin = -3.14;
	if(xbins[i]==1.25) xMin = -0.25;
        double xMax = xbins[i];
        double yMin = ybins[i]==0.7?-0.7:0;
	//cerr << "got some bounds" << endl;
	if(ybins[i]==3.14) yMin = -3.14;
	if(ybins[i]==1.25) yMin = -0.25;
        double yMax = ybins[i];
	int nx = 100;
	int ny = 100;
	if(xbins[i] == 64) nx = 64;
        // Create the histogram
	//cerr << "create histo" <<i << " " << j << endl;
        allhists[i][j] = new TH2F(name.c_str(), title.c_str(), nx, xMin, xMax, ny, yMin, yMax);
	//cerr << "created" << endl;
      }
    }
    //std:://cerr << "made hists" << endl;
    //cout << "made hists" << endl;
    // Loop over entries in the tree
    const int nRatio = 7;
    TH1F* forRatio[nRatio];
    const int nbin = 9;
    float bins[nbin+1] = {12,16,21,26,32,38,45,52,60,70};
    //for(int i=-2; i<nbin-2; ++i) bins[i+2] = 15*pow(6,((float)i)/(nbin-4));
    for(int i=0; i<nRatio; ++i)
      {
	forRatio[i] = new TH1F(("h1_forRatio_"+to_string(i)).c_str(),"",i<8?1000:100,0,i<5?100:1);
      }
    TH1F* dijetCheckRat = new TH1F("dijetCheckRat","",1000,0,100);
    TH1F* dijetCheckRatFull = new TH1F("dijetCheckRatFull","",1000,0,100);
    const int nSpectra = 33;
    TH1F* jetSpectra[nSpectra];
    TH1F* xJ[4];
    xJ[0] = new TH1F("xJ0","",100,0,1);
    xJ[1] = new TH1F("xJ1","",100,0,1);
    xJ[2] = new TH1F("new_xJ2","",100,0,1);
    xJ[3] = new TH1F("new_xJ3","",100,0,1);
    
    for(int i=0; i<nSpectra; ++i)
      {
	jetSpectra[i] = new TH1F(("h1_jetSpectra_"+to_string(i)).c_str(),"",nbin,bins);
	//cout << jetSpectra[i] << endl;
      }
    bool cutArr[nSpectra];
    
    const int nzhist = 6;
    TH1F* zhists[nzhist];
    for(int i=0; i<nzhist; ++i)
      {
	zhists[i] = new TH1F(("data_zhist"+to_string(i)).c_str(),"",300,-150,150);
      }


    Long64_t nEntries = jet_tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        jet_tree->GetEntry(i);
	if(abs(zvtx) > 30) continue;
	if(check_bad_jet_eta(eta,zvtx,0.4)) continue;
	bool dPhiCut = (dphi < 3*M_PI/4 && isdijet); //(1-frcem-frcoh) > ((2.0/3.0)*frcoh);//((elmbgvec >> 4) & 1);
	if(bbfqavec != 0) cout << bbfqavec << endl;
	bool dhCut = ((bbfqavec >> 5) & 1);
	bool ihCut = (frcem+frcoh) < 0.65;
	bool loETCut = ((frcem < 0.1) && (jet_ET > (50*frcem+20))) && (dPhiCut || !isdijet);
	bool hiETCut = ((frcem > 0.9) && (jet_ET > (-50*frcem+75))) && (dPhiCut || !isdijet);
	bool chi2cut = jet_ET > 25 && maxETowChi2 < 10;
	bool specialLoETCut = (frcem < 0.1) && (jet_ET > (50*frcem+20));
	bool specialHiETCut = (frcem > 0.9) && (jet_ET > (-50*frcem+75));
	zhists[0]->Fill(zvtx);
	if(dhCut) zhists[1]->Fill(zvtx);
	if(ihCut) zhists[2]->Fill(zvtx);
	if(loETCut) zhists[3]->Fill(zvtx);
	if(hiETCut) zhists[4]->Fill(zvtx);
	cutArr[0]=false;
	cutArr[1]=dPhiCut;
	cutArr[2]=dhCut;
	cutArr[3]=ihCut;
	cutArr[4]=loETCut;
	cutArr[5]=hiETCut;
	cutArr[6]=(dPhiCut || dhCut);
	cutArr[7]=(dPhiCut || ihCut);
	cutArr[8]=(dPhiCut || loETCut);
	cutArr[9]=(dPhiCut || hiETCut);
	cutArr[10]=(dhCut || ihCut);
	cutArr[11]=(dhCut || loETCut);
	cutArr[12]=(dhCut || hiETCut);
	cutArr[13]=(ihCut || loETCut);
	cutArr[14]=(ihCut || hiETCut);
	cutArr[15]=(loETCut || hiETCut);
	cutArr[16]=(dPhiCut || dhCut || ihCut);
	cutArr[17]=(dPhiCut || dhCut || loETCut);
	cutArr[18]=(dPhiCut || dhCut || hiETCut);
	cutArr[19]=(dPhiCut || ihCut || loETCut);
	cutArr[20]=(dPhiCut || ihCut || hiETCut);
	cutArr[21]=(dPhiCut || loETCut || hiETCut);
	cutArr[22]=(dhCut || ihCut || loETCut);
	cutArr[23]=(dhCut || ihCut || hiETCut);
	cutArr[24]=(dhCut || loETCut || hiETCut);
	cutArr[25]=(ihCut || loETCut || hiETCut);
	cutArr[26]=(dPhiCut || dhCut || ihCut || loETCut);
	cutArr[27]=(dPhiCut || dhCut || ihCut || hiETCut);
	cutArr[28]=(dPhiCut || dhCut || loETCut || hiETCut);
	cutArr[29]=(dPhiCut || ihCut || loETCut || hiETCut);
	cutArr[30]=(dhCut || ihCut || specialLoETCut || specialHiETCut);
	cutArr[31]=(dhCut || ihCut || loETCut || hiETCut);// || chi2cut);
	cutArr[32]=chi2cut;
	if(isdijet && (specialHiETCut || specialLoETCut) && !cutArr[31])
	  {
	    asdich_fail->Fill((jet_ET-subjet_ET)/(jet_ET+subjet_ET),jet_ET);
	  }
	if(!cutArr[31] && isdijet)
	  {
	    asdich_all->Fill((jet_ET-subjet_ET)/(jet_ET+subjet_ET),jet_ET);
	    if(dphi > 3*M_PI/4) asdich_dphi->Fill((jet_ET-subjet_ET)/(jet_ET+subjet_ET),jet_ET);
	  }
	if(!cutArr[31]) zhists[5]->Fill(zvtx);
	float closejetdphi = M_PI;
	for(int j=0; j<jet_n; ++j)
	  {
	    float testdphi = phi - jet_phi[j];
	    if(testdphi < 0.05) continue;
	    if(check_bad_jet_eta(jet_eta[j],zvtx,0.4)) continue;
	    if(testdphi > M_PI) testdphi = 2*M_PI - testdphi;
	    if(testdphi < closejetdphi) closejetdphi = testdphi;
	  }
	for(int j=0; j<nSpectra; ++j)
	  {
	    if(!cutArr[j])
	      {
		for(int k=0; k<jet_n; ++k)
		  {
		    if(check_bad_jet_eta(jet_eta[k],zvtx,0.4)) continue;
		    if(jet_et[k] > 8)
		      {
			jetSpectra[j]->Fill(jet_et[k]);
			
		      }
		  }
	      }
	  }
	//cout << "test" << endl;
	for(int j=0; j<n2pc; ++j)
	  {
	    h2_n2pc[0]->Fill(dEta2pc[j],dPhi2pc[j]);
	    if(l2pcEta > 0) h2_n2pcPlus[0]->Fill(dEta2pc[j],dPhi2pc[j]);
	    else h2_n2pcMinus[0]->Fill(dEta2pc[j],dPhi2pc[j]);
	    if(!cutArr[31]) 
	      {
		h2_n2pc[1]->Fill(dEta2pc[j],dPhi2pc[j]);
		if(l2pcEta > 0) h2_n2pcPlus[1]->Fill(dEta2pc[j],dPhi2pc[j]);
		else h2_n2pcMinus[1]->Fill(dEta2pc[j],dPhi2pc[j]);
	      }
	  }
	
	for(int j=0; j<jet_n; ++j)
	  {
	    h2_dPhiLayer[0]->Fill(dPhiLayer[j],alljetfrcem[j]);
	    if(l2pcEta > 0) h2_dPhiLayerPlus[0]->Fill(dPhiLayer[j],alljetfrcem[j]);
	    else h2_dPhiLayerMinus[0]->Fill(dPhiLayer[j],alljetfrcem[j]);
	    if(!cutArr[31])
	      {
		h2_dPhiLayer[1]->Fill(dPhiLayer[j],alljetfrcem[j]);
		if(l2pcEta > 0) h2_dPhiLayerPlus[1]->Fill(dPhiLayer[j],alljetfrcem[j]);
		else h2_dPhiLayerMinus[1]->Fill(dPhiLayer[j],alljetfrcem[j]);
	      }
	  }
	//cout << "test2" << endl;
	forRatio[0]->Fill(jet_ET);
	if(isdijet && dphi > 3*M_PI/4) forRatio[1]->Fill(jet_ET);
	if(!cutArr[1] && isdijet) forRatio[2]->Fill(jet_ET);
	if(!cutArr[31]) forRatio[3]->Fill(jet_ET);
	if(!cutArr[31] && isdijet) forRatio[4]->Fill(jet_ET);
	if(!cutArr[30]) dijetCheckRatFull->Fill(jet_ET);
	if(!cutArr[30] && isdijet) dijetCheckRat->Fill(jet_ET);
	if(isdijet && jet_ET > 14 && subjet_ET > 10)
	  {
	    forRatio[5]->Fill((jet_ET-subjet_ET)/(jet_ET+subjet_ET));
	    if(jet_ET > 30 && subjet_ET > 20) xJ[2]->Fill(subjet_ET/jet_ET);
	    xJ[0]->Fill(subjet_ET/jet_ET);
	  }
	if(isdijet && !cutArr[31] && jet_ET > 14 && subjet_ET > 10)
	  {
	    forRatio[6]->Fill((jet_ET-subjet_ET)/(jet_ET+subjet_ET));
	    if(jet_ET > 30 && subjet_ET > 20) xJ[3]->Fill(subjet_ET/jet_ET);
	    xJ[1]->Fill(subjet_ET/jet_ET);
	  }
	//cout << "filling (got entry). isdijet = " << isdijet << endl;
	float chi2 = maxTowChi2[0];
	if(maxTowChi2[2] > chi2) chi2 = maxTowChi2[2];
	string det;
	if(maxETowChi2Det == 0) det = "EMCal";
	else if(maxETowChi2Det == 1) det = "IHCal";
	else if(maxETowChi2Det == 2) det = "OHCal";
	else det = "Error";
	//int whichhist = isdijet%2 + 2*maxETowIsZS;
	int whichhist[3];
	int threshes[numHistograms/6] = {8,15,20,25,35,40};
	int slt[numHistograms/6] = {8,8,10,10,15,15};
	for(int j=0; j<numHistograms/6; ++j)
	  {
	    if(jet_ET < threshes[j]) continue;
	    if(cutArr[31]) continue;
	    whichhist[2] = j+ 2*numHistograms/3;
	    if(isdijet)
	      {
		whichhist[0] = j + (dphi>3*M_PI/4?numHistograms/6:0);
		whichhist[1] = j+ numHistograms/3 + (closejetdphi < M_PI/4? numHistograms/6:0);

	      }
	    else
	      {
		whichhist[0] = -1;
		whichhist[1] = -1;
	      }
	    for(int k=0; k<3; ++k)
	      {
		if(whichhist[k] < 0 || whichhist[k] > numHistograms-1) continue;


		if(isdijet)
		  {
		    h2_AJ_dphi[whichhist[k]]->Fill((jet_ET-subjet_ET)/(jet_ET+subjet_ET),dphi);
		  }
	//if(maxTowE > 15 && maxETowChi2 < 100) //cout << "Detector: " << det << " tower ET: " << maxTowE << " chi2: " << maxETowChi2 << " is ZS: " << maxETowIsZS <<endl;
	//if(cutArr[31]) continue;
		h2_maxETowChi2_nBadChi2[whichhist[k]]->Fill(maxETowChi2, nBadChi2);
		h2_maxETowChi2_maxTowDiff[whichhist[k]]->Fill(maxETowChi2, maxTowDiff);
		h2_maxETowChi2_subTowE[whichhist[k]]->Fill(maxETowChi2, subTowE);
		h2_maxETowChi2_maxTowE[whichhist[k]]->Fill(maxETowChi2, maxTowE);
		h2_maxETowChi2_ecc[whichhist[k]]->Fill(maxETowChi2, ecc);
		h2_maxETowChi2_chi2[whichhist[k]]->Fill(maxETowChi2, chi2);
		h2_maxETowChi2_frcoh[whichhist[k]]->Fill(maxETowChi2, frcoh);
		h2_maxETowChi2_frcem[whichhist[k]]->Fill(maxETowChi2, frcem);
		h2_maxETowChi2_eta[whichhist[k]]->Fill(maxETowChi2, eta);
		h2_maxETowChi2_phi[whichhist[k]]->Fill(maxETowChi2, phi);
		h2_maxETowChi2_jet_ET[whichhist[k]]->Fill(maxETowChi2, jet_ET);
		h2_maxETowChi2_dphi[whichhist[k]]->Fill(maxETowChi2, dphi);
		h2_maxETowChi2_subjet_ET[whichhist[k]]->Fill(maxETowChi2, subjet_ET);
		
		h2_nBadChi2_maxTowDiff[whichhist[k]]->Fill(nBadChi2, maxTowDiff);
		h2_nBadChi2_subTowE[whichhist[k]]->Fill(nBadChi2, subTowE);
		h2_nBadChi2_maxTowE[whichhist[k]]->Fill(nBadChi2, maxTowE);
		h2_nBadChi2_ecc[whichhist[k]]->Fill(nBadChi2, ecc);
		h2_nBadChi2_chi2[whichhist[k]]->Fill(nBadChi2, chi2);
		h2_nBadChi2_frcoh[whichhist[k]]->Fill(nBadChi2, frcoh);
		h2_nBadChi2_frcem[whichhist[k]]->Fill(nBadChi2, frcem);
		h2_nBadChi2_eta[whichhist[k]]->Fill(nBadChi2, eta);
		h2_nBadChi2_phi[whichhist[k]]->Fill(nBadChi2, phi);
		h2_nBadChi2_jet_ET[whichhist[k]]->Fill(nBadChi2, jet_ET);
		h2_nBadChi2_dphi[whichhist[k]]->Fill(nBadChi2, dphi);
		h2_nBadChi2_subjet_ET[whichhist[k]]->Fill(nBadChi2, subjet_ET);
		
		h2_maxTowDiff_subTowE[whichhist[k]]->Fill(maxTowDiff, subTowE);
		h2_maxTowDiff_maxTowE[whichhist[k]]->Fill(maxTowDiff, maxTowE);
		h2_maxTowDiff_ecc[whichhist[k]]->Fill(maxTowDiff, ecc);
		h2_maxTowDiff_chi2[whichhist[k]]->Fill(maxTowDiff, chi2);
		h2_maxTowDiff_frcoh[whichhist[k]]->Fill(maxTowDiff, frcoh);
		h2_maxTowDiff_frcem[whichhist[k]]->Fill(maxTowDiff, frcem);
		h2_maxTowDiff_eta[whichhist[k]]->Fill(maxTowDiff, eta);
		h2_maxTowDiff_phi[whichhist[k]]->Fill(maxTowDiff, phi);
		h2_maxTowDiff_jet_ET[whichhist[k]]->Fill(maxTowDiff, jet_ET);
		h2_maxTowDiff_dphi[whichhist[k]]->Fill(maxTowDiff, dphi);
		h2_maxTowDiff_subjet_ET[whichhist[k]]->Fill(maxTowDiff, subjet_ET);
		
		h2_subTowE_maxTowE[whichhist[k]]->Fill(subTowE, maxTowE);
		h2_subTowE_ecc[whichhist[k]]->Fill(subTowE, ecc);
		h2_subTowE_chi2[whichhist[k]]->Fill(subTowE, chi2);
		h2_subTowE_frcoh[whichhist[k]]->Fill(subTowE, frcoh);
		h2_subTowE_frcem[whichhist[k]]->Fill(subTowE, frcem);
		h2_subTowE_eta[whichhist[k]]->Fill(subTowE, eta);
		h2_subTowE_phi[whichhist[k]]->Fill(subTowE, phi);
		h2_subTowE_jet_ET[whichhist[k]]->Fill(subTowE, jet_ET);
		h2_subTowE_dphi[whichhist[k]]->Fill(subTowE, dphi);
		h2_subTowE_subjet_ET[whichhist[k]]->Fill(subTowE, subjet_ET);
		
		h2_maxTowE_ecc[whichhist[k]]->Fill(maxTowE, ecc);
		h2_maxTowE_chi2[whichhist[k]]->Fill(maxTowE, chi2);
		h2_maxTowE_frcoh[whichhist[k]]->Fill(maxTowE, frcoh);
		h2_maxTowE_frcem[whichhist[k]]->Fill(maxTowE, frcem);
		h2_maxTowE_eta[whichhist[k]]->Fill(maxTowE, eta);
		h2_maxTowE_phi[whichhist[k]]->Fill(maxTowE, phi);
		h2_maxTowE_jet_ET[whichhist[k]]->Fill(maxTowE, jet_ET);
		h2_maxTowE_dphi[whichhist[k]]->Fill(maxTowE, dphi);
		h2_maxTowE_subjet_ET[whichhist[k]]->Fill(maxTowE, subjet_ET);
		
		h2_ecc_chi2[whichhist[k]]->Fill(ecc, chi2);
		h2_ecc_frcoh[whichhist[k]]->Fill(ecc, frcoh);
		h2_ecc_frcem[whichhist[k]]->Fill(ecc, frcem);
		h2_ecc_eta[whichhist[k]]->Fill(ecc, eta);
		h2_ecc_phi[whichhist[k]]->Fill(ecc, phi);
		h2_ecc_jet_ET[whichhist[k]]->Fill(ecc, jet_ET);
		h2_ecc_dphi[whichhist[k]]->Fill(ecc, dphi);
		h2_ecc_subjet_ET[whichhist[k]]->Fill(ecc, subjet_ET);
		//cout << "filled first block" << endl;
		h2_chi2_frcoh[whichhist[k]]->Fill(chi2, frcoh);
		h2_chi2_frcem[whichhist[k]]->Fill(chi2, frcem);
		h2_chi2_eta[whichhist[k]]->Fill(chi2, eta);
		h2_chi2_phi[whichhist[k]]->Fill(chi2, phi);
		h2_chi2_jet_ET[whichhist[k]]->Fill(chi2, jet_ET);
		h2_chi2_dphi[whichhist[k]]->Fill(chi2, dphi);
		h2_chi2_subjet_ET[whichhist[k]]->Fill(chi2, subjet_ET);
		//cout << "filled second block" << endl;
		h2_frcoh_frcem[whichhist[k]]->Fill(frcoh, frcem);
		h2_frcoh_eta[whichhist[k]]->Fill(frcoh, eta);
		h2_frcoh_phi[whichhist[k]]->Fill(frcoh, phi);
		h2_frcoh_jet_ET[whichhist[k]]->Fill(frcoh, jet_ET);
		h2_frcoh_dphi[whichhist[k]]->Fill(frcoh, dphi);
		h2_frcoh_subjet_ET[whichhist[k]]->Fill(frcoh, subjet_ET);
		//cout << "filled third block" << endl;
		h2_frcem_eta[whichhist[k]]->Fill(frcem, eta);
		h2_frcem_phi[whichhist[k]]->Fill(frcem, phi);
		h2_frcem_jet_ET[whichhist[k]]->Fill(frcem, jet_ET);
		h2_frcem_dphi[whichhist[k]]->Fill(frcem, dphi);
		h2_frcem_subjet_ET[whichhist[k]]->Fill(frcem, subjet_ET);
		//cout << "filled fourth block" << endl;
		h2_eta_phi[whichhist[k]]->Fill(eta, phi);
		h2_eta_jet_ET[whichhist[k]]->Fill(eta, jet_ET);
		h2_eta_dphi[whichhist[k]]->Fill(eta, dphi);
		h2_eta_subjet_ET[whichhist[k]]->Fill(eta, subjet_ET);
		//cout << "filled fifth block" << endl;
		h2_phi_jet_ET[whichhist[k]]->Fill(phi, jet_ET);
		h2_phi_dphi[whichhist[k]]->Fill(phi, dphi);
		h2_phi_subjet_ET[whichhist[k]]->Fill(phi, subjet_ET);
		//cout << "filled sixth block" << endl;
		h2_jet_ET_dphi[whichhist[k]]->Fill(jet_ET, dphi);
		h2_jet_ET_subjet_ET[whichhist[k]]->Fill(jet_ET, subjet_ET);
		//cout << "filled seventh block" << endl;
		h2_subjet_ET_dphi[whichhist[k]]->Fill(subjet_ET, dphi);
	      }
	  }
    }
    //std:://cerr << "filled" << endl;
    // Save histograms to a file
    TFile *outputFile = TFile::Open((filebase+"_hist.root").c_str(), "RECREATE");
    // Write additional 2D histograms
    for(int i=0; i<numHistograms; ++i)
      {
	for(int j=0; j<numTypes; ++j)
	  {
	    //cerr << "write file" << i << " " << j << endl;
	    allhists[j][i]->Write();
	  }
      }

    int sd17 = get_scaledown17(runnumber);
    int sd10 = get_scaledown10(runnumber);
    int nmb = get_nmb(runnumber);
    //cout << sd17 << " " << sd10 << " " << nmb << endl;
    for(int i=0; i<nSpectra; ++i)
      {
	//cout << "jetSpectrum: " << jetSpectra[i] << endl;
	if(sd17 >= 0 && sd10 >= 0 && nmb > 0 && std::isfinite(sd17) && std::isfinite(sd10) && std::isfinite(nmb))
	  {
	    float effevt = (((1.*(1+sd10))/(1+sd17))*nmb);
	    //jetSpectra[i]->Scale(1./effevt);
	    //cout << "effevt: " << effevt << endl;
	    jetSpectra[i]->Write();
	  }
	else
	  {
	    jetSpectra[i]->Clear();
	    //cout << "PANIC IN " << runnumber << " " << filebase << " " << i << endl;
	  }
      }
    
    for(int i=0; i<nRatio; ++i)
      {
	forRatio[i]->Write();
      }
    xJ[0]->Write();
    xJ[1]->Write();
    xJ[2]->Write();
    xJ[3]->Write();
    for(int i=0; i<nzhist; ++i)
      {
	zhists[i]->Write();
      }
    for(int i=0; i<2; ++i)
      {
	h2_n2pc[i]->Write();
	h2_dPhiLayer[i]->Write();
      }


    for(int i=0; i<2; ++i)
      {
	h2_n2pcMinus[i]->Write();
	h2_n2pcPlus[i]->Write();
	h2_dPhiLayerMinus[i]->Write();
	h2_dPhiLayerPlus[i]->Write();
      }
    dijetCheckRat->Write();
    dijetCheckRatFull->Write();
    asdich_all->Write();
    asdich_fail->Write();
    asdich_dphi->Write();
    outputFile->Close();

    //cout << "wrote" << endl;
    // Clean up
    file->Close();
    //delete file;
    return 0;
  }
