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

  cout << runnumber << endl;
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
  cout << endl;
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
    Float_t maxTowChi2[3];
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
    //std:://cerr << "set branches" << endl;
    // Create 2D histograms for all combinations of variables
    const int numHistograms = 24; // Number of histograms in each array
    const int numTypes = 92;
    // Arrays to hold the histograms
    cout << "got branches" << endl;
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
    cout << "made all hist pointers" << endl;
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
    cout << "put hists in array" << endl;
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
    std::cout << "got all arrays for making hists" << endl;
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
    std::cout << "made hists" << endl;
    // Loop over entries in the tree
    const int nRatio = 7;
    TH1F* forRatio[nRatio];
    for(int i=0; i<nRatio; ++i)
      {
	forRatio[i] = new TH1F(("h1_forRatio_"+to_string(i)).c_str(),"",i<5?1000:100,0,i<5?100:1);
      }
    const int nSpectra = 33;
    TH1F* jetSpectra[nSpectra];
    TH1F* xJ[2];
    xJ[0] = new TH1F("xJ0","",100,0,1);
    xJ[1] = new TH1F("xJ1","",100,0,1);
    for(int i=0; i<nSpectra; ++i)
      {
	jetSpectra[i] = new TH1F(("h1_jetSpectra_"+to_string(i)).c_str(),"",1000,0,100);
	cout << jetSpectra[i] << endl;
      }
    bool cutArr[nSpectra];
    Long64_t nEntries = jet_tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        jet_tree->GetEntry(i);
	if(abs(zvtx) > 30) continue;
	if(abs(eta) > 0.7) continue;
	bool stripCut = (dphi < 3*M_PI/4 && isdijet); //(1-frcem-frcoh) > ((2.0/3.0)*frcoh);//((elmbgvec >> 4) & 1);
	bool dhCut = ((bbfqavec >> 5) & 1);
	bool dPhiCut = (frcem+frcoh) < 0.7;
	bool ETCut = ((frcem < 0.1) && (jet_ET > (50*frcem+20))) && (stripCut || !isdijet);
	bool ZSCut = ((frcem > 0.9) && (jet_ET > (-50*frcem+75))) && (stripCut || !isdijet);
	bool chi2cut = jet_ET > 25 && maxETowChi2 < 10;
	cutArr[0]=false;
	cutArr[1]=stripCut;
	cutArr[2]=dhCut;
	cutArr[3]=dPhiCut;
	cutArr[4]=ETCut;
	cutArr[5]=ZSCut;
	cutArr[6]=(stripCut || dhCut);
	cutArr[7]=(stripCut || dPhiCut);
	cutArr[8]=(stripCut || ETCut);
	cutArr[9]=(stripCut || ZSCut);
	cutArr[10]=(dhCut || dPhiCut);
	cutArr[11]=(dhCut || ETCut);
	cutArr[12]=(dhCut || ZSCut);
	cutArr[13]=(dPhiCut || ETCut);
	cutArr[14]=(dPhiCut || ZSCut);
	cutArr[15]=(ETCut || ZSCut);
	cutArr[16]=(stripCut || dhCut || dPhiCut);
	cutArr[17]=(stripCut || dhCut || ETCut);
	cutArr[18]=(stripCut || dhCut || ZSCut);
	cutArr[19]=(stripCut || dPhiCut || ETCut);
	cutArr[20]=(stripCut || dPhiCut || ZSCut);
	cutArr[21]=(stripCut || ETCut || ZSCut);
	cutArr[22]=(dhCut || dPhiCut || ETCut);
	cutArr[23]=(dhCut || dPhiCut || ZSCut);
	cutArr[24]=(dhCut || ETCut || ZSCut);
	cutArr[25]=(dPhiCut || ETCut || ZSCut);
	cutArr[26]=(stripCut || dhCut || dPhiCut || ETCut);
	cutArr[27]=(stripCut || dhCut || dPhiCut || ZSCut);
	cutArr[28]=(stripCut || dhCut || ETCut || ZSCut);
	cutArr[29]=(stripCut || dPhiCut || ETCut || ZSCut);
	cutArr[30]=(dhCut || dPhiCut || ETCut || ZSCut || stripCut);
	cutArr[31]=(stripCut || dhCut || dPhiCut || ETCut || ZSCut);// || chi2cut);
	cutArr[32]=chi2cut;
	for(int j=0; j<nSpectra; ++j)
	  {
	    if(!cutArr[j])
	      {
		jetSpectra[j]->Fill(jet_ET);
		if(isdijet) jetSpectra[j]->Fill(subjet_ET);
		
	      }
	  }
	forRatio[0]->Fill(jet_ET);
	if(isdijet) forRatio[1]->Fill(jet_ET);
	if(!cutArr[1] && isdijet) forRatio[2]->Fill(jet_ET);
	if(!cutArr[31]) forRatio[3]->Fill(jet_ET);
	if(!cutArr[31] && isdijet) forRatio[4]->Fill(jet_ET);
	if(isdijet && jet_ET > 20 && subjet_ET > 10)
	  {
	    forRatio[5]->Fill((jet_ET-subjet_ET)/(jet_ET+subjet_ET));
	    xJ[0]->Fill(subjet_ET/jet_ET);
	  }
	if(isdijet && !cutArr[31] && jet_ET > 20 && subjet_ET > 10)
	  {
	    forRatio[6]->Fill((jet_ET-subjet_ET)/(jet_ET+subjet_ET));
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
	int whichhist = -1;
	int threshes[numHistograms/4] = {8,15,20,25,35,40};
	int slt[numHistograms/4] = {8,8,10,10,15,15};
	for(int j=0; j<numHistograms/4; ++j)
	  {
	    if(jet_ET < threshes[j] && j < numHistograms/4) continue;
	    whichhist = j + (isdijet?numHistograms/4:0); //(cutArr[31]?numHistograms/4:0); //((subjet_ET > slt[j] && isdijet)?numHistograms/4:0);// + ();

	if(whichhist < 0 || whichhist > numHistograms-1) continue;


	if(isdijet)
	  {
	    h2_AJ_dphi[whichhist]->Fill((jet_ET-subjet_ET)/(jet_ET+subjet_ET),dphi);
	  }
	//if(maxTowE > 15 && maxETowChi2 < 100) cout << "Detector: " << det << " tower ET: " << maxTowE << " chi2: " << maxETowChi2 << " is ZS: " << maxETowIsZS <<endl;
	//if(cutArr[31]) continue;
	h2_maxETowChi2_nBadChi2[whichhist]->Fill(maxETowChi2, nBadChi2);
	h2_maxETowChi2_maxTowDiff[whichhist]->Fill(maxETowChi2, maxTowDiff);
	h2_maxETowChi2_subTowE[whichhist]->Fill(maxETowChi2, subTowE);
	h2_maxETowChi2_maxTowE[whichhist]->Fill(maxETowChi2, maxTowE);
	h2_maxETowChi2_ecc[whichhist]->Fill(maxETowChi2, ecc);
	h2_maxETowChi2_chi2[whichhist]->Fill(maxETowChi2, chi2);
        h2_maxETowChi2_frcoh[whichhist]->Fill(maxETowChi2, frcoh);
        h2_maxETowChi2_frcem[whichhist]->Fill(maxETowChi2, frcem);
        h2_maxETowChi2_eta[whichhist]->Fill(maxETowChi2, eta);
        h2_maxETowChi2_phi[whichhist]->Fill(maxETowChi2, phi);
        h2_maxETowChi2_jet_ET[whichhist]->Fill(maxETowChi2, jet_ET);
        h2_maxETowChi2_dphi[whichhist]->Fill(maxETowChi2, dphi);
        h2_maxETowChi2_subjet_ET[whichhist]->Fill(maxETowChi2, subjet_ET);
	
	h2_nBadChi2_maxTowDiff[whichhist]->Fill(nBadChi2, maxTowDiff);
	h2_nBadChi2_subTowE[whichhist]->Fill(nBadChi2, subTowE);
	h2_nBadChi2_maxTowE[whichhist]->Fill(nBadChi2, maxTowE);
	h2_nBadChi2_ecc[whichhist]->Fill(nBadChi2, ecc);
	h2_nBadChi2_chi2[whichhist]->Fill(nBadChi2, chi2);
        h2_nBadChi2_frcoh[whichhist]->Fill(nBadChi2, frcoh);
        h2_nBadChi2_frcem[whichhist]->Fill(nBadChi2, frcem);
        h2_nBadChi2_eta[whichhist]->Fill(nBadChi2, eta);
        h2_nBadChi2_phi[whichhist]->Fill(nBadChi2, phi);
        h2_nBadChi2_jet_ET[whichhist]->Fill(nBadChi2, jet_ET);
        h2_nBadChi2_dphi[whichhist]->Fill(nBadChi2, dphi);
        h2_nBadChi2_subjet_ET[whichhist]->Fill(nBadChi2, subjet_ET);
	
	h2_maxTowDiff_subTowE[whichhist]->Fill(maxTowDiff, subTowE);
	h2_maxTowDiff_maxTowE[whichhist]->Fill(maxTowDiff, maxTowE);
	h2_maxTowDiff_ecc[whichhist]->Fill(maxTowDiff, ecc);
	h2_maxTowDiff_chi2[whichhist]->Fill(maxTowDiff, chi2);
        h2_maxTowDiff_frcoh[whichhist]->Fill(maxTowDiff, frcoh);
        h2_maxTowDiff_frcem[whichhist]->Fill(maxTowDiff, frcem);
        h2_maxTowDiff_eta[whichhist]->Fill(maxTowDiff, eta);
        h2_maxTowDiff_phi[whichhist]->Fill(maxTowDiff, phi);
        h2_maxTowDiff_jet_ET[whichhist]->Fill(maxTowDiff, jet_ET);
        h2_maxTowDiff_dphi[whichhist]->Fill(maxTowDiff, dphi);
        h2_maxTowDiff_subjet_ET[whichhist]->Fill(maxTowDiff, subjet_ET);

	h2_subTowE_maxTowE[whichhist]->Fill(subTowE, maxTowE);
	h2_subTowE_ecc[whichhist]->Fill(subTowE, ecc);
	h2_subTowE_chi2[whichhist]->Fill(subTowE, chi2);
        h2_subTowE_frcoh[whichhist]->Fill(subTowE, frcoh);
        h2_subTowE_frcem[whichhist]->Fill(subTowE, frcem);
        h2_subTowE_eta[whichhist]->Fill(subTowE, eta);
        h2_subTowE_phi[whichhist]->Fill(subTowE, phi);
        h2_subTowE_jet_ET[whichhist]->Fill(subTowE, jet_ET);
        h2_subTowE_dphi[whichhist]->Fill(subTowE, dphi);
        h2_subTowE_subjet_ET[whichhist]->Fill(subTowE, subjet_ET);
	
	h2_maxTowE_ecc[whichhist]->Fill(maxTowE, ecc);
	h2_maxTowE_chi2[whichhist]->Fill(maxTowE, chi2);
        h2_maxTowE_frcoh[whichhist]->Fill(maxTowE, frcoh);
        h2_maxTowE_frcem[whichhist]->Fill(maxTowE, frcem);
        h2_maxTowE_eta[whichhist]->Fill(maxTowE, eta);
        h2_maxTowE_phi[whichhist]->Fill(maxTowE, phi);
        h2_maxTowE_jet_ET[whichhist]->Fill(maxTowE, jet_ET);
        h2_maxTowE_dphi[whichhist]->Fill(maxTowE, dphi);
        h2_maxTowE_subjet_ET[whichhist]->Fill(maxTowE, subjet_ET);

        h2_ecc_chi2[whichhist]->Fill(ecc, chi2);
        h2_ecc_frcoh[whichhist]->Fill(ecc, frcoh);
        h2_ecc_frcem[whichhist]->Fill(ecc, frcem);
        h2_ecc_eta[whichhist]->Fill(ecc, eta);
        h2_ecc_phi[whichhist]->Fill(ecc, phi);
        h2_ecc_jet_ET[whichhist]->Fill(ecc, jet_ET);
        h2_ecc_dphi[whichhist]->Fill(ecc, dphi);
        h2_ecc_subjet_ET[whichhist]->Fill(ecc, subjet_ET);
	//cout << "filled first block" << endl;
        h2_chi2_frcoh[whichhist]->Fill(chi2, frcoh);
        h2_chi2_frcem[whichhist]->Fill(chi2, frcem);
        h2_chi2_eta[whichhist]->Fill(chi2, eta);
        h2_chi2_phi[whichhist]->Fill(chi2, phi);
        h2_chi2_jet_ET[whichhist]->Fill(chi2, jet_ET);
        h2_chi2_dphi[whichhist]->Fill(chi2, dphi);
        h2_chi2_subjet_ET[whichhist]->Fill(chi2, subjet_ET);
	//cout << "filled second block" << endl;
        h2_frcoh_frcem[whichhist]->Fill(frcoh, frcem);
	if(whichhist < numHistograms/4) h2_frcoh_frcem[whichhist+numHistograms/4]->Fill(frcoh,frcem);
        h2_frcoh_eta[whichhist]->Fill(frcoh, eta);
        h2_frcoh_phi[whichhist]->Fill(frcoh, phi);
        h2_frcoh_jet_ET[whichhist]->Fill(frcoh, jet_ET);
        h2_frcoh_dphi[whichhist]->Fill(frcoh, dphi);
        h2_frcoh_subjet_ET[whichhist]->Fill(frcoh, subjet_ET);
	//cout << "filled third block" << endl;
        h2_frcem_eta[whichhist]->Fill(frcem, eta);
        h2_frcem_phi[whichhist]->Fill(frcem, phi);
        h2_frcem_jet_ET[whichhist]->Fill(frcem, jet_ET);
        h2_frcem_dphi[whichhist]->Fill(frcem, dphi);
        h2_frcem_subjet_ET[whichhist]->Fill(frcem, subjet_ET);
	//cout << "filled fourth block" << endl;
        h2_eta_phi[whichhist]->Fill(eta, phi);
	if(whichhist < numHistograms/4) h2_eta_phi[whichhist+numHistograms/4]->Fill(eta,phi);
        h2_eta_jet_ET[whichhist]->Fill(eta, jet_ET);
        h2_eta_dphi[whichhist]->Fill(eta, dphi);
        h2_eta_subjet_ET[whichhist]->Fill(eta, subjet_ET);
	//cout << "filled fifth block" << endl;
        h2_phi_jet_ET[whichhist]->Fill(phi, jet_ET);
        h2_phi_dphi[whichhist]->Fill(phi, dphi);
        h2_phi_subjet_ET[whichhist]->Fill(phi, subjet_ET);
	//cout << "filled sixth block" << endl;
        h2_jet_ET_dphi[whichhist]->Fill(jet_ET, dphi);
	if(whichhist < numHistograms/4) h2_jet_ET_dphi[whichhist+numHistograms/4]->Fill(jet_ET,dphi);
        h2_jet_ET_subjet_ET[whichhist]->Fill(jet_ET, subjet_ET);
	//cout << "filled seventh block" << endl;
	h2_subjet_ET_dphi[whichhist]->Fill(subjet_ET, dphi);
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
    cout << sd17 << " " << sd10 << " " << nmb << endl;
    for(int i=0; i<nSpectra; ++i)
      {
	cout << "jetSpectrum: " << jetSpectra[i] << endl;
	if(sd17 >= 0 && sd10 >= 0 && nmb > 0 && std::isfinite(sd17) && std::isfinite(sd10) && std::isfinite(nmb))
	  {
	    float effevt = (((1.*(1+sd10))/(1+sd17))*nmb);
	    //jetSpectra[i]->Scale(1./effevt);
	    cout << "effevt: " << effevt << endl;
	    jetSpectra[i]->Write();
	  }
	else
	  {
	    jetSpectra[i]->Clear();
	    cout << "PANIC IN " << runnumber << " " << filebase << " " << i << endl;
	  }
      }
    
    for(int i=0; i<nRatio; ++i)
      {
	forRatio[i]->Write();
      }
    xJ[0]->Write();
    xJ[1]->Write();
    outputFile->Close();

    //cout << "wrote" << endl;
    // Clean up
    file->Close();
    delete file;
    return 0;
  }
