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

int build_chi2hists(string filebase)
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
    Float_t ecc, theta, frcoh, frcem, eta, phi, jet_ET, dphi, subjet_ET, maxTowE, subTowE, maxTowDiff, maxETowChi2;
    int isdijet, nBadChi2, maxETowChi2Det, maxETowIsZS;
    Float_t jetcompE[3][512], jetcompEta[3][512], jetcompPhi[3][512];
    Float_t maxTowChi2[3];
    
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
    //std:://cerr << "set branches" << endl;
    // Create 2D histograms for all combinations of variables
    const int numHistograms = 2; // Number of histograms in each array
    const int numTypes = 90;
    // Arrays to hold the histograms

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
      h2_jet_ET_dphi, h2_jet_ET_subjet_ET
    };

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
      "jet_ET_dphi", "jet_ET_subjet_ET"
    };

    const double xbins[] = {
      1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4,
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
      100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
      100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
      100, 100, 100, 100, 100, 100, 100, 100, 100,
      1, 1, 1, 1, 1, 1, 1, 1,
      100, 100, 100, 100, 100, 100, 100,
      1.25, 1.25, 1.25, 1.25, 1.25, 1.25,
      1.25, 1.25, 1.25, 1.25, 1.25,
      1.5, 1.5, 1.5, 1.5,
      3.14, 3.14, 3.14,
      100, 100
    };

    const double ybins[] = {
      64, 100, 100, 1.25, 1.25, 1.5, 3.14, 100, 3.15, 100, 1, 100, 100,
      100, 100, 1.25, 1.25, 1.5, 3.14, 100, 3.15, 100, 1, 100, 100,
      100, 1.25, 1.25, 1.5, 3.14, 100, 3.15, 100, 1, 100, 100,
      100, 1.25, 1.25, 1.5, 3.14, 100, 3.15, 100, 1, 100,
      100, 1.25, 1.25, 1.5, 3.14, 100, 3.15, 100, 1,
      100, 1.25, 1.25, 1.5, 3.14, 100, 3.15, 100,
      1.25, 1.25, 1.5, 3.14, 100, 3.15, 100,
      1.25, 1.5, 3.14, 100, 3.15, 100,
      1.5, 3.14, 100, 3.15, 100,
      3.14, 100, 3.15, 100,
      100, 3.15, 100,
      3.15, 100
    };
    //std:://cerr << "prep to make hists" << std::endl;
    //std:://cout << "got all arrays for making hists" << endl;
    // Loop to create histograms
    for (int i = 0; i < numTypes; ++i) {
      for (int j = 0; j < numHistograms; ++j) {
	std::string name = "h2_" + std::string(names[i]) + "_" + std::to_string(j);
	std::string title = ""; //std::string(titles[i]) + " " + std::to_string(j);

        // Set ranges based on index
        double xMin = xbins[i]==1.5?-1.5:0;
	if(xbins[i]==3.14) xMin = -3.14;
	if(xbins[i]==1.25) xMin = -0.25;
        double xMax = xbins[i];
        double yMin = ybins[i]==1.5?-1.5:0;
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
    //std:://cout << "made hists" << endl;
    // Loop over entries in the tree
    Long64_t nEntries = jet_tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        jet_tree->GetEntry(i);
	//cout << "filling (got entry). isdijet = " << isdijet << endl;
	float chi2 = maxTowChi2[0];
	if(maxTowChi2[2] > chi2) chi2 = maxTowChi2[2];
	string det;
	if(maxETowChi2Det == 0) det = "EMCal";
	else if(maxETowChi2Det == 1) det = "IHCal";
	else if(maxETowChi2Det == 2) det = "OHCal";
	else det = "Error";
	if(maxTowE > 15 && maxETowChi2 < 100) cout << "Detector: " << det << " tower ET: " << maxTowE << " chi2: " << maxETowChi2 << " is ZS: " << maxETowIsZS <<endl;
	h2_maxETowChi2_nBadChi2[isdijet]->Fill(maxETowChi2, nBadChi2);
	h2_maxETowChi2_maxTowDiff[isdijet]->Fill(maxETowChi2, maxTowDiff);
	h2_maxETowChi2_subTowE[isdijet]->Fill(maxETowChi2, subTowE);
	h2_maxETowChi2_maxTowE[isdijet]->Fill(maxETowChi2, maxTowE);
	h2_maxETowChi2_ecc[isdijet]->Fill(maxETowChi2, ecc);
	h2_maxETowChi2_chi2[isdijet]->Fill(maxETowChi2, chi2);
        h2_maxETowChi2_frcoh[isdijet]->Fill(maxETowChi2, frcoh);
        h2_maxETowChi2_frcem[isdijet]->Fill(maxETowChi2, frcem);
        h2_maxETowChi2_eta[isdijet]->Fill(maxETowChi2, eta);
        h2_maxETowChi2_phi[isdijet]->Fill(maxETowChi2, phi);
        h2_maxETowChi2_jet_ET[isdijet]->Fill(maxETowChi2, jet_ET);
        h2_maxETowChi2_dphi[isdijet]->Fill(maxETowChi2, dphi);
        h2_maxETowChi2_subjet_ET[isdijet]->Fill(maxETowChi2, subjet_ET);
	
	h2_nBadChi2_maxTowDiff[isdijet]->Fill(nBadChi2, maxTowDiff);
	h2_nBadChi2_subTowE[isdijet]->Fill(nBadChi2, subTowE);
	h2_nBadChi2_maxTowE[isdijet]->Fill(nBadChi2, maxTowE);
	h2_nBadChi2_ecc[isdijet]->Fill(nBadChi2, ecc);
	h2_nBadChi2_chi2[isdijet]->Fill(nBadChi2, chi2);
        h2_nBadChi2_frcoh[isdijet]->Fill(nBadChi2, frcoh);
        h2_nBadChi2_frcem[isdijet]->Fill(nBadChi2, frcem);
        h2_nBadChi2_eta[isdijet]->Fill(nBadChi2, eta);
        h2_nBadChi2_phi[isdijet]->Fill(nBadChi2, phi);
        h2_nBadChi2_jet_ET[isdijet]->Fill(nBadChi2, jet_ET);
        h2_nBadChi2_dphi[isdijet]->Fill(nBadChi2, dphi);
        h2_nBadChi2_subjet_ET[isdijet]->Fill(nBadChi2, subjet_ET);
	
	h2_maxTowDiff_subTowE[isdijet]->Fill(maxTowDiff, subTowE);
	h2_maxTowDiff_maxTowE[isdijet]->Fill(maxTowDiff, maxTowE);
	h2_maxTowDiff_ecc[isdijet]->Fill(maxTowDiff, ecc);
	h2_maxTowDiff_chi2[isdijet]->Fill(maxTowDiff, chi2);
        h2_maxTowDiff_frcoh[isdijet]->Fill(maxTowDiff, frcoh);
        h2_maxTowDiff_frcem[isdijet]->Fill(maxTowDiff, frcem);
        h2_maxTowDiff_eta[isdijet]->Fill(maxTowDiff, eta);
        h2_maxTowDiff_phi[isdijet]->Fill(maxTowDiff, phi);
        h2_maxTowDiff_jet_ET[isdijet]->Fill(maxTowDiff, jet_ET);
        h2_maxTowDiff_dphi[isdijet]->Fill(maxTowDiff, dphi);
        h2_maxTowDiff_subjet_ET[isdijet]->Fill(maxTowDiff, subjet_ET);

	h2_subTowE_maxTowE[isdijet]->Fill(subTowE, maxTowE);
	h2_subTowE_ecc[isdijet]->Fill(subTowE, ecc);
	h2_subTowE_chi2[isdijet]->Fill(subTowE, chi2);
        h2_subTowE_frcoh[isdijet]->Fill(subTowE, frcoh);
        h2_subTowE_frcem[isdijet]->Fill(subTowE, frcem);
        h2_subTowE_eta[isdijet]->Fill(subTowE, eta);
        h2_subTowE_phi[isdijet]->Fill(subTowE, phi);
        h2_subTowE_jet_ET[isdijet]->Fill(subTowE, jet_ET);
        h2_subTowE_dphi[isdijet]->Fill(subTowE, dphi);
        h2_subTowE_subjet_ET[isdijet]->Fill(subTowE, subjet_ET);
	
	h2_maxTowE_ecc[isdijet]->Fill(maxTowE, ecc);
	h2_maxTowE_chi2[isdijet]->Fill(maxTowE, chi2);
        h2_maxTowE_frcoh[isdijet]->Fill(maxTowE, frcoh);
        h2_maxTowE_frcem[isdijet]->Fill(maxTowE, frcem);
        h2_maxTowE_eta[isdijet]->Fill(maxTowE, eta);
        h2_maxTowE_phi[isdijet]->Fill(maxTowE, phi);
        h2_maxTowE_jet_ET[isdijet]->Fill(maxTowE, jet_ET);
        h2_maxTowE_dphi[isdijet]->Fill(maxTowE, dphi);
        h2_maxTowE_subjet_ET[isdijet]->Fill(maxTowE, subjet_ET);

        h2_ecc_chi2[isdijet]->Fill(ecc, chi2);
        h2_ecc_frcoh[isdijet]->Fill(ecc, frcoh);
        h2_ecc_frcem[isdijet]->Fill(ecc, frcem);
        h2_ecc_eta[isdijet]->Fill(ecc, eta);
        h2_ecc_phi[isdijet]->Fill(ecc, phi);
        h2_ecc_jet_ET[isdijet]->Fill(ecc, jet_ET);
        h2_ecc_dphi[isdijet]->Fill(ecc, dphi);
        h2_ecc_subjet_ET[isdijet]->Fill(ecc, subjet_ET);
	//cout << "filled first block" << endl;
        h2_chi2_frcoh[isdijet]->Fill(chi2, frcoh);
        h2_chi2_frcem[isdijet]->Fill(chi2, frcem);
        h2_chi2_eta[isdijet]->Fill(chi2, eta);
        h2_chi2_phi[isdijet]->Fill(chi2, phi);
        h2_chi2_jet_ET[isdijet]->Fill(chi2, jet_ET);
        h2_chi2_dphi[isdijet]->Fill(chi2, dphi);
        h2_chi2_subjet_ET[isdijet]->Fill(chi2, subjet_ET);
	//cout << "filled second block" << endl;
        h2_frcoh_frcem[isdijet]->Fill(frcoh, frcem);
        h2_frcoh_eta[isdijet]->Fill(frcoh, eta);
        h2_frcoh_phi[isdijet]->Fill(frcoh, phi);
        h2_frcoh_jet_ET[isdijet]->Fill(frcoh, jet_ET);
        h2_frcoh_dphi[isdijet]->Fill(frcoh, dphi);
        h2_frcoh_subjet_ET[isdijet]->Fill(frcoh, subjet_ET);
	//cout << "filled third block" << endl;
        h2_frcem_eta[isdijet]->Fill(frcem, eta);
        h2_frcem_phi[isdijet]->Fill(frcem, phi);
        h2_frcem_jet_ET[isdijet]->Fill(frcem, jet_ET);
        h2_frcem_dphi[isdijet]->Fill(frcem, dphi);
        h2_frcem_subjet_ET[isdijet]->Fill(frcem, subjet_ET);
	//cout << "filled fourth block" << endl;
        h2_eta_phi[isdijet]->Fill(eta, phi);
        h2_eta_jet_ET[isdijet]->Fill(eta, jet_ET);
        h2_eta_dphi[isdijet]->Fill(eta, dphi);
        h2_eta_subjet_ET[isdijet]->Fill(eta, subjet_ET);
	//cout << "filled fifth block" << endl;
        h2_phi_jet_ET[isdijet]->Fill(phi, jet_ET);
        h2_phi_dphi[isdijet]->Fill(phi, dphi);
        h2_phi_subjet_ET[isdijet]->Fill(phi, subjet_ET);
	//cout << "filled sixth block" << endl;
        h2_jet_ET_dphi[isdijet]->Fill(jet_ET, dphi);
        h2_jet_ET_subjet_ET[isdijet]->Fill(jet_ET, subjet_ET);
	//cout << "filled seventh block" << endl;
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
    outputFile->Close();
    //cout << "wrote" << endl;
    // Clean up
    file->Close();
    delete file;
    return 0;
  }
