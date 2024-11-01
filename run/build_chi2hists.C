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
      std::cerr << "Error opening file!" << std::endl;
      return 1;
    }

    // Get the tree
    TTree *jet_tree = (TTree*)file->Get("jet_tree");
    if (!jet_tree) {
      std::cerr << "Error getting jet_tree!" << std::endl;
      return 1;
    }

    // Declare variables for branches
    Float_t ecc, theta, frcoh, frcem, eta, phi, jet_ET, dphi, subjet_ET;
    int isdijet;
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

    // Create arrays for histogram types
    const char* baseNames2D[] = {"h2_ecc_layer", "h2_ecc_angle", "h2_ecc_E"};
    const char* baseNamesG20[] = {"h2_g20_ecc_angle", "h2_g20_ecc_frcoh", "h2_g20_ecc_frcem"};
    const char* baseNames1D[] = {"h1_jet_eta", "h1_jet_phi"};

    /*
    std::vector<TH2F*> histograms2D;
    std::vector<TH2F*> histogramsG20;
    std::vector<TH1F*> histograms1D;

    // Retrieve histograms from the file
    for (int h = 0; h < 3; h++) {
      for (int i = 0; i < 6; i++) {
	std::string histName = std::string(baseNames2D[h]) + std::to_string(h) + "_" + std::to_string(i);
	TH2F *hist = (TH2F*)file->Get(histName.c_str());
	if (hist) {
	  histograms2D.push_back(hist);
	}
      }
    }

    for (int h = 0; h < 3; h++) {
      for (int i = 0; i < 6; i++) {
	std::string histName = std::string(baseNamesG20[h]) + std::to_string(h) + "_" + std::to_string(i);
	TH2F *hist = (TH2F*)file->Get(histName.c_str());
	if (hist) {
	  histogramsG20.push_back(hist);
	}
      }
    }

    for (int h = 0; h < 2; h++) {
      for (int i = 0; i < 6; i++) {
	std::string histName = std::string(baseNames1D[h]) + std::to_string(h) + "_" + std::to_string(i);
	TH1F *hist = (TH1F*)file->Get(histName.c_str());
	if (hist) {
	  histograms1D.push_back(hist);
	}
      }
    }
    */
    // Create 2D histograms for all combinations of variables
    const int numHistograms = 2; // Number of histograms in each array
    const int numTypes = 35;
    // Arrays to hold the histograms
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
      "ecc_chi2", "ecc_frcoh", "ecc_frcem", "ecc_eta", "ecc_phi", "ecc_jet_ET", "ecc_dphi", "ecc_subjet_ET",
      "chi2_frcoh", "chi2_frcem", "chi2_eta", "chi2_phi", "chi2_jet_ET", "chi2_dphi", "chi2_subjet_ET",
      "frcoh_frcem", "frcoh_eta", "frcoh_phi", "frcoh_jet_ET", "frcoh_dphi", "frcoh_subjet_ET",
      "frcem_eta", "frcem_phi", "frcem_jet_ET", "frcem_dphi", "frcem_subjet_ET",
      "eta_phi", "eta_jet_ET", "eta_dphi", "eta_subjet_ET",
      "phi_jet_ET", "phi_dphi", "phi_subjet_ET",
      "jet_ET_dphi", "jet_ET_subjet_ET"
    };

    const char* titles[] = {
      "ecc vs chi2", "ecc vs frcoh", "ecc vs frcem", "ecc vs eta", "ecc vs phi", "ecc vs jet_ET", "ecc vs dphi", "ecc vs subjet_ET",
      "chi2 vs frcoh", "chi2 vs frcem", "chi2 vs eta", "chi2 vs phi", "chi2 vs jet_ET", "chi2 vs dphi", "chi2 vs subjet_ET",
      "frcoh vs frcem", "frcoh vs eta", "frcoh vs phi", "frcoh vs jet_ET", "frcoh vs dphi", "frcoh vs subjet_ET",
      "frcem vs eta", "frcem vs phi", "frcem vs jet_ET", "frcem vs dphi", "frcem vs subjet_ET",
      "eta vs phi", "eta vs jet_ET", "eta vs dphi", "eta vs subjet_ET",
      "phi vs jet_ET", "phi vs dphi", "phi vs subjet_ET",
      "jet_ET vs dphi", "jet_ET vs subjet_ET"
    };

    const double xbins[] = {
      1, 1, 1, 1, 1, 1, 1, 1,
      1e4, 1e4, 1e4, 1e4, 1e4, 1e4, 1e4,
      1.2, 1.2, 1.2, 1.2, 1.2, 1.2,
      1.2, 1.2, 1.2, 1.2, 1.2,
      1.5, 1.5, 1.5, 1.5,
      3.14, 3.14, 3.14,
      100, 100
    };

    const double ybins[] = {
      1e4, 1.2, 1.2, 1.5, 3.14, 100, 3.15, 100,
      1.2, 1.2, 1.5, 3.14, 100, 3.15, 100,
      1.2, 1.5, 3.14, 100, 3.15, 100,
      1.5, 3.14, 100, 3.15, 100,
      3.14, 100, 3.15, 100,
      100, 3.15, 100,
      3.15, 100
    };
    //std:://cout << "got all arrays for making hists" << endl;
    // Loop to create histograms
    for (int i = 0; i < numTypes; ++i) {
      for (int j = 0; j < numHistograms; ++j) {
	std::string name = "h2_" + std::string(names[i]) + "_" + std::to_string(j);
	std::string title = std::string(titles[i]) + " " + std::to_string(j);

        // Set ranges based on index
        double xMin = xbins[i]==1.5?-1.5:0;
	if(xbins[i]==3.14) xMin = -3.14;
        double xMax = xbins[i];
        double yMin = ybins[i]==1.5?-1.5:0;
	if(ybins[i]==3.14) yMin = -3.14;
        double yMax = ybins[i];

        // Create the histogram
        allhists[i][j] = new TH2F(name.c_str(), title.c_str(), 100, xMin, xMax, 100, yMin, yMax);
      }
    }
    //std:://cout << "made hists" << endl;
    // Loop over entries in the tree
    Long64_t nEntries = jet_tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        jet_tree->GetEntry(i);
        /*
        // Fill histograms based on your analysis
        for (auto hist : histograms2D) {
            hist->Fill(ecc, jet_ET); // Adjust as needed
        }
        for (auto hist : histogramsG20) {
            hist->Fill(frcoh, frcem); // Adjust as needed
        }
        for (auto hist : histograms1D) {
            if (hist->GetName() == std::string("h1_jet_eta")) {
                hist->Fill(eta);
            } else if (hist->GetName() == std::string("h1_jet_phi")) {
                hist->Fill(phi);
            }
        }
	*/
        // Fill additional 2D histograms
	//cout << "filling (got entry). isdijet = " << isdijet << endl;
	float chi2 = maxTowChi2[0];
	if(maxTowChi2[1] > chi2) chi2 = maxTowChi2[1];
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
    //cout << "done filling" << endl;
    // Save histograms to a file
    TFile *outputFile = TFile::Open((filebase+"_hist.root").c_str(), "RECREATE");
    /*
    for (auto hist : histograms2D) {
      hist->Write();
    }
    for (auto hist : histogramsG20) {
      hist->Write();
    }
    for (auto hist : histograms1D) {
      hist->Write();
    }
    */
    // Write additional 2D histograms
    for(int i=0; i<numHistograms; ++i)
      {
    h2_ecc_chi2[i]->Write();
    h2_ecc_frcoh[i]->Write();
    h2_ecc_frcem[i]->Write();
    h2_ecc_eta[i]->Write();
    h2_ecc_phi[i]->Write();
    h2_ecc_jet_ET[i]->Write();
    h2_ecc_dphi[i]->Write();
    h2_ecc_subjet_ET[i]->Write();

    h2_chi2_frcoh[i]->Write();
    h2_chi2_frcem[i]->Write();
    h2_chi2_eta[i]->Write();
    h2_chi2_phi[i]->Write();
    h2_chi2_jet_ET[i]->Write();
    h2_chi2_dphi[i]->Write();
    h2_chi2_subjet_ET[i]->Write();

    h2_frcoh_frcem[i]->Write();
    h2_frcoh_eta[i]->Write();
    h2_frcoh_phi[i]->Write();
    h2_frcoh_jet_ET[i]->Write();
    h2_frcoh_dphi[i]->Write();
    h2_frcoh_subjet_ET[i]->Write();

    h2_frcem_eta[i]->Write();
    h2_frcem_phi[i]->Write();
    h2_frcem_jet_ET[i]->Write();
    h2_frcem_dphi[i]->Write();
    h2_frcem_subjet_ET[i]->Write();

    h2_eta_phi[i]->Write();
    h2_eta_jet_ET[i]->Write();
    h2_eta_dphi[i]->Write();
    h2_eta_subjet_ET[i]->Write();

    h2_phi_jet_ET[i]->Write();
    h2_phi_dphi[i]->Write();
    h2_phi_subjet_ET[i]->Write();

    h2_jet_ET_dphi[i]->Write();
    h2_jet_ET_subjet_ET[i]->Write();
      }
    outputFile->Close();
    //cout << "wrote" << endl;
    // Clean up
    file->Close();
    delete file;
    return 0;
  }
