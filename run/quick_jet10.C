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
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TAxis.h>
#include <TLine.h>
#include <TFile.h>
#include <algorithm>
#include "dlUtility.h"


int quick_jet10(string filebase="", int njob=0, int dotow = 0)
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
  //gROOT->SetStyle("Plain");
  //SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0);

  string filename=filebase;
  TChain* tree;
  ifstream list;
  string line;
  list.open(filename,ifstream::in);
  string runnum = filename.substr(6,5);
  string simliststr = "lists/sim.imagelist";
  string simstr = "sim";
  string datstr = "dat";
  string idstr;
  if(filename == simliststr)
    {
      runnum = "sim";
      idstr = simstr;
    }
  else
    {
      idstr = datstr;
    }
  int counter = 0;
  tree = new TChain("ttree");
  if(!list)
    {
      cout << "nosimlist" << endl;
      exit(1);
    }
  
  for(int i=0; i<10; ++i)//10*njob; i<10*(njob+1); ++i)
    {
      int breakit = 0;
      getline(list,line);
      if(list.eof()) breakit = 1;
      try
	{
	  tree->Add(line.c_str());
	}
      catch(...)
	{
	  continue;
	}
      if(breakit) break;
    }
  //cout << "test-1" << endl;
  long unsigned int bbfqavec;
  int njet;
  float vtx[3];
  float jet_e[100];
  float jet_et[100];
  float frcem[100];
  float jet_ph[100];
  tree->SetBranchAddress("bbfqavec",&bbfqavec);
  tree->SetBranchAddress("vtx",vtx);
  tree->SetBranchAddress("njet",&njet);
  tree->SetBranchAddress("frcem",frcem);
  tree->SetBranchAddress("jet_ph",jet_ph);
  tree->SetBranchAddress("jet_e",jet_e);
  tree->SetBranchAddress("jet_et",jet_et);
  cout << "end getting branches" << endl;
  
  TFile* jetfile = TFile::Open(("output/simhists/run_jet10_"+idstr+"_"+to_string(njob)+"_simhists.root").c_str(),"RECREATE");
  

  int eventbase = {0};
  int prevt = 0;
  const int nh2 = 6;
  TH2F* hists2[nh2];
  const int nz = 2;
  TH1F* h1_zdist[nz];
  for(int i=0; i<nz; ++i)
    {
      h1_zdist[i] = new TH1F(("h1_zdist"+to_string(i)).c_str(),"",200,-150,150);
    }
  TH1F* h1_forrat[5];
  for(int i=0; i<5; ++i)
    {
      h1_forrat[i] = new TH1F(("h2_forrat_tosee"+to_string(i)).c_str(),"",i<3?1000:100,0,i<3?100:1);
    }

  float xlo[nh2] = {-0.2,-0.2,-0.2,-0.2,-0.2,-0.2};
  float xhi[nh2] = {1.2,1.2,1.2,1.2,1.2,1.2};
  float ylo[nh2] = {0,0,0,0,-0.2,-0.2};
  float yhi[nh2] = {100,100,100,100,1.2,1.2};
  
  for(int i=0; i<nh2; ++i)
    {
      hists2[i] = new TH2F(("hists2"+to_string(i)).c_str(),"",100,xlo[i],xhi[i],100,ylo[i],yhi[i]);
    }
  
  TH1F* h1_ucspec = new TH1F("h1_ucspec","",1000,0,100);
  TH1F* h1_cspec = new TH1F("h1_cspec","",1000,0,100);
  TH1F* xJ[2];
  xJ[0]=new TH1F("xJ0_sim","",100,0,1);
  xJ[1]=new TH1F("xJ1_sim","",100,0,1);

  cout << "start processing: " << endl;
  for(int i=0; i<tree->GetEntries(); ++i)
    {

      tree->GetEntry(i);   
      

      if(abs(vtx[2]) > 150)
	{
	  continue;
	}
      float ljetET = 0;
      float subjetET = 0;
      float ljetph = 0;
      float subjetph = 0;
      float ljetfrcem = -1;
      int isdijet = 0;
      float ljeteta = 0;
      for(int j=0; j<njet; ++j)
	{
	  if(jet_e[j] > ljetET)
	    {
	      subjetET = ljetET;
	      subjetph = ljetph;
	      ljetET = jet_e[j];
	      ljetph = jet_ph[j];
	      ljetfrcem = frcem[j];
	      ljeteta = jet_et[j];
	    }
	  if(jet_e[j] > 8)
	    {
	      hists2[0]->Fill(ljetfrcem,jet_e[j]);
	      //hists[1]->Fill(
	    }
	  else continue;
	  if(abs(jet_et[j]) > 0.7) continue;
	  h1_ucspec->Fill(jet_e[j]);
	  
	}
      if(subjetET > 8) isdijet = 1;
      
      float dphi = abs(ljetph - subjetph);
      if(dphi > M_PI) dphi = 2*M_PI - dphi;

      bool hdPhiCut = dphi < 3*M_PI/4 && isdijet;
      bool bbCut = (bbfqavec >> 5) & 1;
      //bool sdPhiCut = (ljetfrcem < 0.4 && dphi < 0.15) && isdijet;
      bool lETCut = ljetfrcem < 0.1 && (ljetET > (50*ljetfrcem+20)) && (hdPhiCut || !isdijet);
      bool hETCut = ljetfrcem > 0.9 && (ljetET > (-50*ljetfrcem+80)) && (hdPhiCut || !isdijet);
      
      bool fullcut = bbCut || lETCut || hETCut;
      h1_forrat[2]->Fill(ljetET);
      
      if(isdijet && ljetET > 20 && subjetET > 10)
	{
	  h1_forrat[3]->Fill((ljetET-subjetET)/(ljetET+subjetET));
	  xJ[0]->Fill(subjetET/ljetET);
	}
      if(!fullcut)
	{
	  h1_zdist[1]->Fill(vtx[2]);
	  h1_forrat[0]->Fill(ljetET);
	  h1_cspec->Fill(ljetET);
	  for(int j=0; j<njet; ++j)
	    {
	      if(jet_e[j] > ljetET)
		{
		  subjetET = ljetET;
		  subjetph = ljetph;
		  ljetET = jet_e[j];
		  ljetph = jet_ph[j];
		  ljetfrcem = frcem[j];
		  ljeteta = jet_et[j];
		}
	      if(jet_e[j] >8)
		{
		  if(abs(jet_et[j]) < 0.7) continue;
		  hists2[3]->Fill(ljetfrcem,jet_e[j]);
		  h1_cspec->Fill(jet_e[j]);
		  
		}
	    }
	  if(isdijet)
	    {
	      h1_forrat[1]->Fill(ljetET);
	      if(ljetET > 20 && subjetET > 10)
		{
		  h1_forrat[4]->Fill((ljetET-subjetET)/(ljetET+subjetET));
		  xJ[1]->Fill(subjetET/ljetET);
		}
	    }
	}
      h1_zdist[0]->Fill(vtx[2]); 
    }


  h1_zdist[0]->Write();
  h1_ucspec->Scale(4e-5/2.8e6);
  h1_ucspec->Write();
  h1_cspec->Scale(4e-5/2.8e6);
  h1_cspec->Write();
  xJ[0]->Write();
  xJ[1]->Write();
  for(int i=0; i<5; ++i)  h1_forrat[i]->Write();

  jetfile->Write();
  return 0;
}
