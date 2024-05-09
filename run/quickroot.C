#include "TApplication.h"
#include "TROOT.h"
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
float get_eta(float eta)
{
  return (eta+1.1)*24/2.2;
}
float get_phi(float phi)
{
  return (phi+M_PI)*64/(2*M_PI);
}
int quickroot(string filebase="")
{
  gStyle->SetPalette(kOcean);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  float emcalen[25000];
  float emcalchi2[25000];
  float emcalt[25000];
  int sectorrt[3];
  int sectorem;
  int njet;
  float enrt[3][1536];
  int etart[3][1536];
  int phirt[3][1536];
  int etabin[25000];
  int phibin[25000];
  float jet_e[1000];
  float jet_et[1000];
  float jet_ph[1000];
  string filename=filebase;
  TChain* tree = new TChain("ttree");
  ifstream list;
  string line;
  list.open(filename,ifstream::in);
  if(!list)
    {
      cout << "nofile" << endl;
      exit(1);
    }
  while(getline(list,line))
    {
      try
	{
	  tree->Add(line.c_str());
	}
      catch(...)
	{
	  continue;
	}
    }
  TH1D* hist = new TH1D("hist","",200,-20,50);
  TH2D* event_display = new TH2D("event_display","",96,-0.5,95.5,256,-0.5,255.5);
  TH2D* event_diphcal = new TH2D("event_display_hc","",24,0,24,64,0,64);
  TH2D* event_disrt[3];
  TH2D* event_sum = new TH2D("event_sum","",24,0,24,64,0,64);
  for(int i=0; i<3; ++i)
    {
      event_disrt[i] = new TH2D(("event_display_rt"+to_string(i)).c_str(),"",24,0,24,64,0,64);
    }
  TH1D* jetfrac = new TH1D("jetfrac","",100,0,1);
  TH1D* jetE = new TH1D("jetE","",200,0,20);
  TH1D* hejet = new TH1D("hejet","",100,0,10);
  jetfrac->GetYaxis()->SetTitle("Counts");
  jetfrac->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
  jetE->GetYaxis()->SetTitle("Counts");
  jetE->GetXaxis()->SetTitle("E_{jet}");
  hejet->GetXaxis()->SetTitle("E_{jet, EMCal} / E_{jet, HCals}");
  hejet->GetYaxis()->SetTitle("Counts");
  float emetot;
  float seedD[1000];
  float ehjet[1000];
  tree->SetBranchAddress("ehjet",ehjet);
  tree->SetBranchAddress("seedD",seedD);
  tree->SetBranchAddress("njet",&njet);
  tree->SetBranchAddress("sectoroh",&sectorrt[2]);
  tree->SetBranchAddress("ohcaletabin",etart[2]);
  tree->SetBranchAddress("ohcalphibin",phirt[2]);
  tree->SetBranchAddress("ohcalen",enrt[2]);
  tree->SetBranchAddress("sectorih",&sectorrt[1]);
  tree->SetBranchAddress("sector_rtem",&sectorrt[0]);
  tree->SetBranchAddress("rtemen",enrt[0]);
  tree->SetBranchAddress("rtemet",etart[0]);
  tree->SetBranchAddress("rtemph",phirt[0]);
  tree->SetBranchAddress("jet_e",jet_e);
  tree->SetBranchAddress("jet_et",jet_et);
  tree->SetBranchAddress("jet_ph",jet_ph);
  tree->SetBranchAddress("ihcaletabin",etart[1]);
  tree->SetBranchAddress("ihcalphibin",phirt[1]);
  tree->SetBranchAddress("ihcalen",enrt[1]);
  //tree->SetBranchAddress("emcalen",emcalen);
  //tree->SetBranchAddress("emcalt",emcalt);
  //tree->SetBranchAddress("emcalchi2",emcalchi2);
  //tree->SetBranchAddress("sectorem",&sectorem);
  //tree->SetBranchAddress("emcaletabin",etabin);
  //tree->SetBranchAddress("emcalphibin",phibin);
  //tree->SetBranchAddress("emetot",&emetot);
  int cancount = 0;
  TCanvas* c = new TCanvas("","",500,1000);
  TCanvas* d = new TCanvas("","",1000,1000);
  c->Divide(2,2);
  event_disrt[0]->GetXaxis()->SetTitle("EMCal #eta Bin");
  event_disrt[0]->GetYaxis()->SetTitle("EMCal #phi Bin");
  event_disrt[1]->GetXaxis()->SetTitle("IHCal #eta Bin");
  event_disrt[1]->GetYaxis()->SetTitle("IHCal #phi Bin");
  event_disrt[2]->GetXaxis()->SetTitle("OHCal #eta Bin");
  event_disrt[2]->GetYaxis()->SetTitle("OHCal #phi Bin");
  for(int i=0; i<3; ++i)
    {
      event_disrt[i]->GetZaxis()->SetTitleOffset(2);
      event_disrt[i]->GetZaxis()->SetTitle("Tower Energy [GeV]");
    }
  event_sum->GetXaxis()->SetTitle("Tower Sum #eta Bin");
  event_sum->GetYaxis()->SetTitle("Tower Sum #phi Bin");
  event_sum->GetZaxis()->SetTitle("Tower Energy [GeV]");
  event_sum->GetZaxis()->SetTitleOffset(2);
  int ncircle = 64;
  int highejet = 0;
  for(int i=0; i<tree->GetEntries(); ++i)
    {
      //if(i % 1000 == 0) cout << i << endl;
      tree->GetEntry(i);
      int countedjets = 0;
      if(!njet) continue;
      //event_display->Reset();
      for(int j=0; j<njet; ++j)
	{
	  if(ehjet[j] > 10) continue;
	  if(seedD[j] > 0.65) continue;
	  ++countedjets;
	}
      if(!countedjets) continue;
      event_sum->Reset();
      for(int j=0; j<3; ++j)
	{
	  event_disrt[j]->Reset();
	  for(int k=0; k<sectorrt[j]; ++k)
	    {
	      event_disrt[j]->Fill(etart[j][k],phirt[j][k],enrt[j][k]);
	      event_sum->Fill(etart[j][k],phirt[j][k],enrt[j][k]);
	      //if(enrt[j][k] > 0.1) cout << j << " " << k << " " << etart[j][k] << " " << phirt[j][k] << " " << enrt[j][k] << endl;
	    }
	  c->cd(j+1);
	  gPad->SetFrameFillColor(kBlack);
	  gPad->SetLogz();
	  gPad->SetRightMargin(0.2);
	  event_disrt[j]->Draw("COLZ");
	  TMarker* jets[1000];
	  for(int k=0; k<njet; ++k)
	    {
	      if(ehjet[k] > 10) continue;
	      if(seedD[k] > 0.65) continue;
	      if(jet_e[k] > 8) highejet = 1;
	      jets[k] = new TMarker(get_eta(jet_et[k]),get_phi(jet_ph[k]),43);
	      jets[k]->SetMarkerSize(jet_e[k]/3);
	      jets[k]->SetMarkerColor(kRed);
	      //jets[k]->Draw();
	      for(int l=0; l<ncircle; ++l)
		{
		  TMarker* circlemarker = new TMarker(get_eta(jet_et[k]+0.4*cos(2*l*M_PI/ncircle)),get_phi(jet_ph[k]+0.4*sin(2*l*M_PI/ncircle)),20);
		  circlemarker->SetMarkerSize(0.2);
		  circlemarker->SetMarkerColor(kRed);
		  circlemarker->Draw();
		}
	    }
	}
      c->cd(4);
      gPad->SetLogz();
      gPad->SetRightMargin(0.2);
      gPad->SetFrameFillColor(kBlack);
      event_sum->Draw("COLZ");
      for(int j=0; j<24; ++j)
	{
	  for(int k=0; k<64; ++k)
	    {
	      float val = event_sum->GetBinContent(j+1,k+1);
	      //if(val > 0.5) cout << j << " " << k << " " << val << endl;
	    }
	}
      TMarker* jets[1000];
      for(int k=0; k<njet; ++k)
	{
	  jetfrac->Fill(seedD[k]);
	  if(seedD[k] > 0.65) continue;
	  hejet->Fill(ehjet[njet]);
	  jetE->Fill(jet_e[k]);
	  jets[k] = new TMarker(get_eta(jet_et[k]),get_phi(jet_ph[k]),20);
	  jets[k]->SetMarkerSize(jet_e[k]/3);
	  jets[k]->SetMarkerColor(kRed);
	  //jets[k]->Draw();
	  for(int l=0; l<ncircle; ++l)
	    {
	      TMarker* circlemarker = new TMarker(get_eta(jet_et[k]+0.4*cos(2*l*M_PI/ncircle)),get_phi(jet_ph[k]+0.4*sin(2*l*M_PI/ncircle)),20);
	      circlemarker->SetMarkerSize(0.2);
	      circlemarker->SetMarkerColor(kRed);
	      circlemarker->Draw();
	    }
	}
      if(highejet)
	{
	  c->SaveAs(("./output/img/candidate_"+filebase+"_"+to_string(cancount)+".pdf").c_str());
	  highejet = 0;
	}
      ++cancount;
    }
  d->cd();
  d->SetLogy();
  jetfrac->Draw();
  d->SaveAs(("output/img/"+filebase+"_jetfrac.pdf").c_str());
  jetE->Draw();
  d->SaveAs(("output/img/"+filebase+"_jetE.pdf").c_str());
  hejet->Draw();
  d->SaveAs(("output/img/"+filebase+"_hejet.pdf").c_str());
  return 0;
}
