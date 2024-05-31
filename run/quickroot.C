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
float get_eta(float eta)
{
  return (eta+1.1)*24/2.2;
}
float get_phi(float phi)
{
  return (phi+M_PI)*64/(2*M_PI);
}
float bintoeta_hc(int etabin)
{
  return (2.2*etabin)/24 - 1.1;
}
float bintophi_hc(int phibin)
{
  return (2*M_PI*phibin)/64;
}
float bintoeta_em(int etabin)
{
  return (2.2*etabin)/96 - 1.1;
}
float bintophi_em(int phibin)
{
  return (2*M_PI*phibin)/256;
}
int save_etaphiE(int etabins[], int phibins[], float energies[], int nsec, string filename, int em)
{
  ofstream thefile;
  thefile.open(filename);
  for(int i=0; i<nsec; ++i)
    {
      thefile << (i>0?",":"") <<"{\"eta\": " << (em?bintoeta_em(etabins[i]):bintoeta_hc(etabins[i])) << ", \"phi\": " << (em?bintophi_em(phibins[i]):bintophi_hc(phibins[i])) << ", \"e\": " << energies[i] << ", \"event\": 0}" << endl;
    }
  thefile.close();
  return 0;
}

int quickroot(string filebase="")
{
  gStyle->SetPalette(kOcean);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  //gStyle->SetOptTitle(0);
  float emcalen[2][25000];
  float emcalchi2[2][25000];
  float emcalt[2][25000];
  int sectorrt[2][3];
  int sectorem[2];
  int njet[2];
  float enrt[2][3][1536];
  int etart[2][3][1536];
  int phirt[2][3][1536];
  int etabin[2][25000];
  int phibin[2][25000];
  float jet_e[2][1000];
  float jet_et[2][1000];
  float jet_ph[2][1000];
  string filename=filebase;
  TChain* tree[2];
  TChain* tree2[2];
  ifstream list[2];
  string line;
  list[0].open("sim.list",ifstream::in);
  list[1].open(filename,ifstream::in);

  for(int h=0; h < 2; ++h)
    {
      tree[h] = new TChain("ttree");
      tree2[h] = new TChain("ttree2");
      if(!list[h])
	{
	  cout << "nosimlist" << endl;
	  exit(1);
	}
      while(getline(list[h],line))
	{
	  try
	    {
	      tree[h]->Add(line.c_str());
	      tree2[h]->Add(line.c_str());
	    }
	  catch(...)
	    {
	      continue;
	    }
	}
    }
  long long unsigned int nevents[2] = {0};
  int evtct[2];
  int mbevt[2];
  long long unsigned int nmb[2] = {0};
  for(int h=0; h<2; ++h)
    {
      tree2[h]->SetBranchAddress("_evtct",&evtct[h]);
      tree2[h]->SetBranchAddress("mbevt",&mbevt[h]);
      for(int i=0; i<tree2[h]->GetEntries(); ++i)
	{
	  tree2[h]->GetEntry(i);
	  nevents[h] += evtct[h];
	  nmb[h] += mbevt[h];
	  //cout << mbevt << endl;
	}
      cout << "Total " << nevents[h] << " events." << endl;
      if(h==0) nmb[h] = nevents[h];
      cout << "with " << nmb[h] << " minbias events." << endl;
    }

  const int ncol = 9;
  double red[ncol] = {1, 1, 1, 1, 1, 1, 1, 1};
  double grn[ncol] = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double blu[ncol] = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double stp[ncol] = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
  TColor::CreateGradientColorTable(ncol, stp, red, grn, blu, ncol);
  TH1D* h1_rej[2][4];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  h1_rej[h][i] = new TH1D(("h1_rej"+to_string(h)+"_"+to_string(i)).c_str(),"",150,0,15);
	  h1_rej[h][i]->GetYaxis()->SetRangeUser(0.5,200000);
	}
    }
  TH2D* event_display = new TH2D("event_display","",96,-0.5,95.5,256,-0.5,255.5);
  TH2D* event_diphcal = new TH2D("event_display_hc","",24,0,24,64,0,64);
  TH1D* h1_dphi[2][4];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  h1_dphi[h][i] = new TH1D(("h1_dphi"+to_string(h)+"_"+to_string(i)).c_str(),"",32,0,2*M_PI);
	}
    }
  TH2D* h2_dphi_deta[2];
  TH1D* h1_AJ[2][4];
    for(int h=0; h<2; ++h)
    {
      h2_dphi_deta[h] = new TH2D(("h2_dphi_deta"+to_string(h)).c_str(),"",20,-2.2,2.2,32,0,2*M_PI);
      for(int i=0; i<4; ++i)
	{
	  h1_AJ[h][i] = new TH1D(("h1_AJ"+to_string(h)+"_"+to_string(i)).c_str(),"",100,0,1);
	}
    }
  TH2D* event_disrt[3];
  TH2D* event_sum = new TH2D("event_sum","Calorimeter Sum",24,0,24,64,0,64);
  TH1D* jetfrac[2][3];
  TH1D* jetE[2][4];
  TH1D* hejet[2];
  for(int h=0; h<2; ++h)
    {
      hejet[h] = new TH1D(("hejet"+to_string(h)).c_str(),"",100,0,10);
      hejet[h]->SetMarkerStyle(43);
      hejet[h]->SetMarkerSize(3);
      hejet[h]->SetMarkerColor(kMagenta+3);
      for(int i=0; i<3; ++i)
	{
	  string calo = "EMCal";
	  if(i==1) calo = "IHCal";
	  if(i==2) calo = "OHCal";
	  jetE[h][i] = new TH1D(("jetE"+to_string(h)+"_"+to_string(i)).c_str(),"",200,0,20);
	  jetfrac[h][i] = new TH1D(("jetfrac"+to_string(h)+"_"+to_string(i)).c_str(),"",100,0,1);
	  event_disrt[i] = new TH2D(("event_display_rt"+to_string(h)+"_"+to_string(i)).c_str(),calo.c_str(),24,0,24,64,0,64);
	  jetfrac[h][i]->SetMarkerStyle(43);
	  jetfrac[h][i]->SetMarkerColor(kMagenta+3);
	  jetfrac[h][i]->SetMarkerSize(3);
	  jetfrac[h][i]->GetYaxis()->SetLabelSize(0.02);
	  jetfrac[h][i]->GetYaxis()->SetTitleSize(0.03);
	  jetfrac[h][i]->GetYaxis()->SetTitleOffset(1.7);
	  jetfrac[h][i]->GetXaxis()->SetLabelSize(0.02);
	  jetfrac[h][i]->GetXaxis()->SetTitleSize(0.03);
	  jetE[h][i]->SetMarkerStyle(43);
	  jetE[h][i]->SetMarkerSize(3);
	  jetE[h][i]->SetMarkerColor(kMagenta+3);
	  jetE[h][i]->SetMarkerStyle(43);
	  jetE[h][i]->SetMarkerColor(kMagenta+3);
	  jetE[h][i]->SetMarkerSize(3);
	  jetE[h][i]->GetYaxis()->SetLabelSize(0.02);
	  jetE[h][i]->GetYaxis()->SetTitleSize(0.03);
	  jetE[h][i]->GetYaxis()->SetTitleOffset(1.7);
	  jetE[h][i]->GetXaxis()->SetLabelSize(0.02);
	  jetE[h][i]->GetXaxis()->SetTitleSize(0.03);
	}
  
      jetE[h][3] = new TH1D(("jetE"+to_string(3)).c_str(),"",200,0,20);
      jetE[h][3]->SetMarkerColor(kMagenta+3);
      jetE[h][3]->SetMarkerSize(3);
      jetE[h][3]->SetMarkerStyle(43);
      jetE[h][3]->SetMarkerColor(kMagenta+3);
      jetE[h][3]->GetYaxis()->SetLabelSize(0.02);
      jetE[h][3]->GetYaxis()->SetTitleSize(0.03);
      jetE[h][3]->GetYaxis()->SetTitleOffset(1.7);
      jetE[h][3]->GetXaxis()->SetLabelSize(0.02);
      jetE[h][3]->GetXaxis()->SetTitleSize(0.03);
      jetfrac[h][0]->GetYaxis()->SetTitle("Event Normalized Counts 4 < E_{jet} < 7.5");
      jetfrac[h][0]->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
      jetfrac[h][1]->GetYaxis()->SetTitle("Event Normalized Counts 7.5 < E_{jet} < 10");
      jetfrac[h][1]->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
      jetfrac[h][2]->GetYaxis()->SetTitle("Event Normalized Counts 10 < E_{jet}");
      jetfrac[h][2]->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
      jetE[h][0]->GetYaxis()->SetTitle("Event Normalized Counts");
      jetE[h][0]->GetXaxis()->SetTitle("E_{jet} No Cuts");
      jetE[h][1]->GetYaxis()->SetTitle("Event Normalized Counts");
      jetE[h][1]->GetXaxis()->SetTitle("E_{jet}  with E_{EMCal} < 10E_{HCal}");
      jetE[h][2]->GetYaxis()->SetTitle("Event Normalized Counts");
      jetE[h][2]->GetXaxis()->SetTitle("E_{jet}  with E_{lead tower} < 0.65E_{jet}");
      jetE[h][3]->GetYaxis()->SetTitle("Event Normalized Counts");
      jetE[h][3]->GetXaxis()->SetTitle("E_{jet}  with E_{EMCal} < 10E_{HCal} \& E_{lead tower} < 0.65E_{jet}");
      hejet[h]->GetXaxis()->SetTitle("E_{jet, EMCal} / E_{jet, HCals}");
      hejet[h]->GetYaxis()->SetTitle("Event Normalized Counts");
      hejet[h]->SetMarkerStyle(43);
      hejet[h]->SetMarkerColor(kMagenta+3);
      hejet[h]->SetMarkerSize(3);
      hejet[h]->GetYaxis()->SetLabelSize(0.02);
      hejet[h]->GetYaxis()->SetTitleSize(0.03);
      hejet[h]->GetYaxis()->SetTitleOffset(1.7);
      hejet[h]->GetXaxis()->SetLabelSize(0.02);
      hejet[h]->GetXaxis()->SetTitleSize(0.03);
    }
  float emetot[2];
  float seedD[2][1000];
  float ehjet[2][1000];
  long unsigned int trigvec;
  tree[1]->SetBranchAddress("triggervec",&trigvec);
  for(int h=0; h<2; ++h)
    {
      tree[h]->SetBranchAddress("ehjet",ehjet[h]);
      tree[h]->SetBranchAddress("seedD",seedD[h]);
      tree[h]->SetBranchAddress("njet",&njet[h]);
      tree[h]->SetBranchAddress("sectoroh",&sectorrt[h][2]);
      tree[h]->SetBranchAddress("ohcaletabin",etart[h][2]);
      tree[h]->SetBranchAddress("ohcalphibin",phirt[h][2]);
      tree[h]->SetBranchAddress("ohcalen",enrt[h][2]);
      tree[h]->SetBranchAddress("sectorih",&sectorrt[h][1]);
      tree[h]->SetBranchAddress("sector_rtem",&sectorrt[h][0]);
      tree[h]->SetBranchAddress("rtemen",enrt[h][0]);
      tree[h]->SetBranchAddress("rtemet",etart[h][0]);
      tree[h]->SetBranchAddress("rtemph",phirt[h][0]);
      tree[h]->SetBranchAddress("jet_e",jet_e[h]);
      tree[h]->SetBranchAddress("jet_et",jet_et[h]);
      tree[h]->SetBranchAddress("jet_ph",jet_ph[h]);
      tree[h]->SetBranchAddress("ihcaletabin",etart[h][1]);
      tree[h]->SetBranchAddress("ihcalphibin",phirt[h][1]);
      tree[h]->SetBranchAddress("ihcalen",enrt[h][1]);
      //tree->SetBranchAddress("emcalen",emcalen);
      //tree->SetBranchAddress("emcalt",emcalt);
      //tree->SetBranchAddress("emcalchi2",emcalchi2);
      //tree->SetBranchAddress("sectorem",&sectorem);
      //tree->SetBranchAddress("emcaletabin",etabin);
      //tree->SetBranchAddress("emcalphibin",phibin);
      //tree->SetBranchAddress("emetot",&emetot);
    }
  int cancount = 0;
  TCanvas* c = new TCanvas("","",800,500);
  TCanvas* d = new TCanvas("","",1000,1000);
  c->Divide(4,1);
  event_disrt[0]->GetXaxis()->SetTitle("EMCal #eta Bin");
  event_disrt[0]->GetYaxis()->SetTitle("EMCal #phi Bin");
  event_disrt[1]->GetXaxis()->SetTitle("IHCal #eta Bin");
  event_disrt[1]->GetYaxis()->SetTitle("IHCal #phi Bin");
  event_disrt[2]->GetXaxis()->SetTitle("OHCal #eta Bin");
  event_disrt[2]->GetYaxis()->SetTitle("OHCal #phi Bin");
  for(int i=0; i<3; ++i)
    {
      event_disrt[i]->GetZaxis()->SetTitleOffset(1.5);
      event_disrt[i]->GetZaxis()->SetTitle("Tower Energy [GeV]");
      event_disrt[i]->GetZaxis()->SetRangeUser(0.1,5);
    }
  event_sum->GetXaxis()->SetTitle("Tower Sum #eta Bin");
  event_sum->GetYaxis()->SetTitle("Tower Sum #phi Bin");
  event_sum->GetZaxis()->SetTitle("Tower Energy [GeV]");
  event_sum->GetZaxis()->SetTitleOffset(1.5);
  event_sum->GetZaxis()->SetRangeUser(0.1,5);
  int ncircle = 64;
  int highejet = 0;
  float fakejets = 0;
  int dispcount = 1;
  int shighe = 0;
  int sshighe = 0;
  int maxeh = 5;
  for(int h=0; h<2; ++h)
    {
      cancount = 0;
  for(int i=0; i<tree[h]->GetEntries(); ++i)
    {
      highejet = 0;
      shighe = 0;
      sshighe = 0;
      int passcut = 1;
      //if(i % 1000 == 0) cout << i << endl;
      tree[h]->GetEntry(i);
      if(!njet[h]) continue;
      if(njet[h] > 4) continue;
      //event_display->Reset();
      //int breakvar = 0;
      float minjet[4] = {9999,9999,9999,9999};
      //float minjetseedD[4] = {0};
      //float minjeteh[4] = {0};
      //if(trigvec >> 10 & 1)
      //{
      for(int j=0; j<njet[h]; ++j)
	{
	  for(int k=0; k<4; ++k)
	    {
	      if(jet_e[h][j] < minjet[k])
		{
		  if(k==0)
		    {
		      minjet[k] = jet_e[h][j];
		      //minjetseedD[k] = seedD[j];
		      //minjeteh[k] = ehjet[j];
		    }
		  else if(k==1 && ehjet[h][j] < maxeh)
		    {
		      minjet[k] = jet_e[h][j];
		      //minjetseedD[k] = seedD[j];
		      //minjeteh[k] = ehjet[j];
		    }
		  else if(k==2 && seedD[h][j] < 0.65)
		    {
		      minjet[k] = jet_e[h][j];
		      //minjetseedD[k] = seedD[j];
		      //minjeteh[k] = ehjet[j];
		    }
		  else if(k==3 && seedD[h][j] <0.65 && ehjet[h][j] < maxeh)
		    {
		      minjet[k] = jet_e[h][j];
		      //minjetseedD[k] = seedD[j];
		      //minjeteh[k] = ehjet[j];
		    }
		}
	    }
	}
      //}
      for(int j=0; j<4; ++j)
	{
	  if(minjet[j] > 9000)
	    {
	      minjet[j] = 0;
	    }
	}
      for(int j=0; j<4; ++j)
	{
	  int it = 40;
	  while(it < 10*minjet[j])
	    {
	      h1_rej[h][j]->Fill(h1_rej[h][j]->GetBinCenter(it+1));
	      ++it;
	    }
	}
      //if(breakvar) continue;
      event_sum->Reset();
      for(int j=0; j<3; ++j)
	{
	  event_disrt[j]->Reset();
	  for(int k=0; k<sectorrt[h][j]; ++k)
	    {
	      event_disrt[j]->Fill(etart[h][j][k],phirt[h][j][k],enrt[h][j][k]);
	      event_sum->Fill(etart[h][j][k],phirt[h][j][k],enrt[h][j][k]);
	      //if(enrt[j][k] > 0.1) cout << j << " " << k << " " << etart[j][k] << " " << phirt[j][k] << " " << enrt[j][k] << endl;
	    }

	  //gPad->SetFrameFillColor(kBlack);
	  gPad->SetLogz();
	  gPad->SetRightMargin(0.2);
	  //gPad->SetTopMargin(0);
	  //gPad->SetBottomMargin(0);
	  //cout << (12*((dispcount%18)/3)+2*(dispcount%3)+(j==2?7:j+1)) << endl;
	  //c->cd(12*((dispcount%18)/3)+2*(dispcount%3)+(j==2?7:j+1));
	  c->cd(j+1);
	  event_disrt[j]->SetTitle(("Event "+to_string(i)).c_str());
	  event_disrt[j]->Draw("COLZ");
	  TMarker* jets[1000];
	  for(int k=0; k<njet[h]; ++k)
	    {
	      /*
	      if(get_phi(jet_ph[k]) > 29.5 && get_phi(jet_ph[k]) < 37.5)
		{
		  continue;
		}
	      */
	      if(ehjet[h][k] > maxeh || ehjet[h][k] < 0) continue;
	      if(seedD[h][k] > 0.65) continue;
	      /*
	      if(jet_e[h][k] > 10)
		{
		  highejet = 1;
		}
	      */
	      if(jet_e[h][k] > 20) shighe = 1;
	      if(jet_e[h][k] > 25)
		{
		  sshighe = 1;
		  save_etaphiE(etart[h][0], phirt[h][0], enrt[h][0], sectorrt[h][0], "output/json/"+filebase+"_"+to_string(h)+"_"+to_string(i)+"_em.txt", 0);
		  save_etaphiE(etart[h][1], phirt[h][1], enrt[h][1], sectorrt[h][1], "output/json/"+filebase+"_"+to_string(h)+"_"+to_string(i)+"_ih.txt", 0);
		  save_etaphiE(etart[h][2], phirt[h][2], enrt[h][2], sectorrt[h][2], "output/json/"+filebase+"_"+to_string(h)+"_"+to_string(i)+"_oh.txt", 0);
		}
	      jets[k] = new TMarker(get_eta(jet_et[h][k]),get_phi(jet_ph[h][k]),43);
	      jets[k]->SetMarkerSize(jet_e[h][k]/3);
	      jets[k]->SetMarkerColor(kRed);
	      //jets[k]->Draw();
	      for(int l=0; l<ncircle; ++l)
		{
		  float eta = get_eta(jet_et[h][k]+0.4*cos(2*l*M_PI/ncircle));
		  float phi = get_phi(jet_ph[h][k]+0.4*sin(2*l*M_PI/ncircle));
		  if(eta > 24 || eta < 0) continue;
		  if(phi > 64) phi -= 64;
		  if(phi < 0) phi += 64;
		  TMarker* circlemarker = new TMarker(eta,phi,20);
		  circlemarker->SetMarkerSize(0.3);
		  circlemarker->SetMarkerColor(kBlue);
		  circlemarker->Draw();
		}
	    }
	}
      gPad->SetLogz();
      gPad->SetRightMargin(0.2);
      //gPad->SetTopMargin(0);
      //gPad->SetBottomMargin(0);
      //gPad->SetFrameFillColor(kBlack);
      //c->cd(12*((dispcount%18)/3)+2*(dispcount%3)+8);
      c->cd(4);
      event_sum->SetTitle(("Event "+to_string(i)).c_str());
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
      int dodraw = 0;
      for(int k=0; k<njet[h]; ++k)
	{
	  

	  /*
	  if(get_phi(jet_ph[k]) > 29.5 && get_phi(jet_ph[k]) < 37.5)
	    {
	      continue;
	    }
	  */
	  if(jet_e[h][k] < 7.5)
	    {
	      jetfrac[h][0]->Fill(seedD[h][k]);
	    }
	  else if(jet_e[h][k] < 10)
	    {
	      jetfrac[h][1]->Fill(seedD[h][k]);
	    }
	  else
	    {
	      jetfrac[h][2]->Fill(seedD[h][k]);
	    }
	  //if(seedD[k] > 0.65) continue;
	  hejet[h]->Fill(ehjet[h][k]);
	  jetE[h][0]->Fill(jet_e[h][k]);
	  if(ehjet[h][k] < maxeh && ehjet[h][k] > 0) jetE[h][1]->Fill(jet_e[h][k]);
	  if(seedD[h][k] < 0.65) jetE[h][2]->Fill(jet_e[h][k]);
	  if(seedD[h][k] < 0.65 && ehjet[h][k] < maxeh && ehjet[h][k] > 0)
	    {
	      jetE[h][3]->Fill(jet_e[h][k]);
	    }
	  
	  for(int l = 0; l<njet[h]; ++l)
	    {
	      if(l==k) continue;
	      if(jet_e[h][k] < jet_e[h][l]) continue;
	      if(jet_e[h][k] > 5 && jet_e[h][l] > 5)
		{
		  float ET1 = jet_e[h][k]/cosh(jet_et[h][k]);
		  float ET2 = jet_e[h][l]/cosh(jet_et[h][l]);
		  float AJ = (ET1-ET2)/(ET1+ET2);
		  h1_AJ[h][0]->Fill(AJ);
		  if(seedD[h][k] < 0.65) h1_AJ[h][2]->Fill(AJ);
		  if(ehjet[h][k] < maxeh) h1_AJ[h][1]->Fill(AJ);
		  if(ehjet[h][k] < maxeh && seedD[h][k] < 0.65) h1_AJ[h][3]->Fill(AJ);
		  float dphi = jet_ph[h][k] - jet_ph[h][l];
		  float deta = jet_et[h][k] - jet_et[h][l];
		  if(sqrt(dphi*dphi+deta*deta) < 0.4 && seedD[h][k] < 0.65 && seedD[h][l] < 0.65 && ehjet[h][k] < maxeh && ehjet[h][k] > 0 && ehjet[h][l] > 0 && ehjet[h][l] < maxeh) dodraw = 1;
		  if(dphi < 0) dphi+=2*M_PI;
		  //cout << dphi << endl;
		  h1_dphi[h][0]->Fill(dphi);
		  if(seedD[h][k] < 0.65) h1_dphi[h][2]->Fill(dphi);
		  if(ehjet[h][k] < maxeh) h1_dphi[h][1]->Fill(dphi);
		  if(ehjet[h][k] < maxeh && seedD[h][k] < 0.65) h1_dphi[h][3]->Fill(dphi);
		  h2_dphi_deta[h]->Fill(deta,dphi);
		}
	    }
	  if(seedD[h][k] > 0.65) continue;
	  if(ehjet[h][k] > maxeh || ehjet[h][k] < 0) continue;
	  jets[k] = new TMarker(get_eta(jet_et[h][k]),get_phi(jet_ph[h][k]),20);
	  jets[k]->SetMarkerSize(jet_e[h][k]/3);
	  jets[k]->SetMarkerColor(kRed);
	  //jets[k]->Draw();
	  for(int l=0; l<ncircle; ++l)
	    {
	      float eta = get_eta(jet_et[h][k]+0.4*cos(2*l*M_PI/ncircle));
	      float phi = get_phi(jet_ph[h][k]+0.4*sin(2*l*M_PI/ncircle));
	      if(eta > 24 || eta < 0) continue;
	      if(phi > 64) phi -= 64;
	      if(phi < 0) phi += 64;
	      TMarker* circlemarker = new TMarker(eta,phi,20);
	      circlemarker->SetMarkerSize(0.3);
	      circlemarker->SetMarkerColor(kBlue);
	      circlemarker->Draw();
	    }
	  std::stringstream e_stream;
	  e_stream << std::fixed << std::setprecision(2) << jet_e[h][k];
	  std::string e_string = e_stream.str();
	  drawText((e_string+" GeV").c_str(),get_eta(jet_et[h][k]),get_phi(jet_ph[h][k]),(get_eta(jet_et[h][k])>15?1:0),kBlack,0.08,42,false);
	  //d->cd();
	  //event_sum->Draw("LEGO");
	}
      //if(dispcount % 18 == 17 && highejet)
      if(dodraw)
	{
	  //cout << "Saving display " << dispcount/18 << endl;
	  //c->SaveAs(("./output/img/candidate_"+filebase+"_"+to_string(dispcount/18)+".pdf").c_str());
	  c->SaveAs(("./output/img/candidate_"+filebase+"_"+to_string(h)+"_"+to_string(cancount)+".png").c_str());
	  //cout << "Saved!" << endl;
	  //c->Clear("D");
	}
      ++cancount;
      if(shighe)
	{
	  c->SaveAs(("./output/hmg/candidate_"+filebase+"_superhighE_"+to_string(h)+"_"+to_string(cancount)+".png").c_str());
	  //d->SaveAs(("./output/hmg/candidate_"+filebase+"_superhighElego_"+to_string(cancount)+".png").c_str());
	}
      if(sshighe) 
	{
	  c->SaveAs(("./output/smg/candidate_"+filebase+"_supersuperhighE_"+to_string(h)+"_"+to_string(cancount)+".png").c_str());
	  //d->SaveAs(("./output/smg/candidate_"+filebase+"_supersuperhighElego_"+to_string(cancount)+".png").c_str());
	}
      if(highejet)
	{
	  //cout << dispcount << endl;
	  ++dispcount;
	}
    }
    }
  //if(dispcount % 18 != 0)
  //  {
  //    c->SaveAs(("./output/img/candidate_"+filebase+"_last.png").c_str());
  //  }
  d->cd();
  float erejf[9] = {14.24,42.76,111.95,281.76,681.85,1466.6,3033.03,6272.51,11086.03};
  float x[9] = {4,5,6,7,8,9,10,11,12};

  for(int h=0; h<2; ++h)
    {
      gPad->SetLogy(0);
      for(int i=0; i<3; ++i)
	{
	  //jetfrac[h][i]->Scale(1./(nevents[h]));
	}
      //h1_dphi->SetMarkerColor(kMagenta+3);
      for(int i=0; i<4; ++i)
	{
	  h1_dphi[h][i]->GetYaxis()->SetLabelSize(0.02);
	  h1_dphi[h][i]->GetYaxis()->SetTitleSize(0.03);
	  h1_dphi[h][i]->GetYaxis()->SetTitleOffset(1.7);
	  h1_dphi[h][i]->GetXaxis()->SetLabelSize(0.02);
	  h1_dphi[h][i]->GetXaxis()->SetTitleSize(0.03);
	  //h1_dphi[h][i]->Scale(1./nevents[h]);
	  h1_dphi[h][i]->GetYaxis()->SetTitle("Event Normalized Counts");
	  h1_dphi[h][i]->GetXaxis()->SetTitle("#Delta#phi [rad]");
	  h1_dphi[h][i]->SetFillColorAlpha(kMagenta,0.4);
	  h1_dphi[h][i]->Draw("HIST");
	  drawText("Leading Jet E > 10 GeV", 0.5, 0.85, 0, kBlack, 0.02);
	  drawText("Subleading Jet E > 5 GeV", 0.5, 0.825, 0, kBlack, 0.02);
	  if(i>1) drawText("E_{lead tower} < 0.65 E_{jet}", 0.5, 0.8, 0, kBlack, 0.02);
	  if(i==1 || i==3) drawText(("E_{EM} < "+to_string(maxeh)+" E_{HAD}").c_str(), 0.5, (i==1?0.8:0.775), 0, kBlack, 0.02);
	  d->SaveAs(("output/gmg/"+filebase+"_"+to_string(h)+"_"+to_string(cancount)+"_"+to_string(i)+"_dphi_1d.png").c_str());
	  
	  h1_AJ[h][i]->GetYaxis()->SetLabelSize(0.02);
	  h1_AJ[h][i]->GetYaxis()->SetTitleSize(0.03);
	  h1_AJ[h][i]->GetYaxis()->SetTitleOffset(1.7);
	  h1_AJ[h][i]->GetXaxis()->SetLabelSize(0.02);
	  h1_AJ[h][i]->GetXaxis()->SetTitleSize(0.03);
	  //h1_AJ[h][i]->Scale(1./nevents[h]);
	  h1_AJ[h][i]->GetYaxis()->SetTitle("Event Normalized Counts");
	  h1_AJ[h][i]->GetXaxis()->SetTitle("A_{J}");
	  h1_AJ[h][i]->SetFillColorAlpha(kMagenta,0.4);
	  h1_AJ[h][i]->Draw("HIST");
	  drawText("Leading Jet E > 10 GeV", 0.5, 0.85, 0, kBlack, 0.02);
	  drawText("Subleading Jet E > 5 GeV", 0.5, 0.825, 0, kBlack, 0.02);
	  if(i>1) drawText("E_{lead tower} < 0.65 E_{jet}", 0.5, 0.8, 0, kBlack, 0.02);
	  if(i==1 || i==3) drawText(("E_{EM} < "+to_string(maxeh)+" E_{HAD}").c_str(), 0.5, (i==1?0.8:0.775), 0, kBlack, 0.02);
	  d->SaveAs(("output/gmg/"+filebase+"_"+to_string(h)+"_"+to_string(cancount)+"_"+to_string(i)+"_AJ_1d.png").c_str());
	}
      //h2_dphi_deta[h]->SetMarkerColor(kMagenta+3);
      h2_dphi_deta[h]->GetYaxis()->SetLabelSize(0.02);
      h2_dphi_deta[h]->GetYaxis()->SetTitleSize(0.03);
      h2_dphi_deta[h]->GetYaxis()->SetTitleOffset(1.7);
      h2_dphi_deta[h]->GetXaxis()->SetLabelSize(0.02);
      h2_dphi_deta[h]->GetXaxis()->SetTitleSize(0.03);
      //h2_dphi_deta[h]->Scale(1./nevents[h]);
      h2_dphi_deta[h]->GetYaxis()->SetTitle("#Delta#phi [rad]");
      h2_dphi_deta[h]->GetXaxis()->SetTitle("#Delta#eta [rad]");
      h2_dphi_deta[h]->GetZaxis()->SetTitle("Event Normalized Counts");
      h2_dphi_deta[h]->Draw("COLZ");
      d->SaveAs(("output/gmg/"+filebase+"_"+to_string(h)+"_"+to_string(cancount)+"_dphi_deta_2d.png").c_str());
      jetfrac[h][0]->Draw();
      d->SaveAs(("output/gmg/"+filebase+"_"+to_string(h)+"_"+to_string(cancount)+"_jetfrac0.png").c_str());
      jetfrac[h][1]->Draw();
      d->SaveAs(("output/gmg/"+filebase+"_"+to_string(h)+"_"+to_string(cancount)+"_jetfrac1.png").c_str());
      jetfrac[h][2]->Draw();
      d->SaveAs(("output/gmg/"+filebase+"_"+to_string(h)+"_"+to_string(cancount)+"_jetfrac2.png").c_str());
    }
  d->SetLogy();
  //float goodfrac = 1 - fakejets/tree->GetEntries();
  TGraph* erejg = new TGraph(9,x,erejf);
  for(int i=0; i<4; ++i)
    {
      //jetE[h][i]->Scale(1./(nevents[h]));
      jetE[0][i]->SetMarkerColor(kGreen);
      for(int h=0; h<2; ++h)
	{
	  jetE[h][i]->Draw();
	}
      d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_jetE"+to_string(i)+".png").c_str());
      for(int h=0; h<2; ++h)
	{
	  for(int j=0; j<160; ++j)
	    {
	      if(h1_rej[h][i]->GetBinContent(j+1) != 0) h1_rej[h][i]->SetBinContent(j+1,nmb[h]/h1_rej[h][i]->GetBinContent(j+1));
	    }
	  //      h1_rej[i]->Scale(nevents[h]);
	  h1_rej[h][i]->GetXaxis()->SetTitle("Threshold E [GeV]");
	  h1_rej[h][i]->GetYaxis()->SetTitle("Rejection Factor");
	  h1_rej[h][i]->SetMarkerStyle((h==1?40:20)+i);
	  h1_rej[h][i]->SetMarkerColor((i<3?2+i:3+i));
	  h1_rej[h][i]->SetMarkerSize(2);
	  h1_rej[h][i]->GetYaxis()->SetLabelSize(0.02);
	  h1_rej[h][i]->GetYaxis()->SetTitleSize(0.03);
	  h1_rej[h][i]->GetYaxis()->SetTitleOffset(1.7);
	  h1_rej[h][i]->GetXaxis()->SetLabelSize(0.02);
	  h1_rej[h][i]->GetXaxis()->SetTitleSize(0.03);
	  
	}
      h1_rej[0][0]->Draw("P");
      h1_rej[1][0]->Draw("SAME P");
      for(int i=1; i<4; ++i)
	{
	  for(int h=0; h<2; ++h)
	    {
	      h1_rej[h][i]->Draw("SAME P");
	    }
	}
      erejg->SetLineWidth(2);
      erejg->Draw("SAME");
    }
  TLegend* leg = new TLegend(0.12,0.6,0.4,0.85);
  for(int h=0; h<2; ++h)
    {
      string simstr = "Sim";
      string datstr = "Data";
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->AddEntry(h1_rej[h][0],((h==0?simstr:datstr)+" Rejection factor (no cuts)").c_str(),"p");
      leg->AddEntry(h1_rej[h][1],((h==0?simstr:datstr)+" Rejection factor (E_{EM} < "+to_string(maxeh)+" E_{had})").c_str(),"p");
      leg->AddEntry(h1_rej[h][2],((h==0?simstr:datstr)+" Rejection factor (E_{lead} < 0.65 E_{jet})").c_str(),"p");
      leg->AddEntry(h1_rej[h][3],((h==0?simstr:datstr)+" Rejection factor (both cuts)").c_str(),"p");
    }
  leg->AddEntry(erejg,"Expected rejection factor","l");
  leg->Draw();
  d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_rej"+".png").c_str());
  //hejet[h]->Scale(1./(nevents[h]));
  for(int h=0; h<2; ++h)
    {
      hejet[h]->Draw();
      d->SaveAs(("output/gmg/"+filebase+"_"+to_string(h)+"_"+to_string(cancount)+"_hejet.png").c_str());
    }
  cout << "cancount/ntree/nevents/nmb: " << cancount << " " << tree[1]->GetEntries() << " " << nevents[1]  << " " << nmb[1] << endl;
  cout << "Rejection factors:" << endl;
  
  
  for(int i=0; i<9; ++i)
    {
      cout << i+4;
      for(int j=0; j<4; ++j)
	{
	  cout << " " << h1_rej[1][j]->GetBinContent(10*(i+4)+1) << " " << h1_rej[1][j]->GetBinContent(10*(i+4)+1)/erejf[i];
	}
      cout << endl;
    }
  return 0;
}
