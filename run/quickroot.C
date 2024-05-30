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
  TChain* tree2 = new TChain("ttree2");
  TChain* simt = new TChain("ttree");
  TChain* simt2 = new TChain("ttree2");
  ifstream list;
  string line;
  ifstream simlist;
  simlist.open("sim.list",ifstream::in);
  if(!simlist)
    {
      cout << "nosimlist" << endl;
      exit(1);
    }
  while(getline(simlist,line))
    {
      try
	{
	  simt->Add(line.c_str());
	  simt2->Add(line.c_str());
	}
      catch(...)
	{
	  continue;
	}
    }
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
	  tree2->Add(line.c_str());
	}
      catch(...)
	{
	  continue;
	}
    }
  long long unsigned int nevents = 0;
  int evtct;
  int mbevt;
  long long unsigned int nmb = 0;
  tree2->SetBranchAddress("_evtct",&evtct);
  tree2->SetBranchAddress("mbevt",&mbevt);
  for(int i=0; i<tree2->GetEntries(); ++i)
    {
      tree2->GetEntry(i);
      nevents += evtct;
      nmb += mbevt;
      //cout << mbevt << endl;
    }
  cout << "Total " << nevents << " events." << endl;
  cout << "with " << nmb << " minbias events." << endl;
  const int ncol = 9;
  double red[ncol] = {1, 1, 1, 1, 1, 1, 1, 1};
  double grn[ncol] = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double blu[ncol] = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double stp[ncol] = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
  TColor::CreateGradientColorTable(ncol, stp, red, grn, blu, ncol);
  TH1D* h1_rej[4];
  for(int i=0; i<4; ++i)
    {
      h1_rej[i] = new TH1D(("h1_rej"+to_string(i)).c_str(),"",150,0,15);
      h1_rej[i]->GetYaxis()->SetRangeUser(0.5,200000);
    }
  TH2D* event_display = new TH2D("event_display","",96,-0.5,95.5,256,-0.5,255.5);
  TH2D* event_diphcal = new TH2D("event_display_hc","",24,0,24,64,0,64);
  TH1D* h1_dphi[4];
  for(int i=0; i<4; ++i)
    {
      h1_dphi[i] = new TH1D(("h1_dphi"+to_string(i)).c_str(),"",32,0,2*M_PI);
    }
  TH2D* h2_dphi_E = new TH2D("h2_dphi_E","",32,0,2*M_PI,50,0,50);
  TH1D* h1_AJ[4];
  for(int i=0; i<4; ++i)
    {
      h1_AJ[i] = new TH1D(("h1_AJ"+to_string(i)).c_str(),"",100,0,1);
    }
  TH2D* event_disrt[3];
  TH2D* event_sum = new TH2D("event_sum","Calorimeter Sum",24,0,24,64,0,64);
  TH1D* jetfrac[3];
  TH1D* jetE[4];
  TH1D* hejet = new TH1D("hejet","",100,0,10);
  hejet->SetMarkerStyle(43);
  hejet->SetMarkerSize(3);
  hejet->SetMarkerColor(kMagenta+3);
  for(int i=0; i<3; ++i)
    {
      string calo = "EMCal";
      if(i==1) calo = "IHCal";
      if(i==2) calo = "OHCal";
      jetE[i] = new TH1D(("jetE"+to_string(i)).c_str(),"",200,0,20);
      jetfrac[i] = new TH1D(("jetfrac"+to_string(i)).c_str(),"",100,0,1);
      event_disrt[i] = new TH2D(("event_display_rt"+to_string(i)).c_str(),calo.c_str(),24,0,24,64,0,64);
      jetfrac[i]->SetMarkerStyle(43);
      jetfrac[i]->SetMarkerColor(kMagenta+3);
      jetfrac[i]->SetMarkerSize(3);
      jetfrac[i]->GetYaxis()->SetLabelSize(0.02);
      jetfrac[i]->GetYaxis()->SetTitleSize(0.03);
      jetfrac[i]->GetYaxis()->SetTitleOffset(1.7);
      jetfrac[i]->GetXaxis()->SetLabelSize(0.02);
      jetfrac[i]->GetXaxis()->SetTitleSize(0.03);
      jetE[i]->SetMarkerStyle(43);
      jetE[i]->SetMarkerSize(3);
      jetE[i]->SetMarkerColor(kMagenta+3);
      jetE[i]->SetMarkerStyle(43);
      jetE[i]->SetMarkerColor(kMagenta+3);
      jetE[i]->SetMarkerSize(3);
      jetE[i]->GetYaxis()->SetLabelSize(0.02);
      jetE[i]->GetYaxis()->SetTitleSize(0.03);
      jetE[i]->GetYaxis()->SetTitleOffset(1.7);
      jetE[i]->GetXaxis()->SetLabelSize(0.02);
      jetE[i]->GetXaxis()->SetTitleSize(0.03);
    }
  
  jetE[3] = new TH1D(("jetE"+to_string(3)).c_str(),"",200,0,20);
  jetE[3]->SetMarkerColor(kMagenta+3);
  jetE[3]->SetMarkerSize(3);
  jetE[3]->SetMarkerStyle(43);
  jetE[3]->SetMarkerColor(kMagenta+3);
  jetE[3]->GetYaxis()->SetLabelSize(0.02);
  jetE[3]->GetYaxis()->SetTitleSize(0.03);
  jetE[3]->GetYaxis()->SetTitleOffset(1.7);
  jetE[3]->GetXaxis()->SetLabelSize(0.02);
  jetE[3]->GetXaxis()->SetTitleSize(0.03);
  jetfrac[0]->GetYaxis()->SetTitle("Event Normalized Counts 4 < E_{jet} < 7.5");
  jetfrac[0]->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
  jetfrac[1]->GetYaxis()->SetTitle("Event Normalized Counts 7.5 < E_{jet} < 10");
  jetfrac[1]->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
  jetfrac[2]->GetYaxis()->SetTitle("Event Normalized Counts 10 < E_{jet}");
  jetfrac[2]->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
  jetE[0]->GetYaxis()->SetTitle("Event Normalized Counts");
  jetE[0]->GetXaxis()->SetTitle("E_{jet} No Cuts");
  jetE[1]->GetYaxis()->SetTitle("Event Normalized Counts");
  jetE[1]->GetXaxis()->SetTitle("E_{jet}  with E_{EMCal} < 10E_{HCal}");
  jetE[2]->GetYaxis()->SetTitle("Event Normalized Counts");
  jetE[2]->GetXaxis()->SetTitle("E_{jet}  with E_{lead tower} < 0.65E_{jet}");
  jetE[3]->GetYaxis()->SetTitle("Event Normalized Counts");
  jetE[3]->GetXaxis()->SetTitle("E_{jet}  with E_{EMCal} < 10E_{HCal} \& E_{lead tower} < 0.65E_{jet}");
  hejet->GetXaxis()->SetTitle("E_{jet, EMCal} / E_{jet, HCals}");
  hejet->GetYaxis()->SetTitle("Event Normalized Counts");
  hejet->SetMarkerStyle(43);
  hejet->SetMarkerColor(kMagenta+3);
  hejet->SetMarkerSize(3);
  hejet->GetYaxis()->SetLabelSize(0.02);
  hejet->GetYaxis()->SetTitleSize(0.03);
  hejet->GetYaxis()->SetTitleOffset(1.7);
  hejet->GetXaxis()->SetLabelSize(0.02);
  hejet->GetXaxis()->SetTitleSize(0.03);
  float emetot;
  float seedD[1000];
  float ehjet[1000];
  long unsigned int trigvec;
  tree->SetBranchAddress("triggervec",&trigvec);
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
  for(int i=0; i<tree->GetEntries(); ++i)
    {
      highejet = 0;
      shighe = 0;
      sshighe = 0;
      int passcut = 1;
      //if(i % 1000 == 0) cout << i << endl;
      tree->GetEntry(i);
      if(!njet) continue;
      if(njet > 4) continue;
      //event_display->Reset();
      //int breakvar = 0;
      float minjet[4] = {9999,9999,9999,9999};
      //float minjetseedD[4] = {0};
      //float minjeteh[4] = {0};
      if(trigvec >> 10 & 1)
	{
      for(int j=0; j<njet; ++j)
	{
	  for(int k=0; k<4; ++k)
	    {
	      if(jet_e[j] < minjet[k])
		{
		  if(k==0)
		    {
		      minjet[k] = jet_e[j];
		      //minjetseedD[k] = seedD[j];
		      //minjeteh[k] = ehjet[j];
		    }
		  else if(k==1 && ehjet[j] < maxeh)
		    {
		      minjet[k] = jet_e[j];
		      //minjetseedD[k] = seedD[j];
		      //minjeteh[k] = ehjet[j];
		    }
		  else if(k==2 && seedD[j] < 0.65)
		    {
		      minjet[k] = jet_e[j];
		      //minjetseedD[k] = seedD[j];
		      //minjeteh[k] = ehjet[j];
		    }
		  else if(k==3 && seedD[j] <0.65 && ehjet[j] < maxeh)
		    {
		      minjet[k] = jet_e[j];
		      //minjetseedD[k] = seedD[j];
		      //minjeteh[k] = ehjet[j];
		    }
		}
	    }
	}
	}
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
	      h1_rej[j]->Fill(h1_rej[j]->GetBinCenter(it+1));
	      ++it;
	    }
	}
      //if(breakvar) continue;
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
	  for(int k=0; k<njet; ++k)
	    {
	      /*
	      if(get_phi(jet_ph[k]) > 29.5 && get_phi(jet_ph[k]) < 37.5)
		{
		  continue;
		}
	      */
	      if(ehjet[k] > maxeh || ehjet[k] < 0) continue;
	      if(seedD[k] > 0.65) continue;
	      if(jet_e[k] > 10)
		{
		  highejet = 1;
		}
	      if(jet_e[k] > 20) shighe = 1;
	      if(jet_e[k] > 25)
		{
		  sshighe = 1;
		  save_etaphiE(etart[0], phirt[0], enrt[0], sectorrt[0], "output/json/"+filebase+"_"+to_string(i)+"_em.txt", 0);
		  save_etaphiE(etart[1], phirt[1], enrt[1], sectorrt[1], "output/json/"+filebase+"_"+to_string(i)+"_ih.txt", 0);
		  save_etaphiE(etart[2], phirt[2], enrt[2], sectorrt[2], "output/json/"+filebase+"_"+to_string(i)+"_oh.txt", 0);
		}
	      jets[k] = new TMarker(get_eta(jet_et[k]),get_phi(jet_ph[k]),43);
	      jets[k]->SetMarkerSize(jet_e[k]/3);
	      jets[k]->SetMarkerColor(kRed);
	      //jets[k]->Draw();
	      for(int l=0; l<ncircle; ++l)
		{
		  float eta = get_eta(jet_et[k]+0.4*cos(2*l*M_PI/ncircle));
		  float phi = get_phi(jet_ph[k]+0.4*sin(2*l*M_PI/ncircle));
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
      for(int k=0; k<njet; ++k)
	{
	  

	  /*
	  if(get_phi(jet_ph[k]) > 29.5 && get_phi(jet_ph[k]) < 37.5)
	    {
	      continue;
	    }
	  */
	  if(jet_e[k] < 7.5)
	    {
	      jetfrac[0]->Fill(seedD[k]);
	    }
	  else if(jet_e[k] < 10)
	    {
	      jetfrac[1]->Fill(seedD[k]);
	    }
	  else
	    {
	      jetfrac[2]->Fill(seedD[k]);
	    }
	  //if(seedD[k] > 0.65) continue;
	  hejet->Fill(ehjet[njet]);
	  jetE[0]->Fill(jet_e[k]);
	  if(ehjet[k] < maxeh && ehjet[k] > 0) jetE[1]->Fill(jet_e[k]);
	  if(seedD[k] < 0.65) jetE[2]->Fill(jet_e[k]);
	  if(seedD[k] < 0.65 && ehjet[k] < maxeh && ehjet[k] > 0)
	    {
	      jetE[3]->Fill(jet_e[k]);
	    }
	  
	  for(int l = 0; l<njet; ++l)
	    {
	      if(l==k) continue;
	      if(jet_e[k] < jet_e[l]) continue;
	      if(jet_e[k] > 10 && jet_e[l] > 5)
		{
		  float ET1 = jet_e[k]/cosh(jet_et[k]);
		  float ET2 = jet_e[l]/cosh(jet_et[l]);
		  float AJ = (ET1-ET2)/(ET1+ET2);
		  h1_AJ[0]->Fill(AJ);
		  if(seedD[k] < 0.65) h1_AJ[2]->Fill(AJ);
		  if(ehjet[k] < maxeh) h1_AJ[1]->Fill(AJ);
		  if(ehjet[k] < maxeh && seedD[k] < 0.65) h1_AJ[3]->Fill(AJ);
		  float dphi = jet_ph[k] - jet_ph[l];
		  if(dphi < 0) dphi+=2*M_PI;
		  //cout << dphi << endl;
		  h1_dphi[0]->Fill(dphi);
		  if(seedD[k] < 0.65) h1_dphi[2]->Fill(AJ);
		  if(ehjet[k] < maxeh) h1_dphi[1]->Fill(AJ);
		  if(ehjet[k] < maxeh && seedD[k] < 0.65) h1_dphi[3]->Fill(AJ);
		  h2_dphi_E->Fill(dphi,jet_e[k]);
		}
	    }
	  if(seedD[k] > 0.65) continue;
	  if(ehjet[k] > maxeh || ehjet[k] < 0) continue;
	  jets[k] = new TMarker(get_eta(jet_et[k]),get_phi(jet_ph[k]),20);
	  jets[k]->SetMarkerSize(jet_e[k]/3);
	  jets[k]->SetMarkerColor(kRed);
	  //jets[k]->Draw();
	  for(int l=0; l<ncircle; ++l)
	    {
	      float eta = get_eta(jet_et[k]+0.4*cos(2*l*M_PI/ncircle));
	      float phi = get_phi(jet_ph[k]+0.4*sin(2*l*M_PI/ncircle));
	      if(eta > 24 || eta < 0) continue;
	      if(phi > 64) phi -= 64;
	      if(phi < 0) phi += 64;
	      TMarker* circlemarker = new TMarker(eta,phi,20);
	      circlemarker->SetMarkerSize(0.3);
	      circlemarker->SetMarkerColor(kBlue);
	      circlemarker->Draw();
	    }
	  std::stringstream e_stream;
	  e_stream << std::fixed << std::setprecision(2) << jet_e[k];
	  std::string e_string = e_stream.str();
	  drawText((e_string+" GeV").c_str(),get_eta(jet_et[k]),get_phi(jet_ph[k]),(get_eta(jet_et[k])>15?1:0),kBlack,0.08,42,false);
	  //d->cd();
	  //event_sum->Draw("LEGO");
	}
      //if(dispcount % 18 == 17 && highejet)
      if(highejet)
	{
	  //cout << "Saving display " << dispcount/18 << endl;
	  //c->SaveAs(("./output/img/candidate_"+filebase+"_"+to_string(dispcount/18)+".pdf").c_str());
	  c->SaveAs(("./output/img/candidate_"+filebase+"_"+to_string(cancount)+".pdf").c_str());
	  //cout << "Saved!" << endl;
	  //c->Clear("D");
	}
      ++cancount;
      if(shighe)
	{
	  c->SaveAs(("./output/hmg/candidate_"+filebase+"_superhighE_"+to_string(cancount)+".pdf").c_str());
	  //d->SaveAs(("./output/hmg/candidate_"+filebase+"_superhighElego_"+to_string(cancount)+".pdf").c_str());
	}
      if(sshighe) 
	{
	  c->SaveAs(("./output/smg/candidate_"+filebase+"_supersuperhighE_"+to_string(cancount)+".pdf").c_str());
	  //d->SaveAs(("./output/smg/candidate_"+filebase+"_supersuperhighElego_"+to_string(cancount)+".pdf").c_str());
	}
      if(highejet)
	{
	  //cout << dispcount << endl;
	  ++dispcount;
	}
    }
  //if(dispcount % 18 != 0)
  //  {
  //    c->SaveAs(("./output/img/candidate_"+filebase+"_last.pdf").c_str());
  //  }
  d->cd();
  for(int i=0; i<3; ++i)
    {
      jetfrac[i]->Scale(1./(nevents));
    }
  //h1_dphi->SetMarkerColor(kMagenta+3);
  for(int i=0; i<4; ++i)
    {
      h1_dphi[i]->GetYaxis()->SetLabelSize(0.02);
      h1_dphi[i]->GetYaxis()->SetTitleSize(0.03);
      h1_dphi[i]->GetYaxis()->SetTitleOffset(1.7);
      h1_dphi[i]->GetXaxis()->SetLabelSize(0.02);
      h1_dphi[i]->GetXaxis()->SetTitleSize(0.03);
      h1_dphi[i]->Scale(1./nevents);
      h1_dphi[i]->GetYaxis()->SetTitle("Event Normalized Counts");
      h1_dphi[i]->GetXaxis()->SetTitle("#Delta#phi [rad]");
      h1_dphi[i]->SetFillColorAlpha(kMagenta,0.4);
      h1_dphi[i]->Draw("HIST");
      drawText("Leading Jet E > 10 GeV", 0.5, 0.85, 0, kBlack, 0.02);
      drawText("Subleading Jet E > 5 GeV", 0.5, 0.825, 0, kBlack, 0.02);
      if(i>1) drawText("E_{lead tower} < 0.65 E_{jet}", 0.5, 0.8, 0, kBlack, 0.02);
      if(i==1 || i==3) drawText(("E_{EM} < "+to_string(maxeh)+" E_{HAD}").c_str(), 0.5, (i==1?0.8:0.775), 0, kBlack, 0.02);
      d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_"+to_string(i)+"_dphi_1d.pdf").c_str());
      
      h1_AJ[i]->GetYaxis()->SetLabelSize(0.02);
      h1_AJ[i]->GetYaxis()->SetTitleSize(0.03);
      h1_AJ[i]->GetYaxis()->SetTitleOffset(1.7);
      h1_AJ[i]->GetXaxis()->SetLabelSize(0.02);
      h1_AJ[i]->GetXaxis()->SetTitleSize(0.03);
      h1_AJ[i]->Scale(1./nevents);
      h1_AJ[i]->GetYaxis()->SetTitle("Event Normalized Counts");
      h1_AJ[i]->GetXaxis()->SetTitle("A_{J}");
      h1_AJ[i]->SetFillColorAlpha(kMagenta,0.4);
      h1_AJ[i]->Draw("HIST");
      drawText("Leading Jet E > 10 GeV", 0.5, 0.85, 0, kBlack, 0.02);
      drawText("Subleading Jet E > 5 GeV", 0.5, 0.825, 0, kBlack, 0.02);
      if(i>1) drawText("E_{lead tower} < 0.65 E_{jet}", 0.5, 0.8, 0, kBlack, 0.02);
      if(i==1 || i==3) drawText(("E_{EM} < "+to_string(maxeh)+" E_{HAD}").c_str(), 0.5, (i==1?0.8:0.775), 0, kBlack, 0.02);
      d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_"+to_string(i)+"_AJ_1d.pdf").c_str());
    }
  //h2_dphi_E->SetMarkerColor(kMagenta+3);
  h2_dphi_E->GetYaxis()->SetLabelSize(0.02);
  h2_dphi_E->GetYaxis()->SetTitleSize(0.03);
  h2_dphi_E->GetYaxis()->SetTitleOffset(1.7);
  h2_dphi_E->GetXaxis()->SetLabelSize(0.02);
  h2_dphi_E->GetXaxis()->SetTitleSize(0.03);
  h2_dphi_E->Scale(1./nevents);
  h2_dphi_E->GetYaxis()->SetTitle("Jet E [GeV]");
  h2_dphi_E->GetXaxis()->SetTitle("#Delta#phi [rad]");
  h2_dphi_E->GetZaxis()->SetTitle("Event Normalized Counts");
  h2_dphi_E->Draw("COLZ");
  d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_dphi_2d.pdf").c_str());
  jetfrac[0]->Draw();
  d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_jetfrac0.pdf").c_str());
  jetfrac[1]->Draw();
  d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_jetfrac1.pdf").c_str());
  jetfrac[2]->Draw();
  d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_jetfrac2.pdf").c_str());
  d->SetLogy();
  //float goodfrac = 1 - fakejets/tree->GetEntries();
  float erejf[9] = {14.24,42.76,111.95,281.76,681.85,1466.6,3033.03,6272.51,11086.03};
  float x[9] = {4,5,6,7,8,9,10,11,12};
  TGraph* erejg = new TGraph(9,x,erejf);
  erejg->SetMarkerStyle(43);
  erejg->SetMarkerColor(kMagenta);
  for(int i=0; i<4; ++i)
    {
      jetE[i]->Scale(1./(nevents));
      jetE[i]->Draw();
      d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_jetE"+to_string(i)+".pdf").c_str());
      for(int j=0; j<160; ++j)
	{
	  if(h1_rej[i]->GetBinContent(j+1) != 0) h1_rej[i]->SetBinContent(j+1,nmb/h1_rej[i]->GetBinContent(j+1));
	}
      //      h1_rej[i]->Scale(nevents);
      h1_rej[i]->GetXaxis()->SetTitle("Threshold E [GeV]");
      h1_rej[i]->GetYaxis()->SetTitle("Rejection Factor");
      h1_rej[i]->SetMarkerStyle(40+i);
      h1_rej[i]->SetMarkerColor((i<3?2+i:3+i));
      h1_rej[i]->SetMarkerSize(2);
      h1_rej[i]->GetYaxis()->SetLabelSize(0.02);
      h1_rej[i]->GetYaxis()->SetTitleSize(0.03);
      h1_rej[i]->GetYaxis()->SetTitleOffset(1.7);
      h1_rej[i]->GetXaxis()->SetLabelSize(0.02);
      h1_rej[i]->GetXaxis()->SetTitleSize(0.03);

    }
  h1_rej[0]->Draw("P");
  for(int i=1; i<4; ++i)
    {
      h1_rej[i]->Draw("SAME P");
    }
  erejg->SetLineWidth(2);
  erejg->Draw("SAME");
  TLegend* leg = new TLegend(0.12,0.6,0.4,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(h1_rej[0],"Rejection factor (no cuts)","p");
  leg->AddEntry(h1_rej[1],("Rejection factor (E_{EM} < "+to_string(maxeh)+" E_{had})").c_str(),"p");
  leg->AddEntry(h1_rej[2],"Rejection factor (E_{lead} < 0.65 E_{jet})","p");
  leg->AddEntry(h1_rej[3],"Rejection factor (both cuts)","p");
  leg->AddEntry(erejg,"Expected rejection factor","l");
  leg->Draw();
  d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_rej"+".pdf").c_str());
  hejet->Scale(1./(nevents));
  hejet->Draw();
  d->SaveAs(("output/gmg/"+filebase+"_"+to_string(cancount)+"_hejet.pdf").c_str());
  cout << "cancount/ntree/nevents/nmb: " << cancount << " " << tree->GetEntries() << " " << nevents  << " " << nmb << endl;
  cout << "Rejection factors:" << endl;
  
  
  for(int i=0; i<9; ++i)
    {
      cout << i+4;
      for(int j=0; j<4; ++j)
	{
	  cout << " " << h1_rej[j]->GetBinContent(10*(i+4)+1) << " " << h1_rej[j]->GetBinContent(10*(i+4)+1)/erejf[i];
	}
      cout << endl;
    }
  return 0;
}
