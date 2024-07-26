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
      thefile << (i>0?",":"") <<"{\"eta\": " << (em?bintoeta_em(etabins[i]):bintoeta_hc(etabins[i])) << ", \"phi\": " << (em?bintophi_em(phibins[i]):bintophi_hc(phibins[i])) << ", \"e\": " << (energies[i]>0?energies[i]:0) << ", \"event\": 0}" << endl;
    }
  thefile.close();
  return 0;
}

int quickroot(string filebase="", int njob=0)
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
  //gStyle->SetPalette(kOcean);
  gROOT->SetStyle("Plain");
  //SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0);
  float emcalen[2][25000];
  float emcalchi2[2][25000];
  float emcalt[2][25000];
  int sectorrt[2][3];
  int sectorem[2];
  int njet[2];
  float enrt[2][3][24576];
  int etart[2][3][24576];
  int phirt[2][3][24576];
  int etabin[2][25000];
  int phibin[2][25000];
  float jet_e[2][1000];
  float jet_et[2][1000];
  float jet_ph[2][1000];
  string filename=filebase;
  TChain* tree[2];
  TChain* tree2[2];
  TChain* jett[2];
  ifstream list[2];
  string line;
  int nosim = 1;

  list[1].open(filename,ifstream::in);
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
  for(int h=nosim; h < 2; ++h)
    {
      int counter = 0;
      jett[h] = new TChain("jett");
      tree[h] = new TChain("ttree");
      tree2[h] = new TChain("ttree2");
      if(!list[h])
	{
	  cout << "nosimlist" << endl;
	  exit(1);
	}
      /*
      while(getline(list[h],line))
	{
	}
      */
      for(int i=0; i<10*njob; ++i)
	{
	  getline(list[h],line);
	}
      for(int i=10*njob; i<10*(njob+1); ++i)
	{
	  int breakit = 0;
	  getline(list[h],line);
	  if(list[h].eof()) breakit = 1;
	  try
	    {
	      jett[h]->Add(line.c_str());
	      tree[h]->Add(line.c_str());
	      tree2[h]->Add(line.c_str());
	    }
	  catch(...)
	    {
	      continue;
	    }
	  if(breakit) break;
	}
    }
  //cout << "test-1" << endl;
  long long unsigned int nevents[2] = {0};
  int evtct[2];
  int mbevt[2];
  long long unsigned int nmb[2] = {0};
  for(int h=nosim; h<2; ++h)
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
      //if(h==0) nmb[h] = nevents[h];
      cout << "with " << nmb[h] << " minbias events." << endl;
    }
  TFile* outfile = TFile::Open(("output/root/run_"+runnum+"_"+idstr+"_"+to_string(njob)+"_fullfile.root").c_str(),"RECREATE");
  cout << outfile->GetName() << endl;
  TTree* outt = new TTree((runnum==simstr?"st":"outt"),"output tree");
  int outnmb = nmb[1];
  int outevt = nevents[1];
  int njmb = nmb[1];
  int nJetTrig[4] = {0};
  outt->Branch("outnmb",&outnmb,"outnmb/I");
  outt->Branch("outevt",&outevt,"outevt/I");
  outt->Branch("njmb",&njmb,"njmb/I");
  outt->Branch("nJetTrig",nJetTrig,"nJetTrig[4]/I");
  const int ncol = 9;
  //cout << "test0" << endl;
  double red[ncol] = {1, 1, 1, 1, 1, 1, 1, 1};
  double grn[ncol] = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double blu[ncol] = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double stp[ncol] = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
  TColor::CreateGradientColorTable(ncol, stp, red, grn, blu, ncol);
  TH1D* h1_rej[2][4];
  TH1D* h1_eta[2][4][5];
  TH1D* h1_phi[2][4][5];
  TH1D* h1_mlt[2][4];
  TH1D* h1_calo_occ[3];
  TH1D* h1_tower_E[3];
  TH1D* h1_calo_E[3];
  TH2D* h2_occ_E[3];
  TH2D* h2_jet_eta_phi[2][5];
  TH2D* h2_jet_eta_e[2][4];
  TH2D* h2_cal_eta_phi[2][3];
  TH1D* h1_cluster_E = new TH1D("h1_cluster_E","",40,0,40);
  TH1D* h1_cluster_eta = new TH1D("h1_cluster_eta","",24,-1.1,1.1);
  TH2D* h2_cluster_eta_E = new TH2D("h2_cluster_eta_E","",24,-1.1,1.1,40,0,40);
  TH2D* h2_cluster_eta_phi = new TH2D("h2_cluster_eta_phi","",24,-1.1,1.1,64,-1.1,1.1);
  for(int h=nosim; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  if(i<3)
	    {
	      h1_calo_occ[i] = new TH1D(("h1_calo_occ"+to_string(i)).c_str(),"",100,0,0.1);
	      h1_tower_E[i] = new TH1D(("h1_tower_E"+to_string(i)).c_str(),"",100,-5,(i==1?5:20));
	      h1_calo_E[i] = new TH1D(("h1_calo_E"+to_string(i)).c_str(),"",100,-5,(i==1?10:50));
	      h2_occ_E[i] = new TH2D(("h2_occ_E"+to_string(i)).c_str(),"",100,0,0.1,100,0,(i==1?10:50));
	    }
	  if(i<3 && i>0) h2_cal_eta_phi[h][i] = new TH2D(("h2_cal_eta_phi"+to_string(h)+to_string(i)).c_str(),"",24,-0.5,23.5,64,-0.5,63.5);
	  if(i==0) h2_cal_eta_phi[h][i] = new TH2D(("h2_cal_eta_phi"+to_string(h)+to_string(i)).c_str(),"",96,-0.5,95.5,256,-0.5,255.5);

	  h2_jet_eta_e[h][i] = new TH2D(("h2_jet_eta_e"+to_string(h)+to_string(i)).c_str(),"",24,-1.1,1.1,160,4,20);
	  h1_rej[h][i] = new TH1D(("h1_rej"+to_string(h)+"_"+to_string(i)).c_str(),"",150,0,15);
	  h1_rej[h][i]->GetYaxis()->SetRangeUser(0.5,200000);
	  for(int j=0; j<5; ++j)
	    {
	      if(i==0) h2_jet_eta_phi[h][j] = new TH2D(("h2_jet_eta_phi"+to_string(h)+to_string(j)).c_str(),"",24,-1.1,1.1,64,-M_PI,M_PI);
	      h1_eta[h][i][j] = new TH1D(("h1_eta"+to_string(h)+"_"+to_string(i)+"_"+to_string(j)).c_str(),"",24,-1.1,1.1);
	      h1_phi[h][i][j] = new TH1D(("h1_phi"+to_string(h)+"_"+to_string(i)+"_"+to_string(j)).c_str(),"",30,-M_PI,M_PI);
	    }
	  h1_mlt[h][i] = new TH1D(("h1_mlt"+to_string(h)+"_"+to_string(i)).c_str(),"",6,-0.5,5.5);
	}
    }
  TH2D* event_display = new TH2D("event_display","",96,-0.5,95.5,256,-0.5,255.5);
  TH2D* event_diphcal = new TH2D("event_display_hc","",24,0,24,64,0,64);
  TH1D* h1_dphi[2][4];
  for(int h=nosim; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  h1_dphi[h][i] = new TH1D(("h1_dphi"+to_string(h)+"_"+to_string(i)).c_str(),"",32,0,2*M_PI);
	}
    }
  TH2D* h2_dphi_deta[2];
  TH1D* h1_AJ[2][12];
  TH2D* h2_tjet_eta_phi = new TH2D("h2_tjet_eta_phi","",22,-1.1,1.1,64,-M_PI,M_PI);
    for(int h=nosim; h<2; ++h)
    {
      h2_dphi_deta[h] = new TH2D(("h2_dphi_deta"+to_string(h)).c_str(),"",20,-2.2,2.2,32,0,2*M_PI);
      for(int i=0; i<4; ++i)
	{
	  h1_AJ[h][i] = new TH1D(("h1_AJ"+to_string(h)+"_"+to_string(i)).c_str(),"",50,0,1);
	  h1_AJ[h][i+4] = new TH1D(("h1_AJ"+to_string(h)+"_"+to_string(i+4)).c_str(),"",50,0,1);
	  h1_AJ[h][i+8] = new TH1D(("h1_AJ"+to_string(h)+"_"+to_string(i+8)).c_str(),"",50,0,1);
	}
    }
  TH2D* event_disrt[3];
  TH2D* event_sum = new TH2D("event_sum","Calorimeter Sum",24,0,24,64,0,64);
  TH1D* jetfrac[2][3];
  TH1D* jetE[2][4];
  TH1D* ljetE = new TH1D("ljetE","",400,0,40);
  TH1D* ljetEta = new TH1D("ljetEta","",20,-1.1,1.1);
  TH1D* tjetE = new TH1D("tjetE","",400,0,40);
  TH1D* hejet[2];
  for(int h=nosim; h<2; ++h)
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
	  jetE[h][i] = new TH1D(("jetE"+to_string(h)+"_"+to_string(i)).c_str(),"",400,0,40);
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
  
      jetE[h][3] = new TH1D(("jetE"+to_string(3)).c_str(),"",400,0,40);
      jetE[h][3]->SetMarkerColor(kMagenta+3);
      jetE[h][3]->SetMarkerSize(3);
      jetE[h][3]->SetMarkerStyle(43);
      jetE[h][3]->SetMarkerColor(kMagenta+3);
      jetE[h][3]->GetYaxis()->SetLabelSize(0.02);
      jetE[h][3]->GetYaxis()->SetTitleSize(0.03);
      jetE[h][3]->GetYaxis()->SetTitleOffset(1.7);
      jetE[h][3]->GetXaxis()->SetLabelSize(0.02);
      jetE[h][3]->GetXaxis()->SetTitleSize(0.03);
      jetfrac[h][0]->GetYaxis()->SetTitle("Counts 4 < E_{jet} < 7.5");
      jetfrac[h][0]->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
      jetfrac[h][1]->GetYaxis()->SetTitle("Counts 7.5 < E_{jet} < 10");
      jetfrac[h][1]->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
      jetfrac[h][2]->GetYaxis()->SetTitle("Counts 10 < E_{jet}");
      jetfrac[h][2]->GetXaxis()->SetTitle("E_{leading tower} / E_{jet}");
      jetE[h][0]->GetYaxis()->SetTitle("Event Normalized Counts");
      jetE[h][0]->GetXaxis()->SetTitle("E_{T,jet} No Cuts");
      jetE[h][1]->GetYaxis()->SetTitle("Event Normalized Counts");
      jetE[h][1]->GetXaxis()->SetTitle("E_{T,jet}  with E_{EMCal} < 10E_{HCal}");
      jetE[h][2]->GetYaxis()->SetTitle("Event Normalized Counts");
      jetE[h][2]->GetXaxis()->SetTitle("E_{T,jet}  with E_{lead tower} < 0.65E_{jet}");
      jetE[h][3]->GetYaxis()->SetTitle("Event Normalized Counts");
      jetE[h][3]->GetXaxis()->SetTitle("E_{T,jet}  with E_{EMCal} < 10E_{HCal} \& E_{lead tower} < 0.65E_{jet}");
      hejet[h]->GetXaxis()->SetTitle("E_{jet, EMCal} / E_{jet, HCals}");
      hejet[h]->GetYaxis()->SetTitle("Counts");
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
  int evtnum[2];
  long unsigned int trigvec = 0;
  int ismb[2];
  float alcet[24000];
  float aceta[1000];
  int njetj = 0;
  int allcomp = 0;
  float vtx[3];
  float mbdq[128];
  int ntj;
  float tjet_e[1000];
  float tjet_eta[1000];
  float tjet_phi[1000];
  float caloE[3];
  float clusterE[10000];
  float clustereta[10000];
  float clusterphi[10000];
  int cluster_n;
  tree[1]->SetBranchAddress("triggervec",&trigvec);
  tree[1]->SetBranchAddress("emetot",&caloE[0]);
  tree[1]->SetBranchAddress("ihetot",&caloE[1]);
  tree[1]->SetBranchAddress("ohetot",&caloE[2]);
  tree[1]->SetBranchAddress("cluster_n",&cluster_n);
  tree[1]->SetBranchAddress("cluster_E",clusterE);
  tree[1]->SetBranchAddress("cluster_eta",clustereta);
  tree[1]->SetBranchAddress("cluster_phi",clusterphi);
  for(int h=nosim; h<2; ++h)
    {
      tree[h]->SetBranchAddress("mbenrgy",mbdq);
      tree[h]->SetBranchAddress("vtx",vtx);
      tree[h]->SetBranchAddress("ismb",&ismb[h]);
      tree[h]->SetBranchAddress("_evtnum",&evtnum[h]);
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
      tree[h]->SetBranchAddress("ntj",&ntj);
      tree[h]->SetBranchAddress("tjet_e",tjet_e);
      tree[h]->SetBranchAddress("tjet_eta",tjet_eta);
      tree[h]->SetBranchAddress("tjet_phi",tjet_phi);
      //tree->SetBranchAddress("emcalen",emcalen);
      //tree->SetBranchAddress("emcalt",emcalt);
      //tree->SetBranchAddress("emcalchi2",emcalchi2);
      //tree->SetBranchAddress("sectorem",&sectorem);
      //tree->SetBranchAddress("emcaletabin",etabin);
      //tree->SetBranchAddress("emcalphibin",phibin);
      //tree->SetBranchAddress("emetot",&emetot);
      jett[h]->SetBranchAddress("alcet",alcet);
      jett[h]->SetBranchAddress("aceta",aceta);
      jett[h]->SetBranchAddress("allcomp",&allcomp);
      jett[h]->SetBranchAddress("njet",&njetj);
    }
  TFile* jetfile = TFile::Open(("output/root/run_"+runnum+"_"+idstr+"_"+to_string(njob)+"_jetfile.root").c_str(),"RECREATE");
  TH1D* alcethist = new TH1D("alcethist","",40,-2,2);
  TH1D* acetahist = new TH1D("acetahist","",40,-2,2);
  TTree* ojt = new TTree("ojt","the duplicated jet tree");
  ojt->Branch("ismb",&ismb[1],"ismb/I");
  ojt->Branch("njet",&njet[1],"njet/I");
  ojt->Branch("jet_e",jet_e[1],"jet_e[njet]/F");
  ojt->Branch("jet_et",jet_et[1],"jet_et[njet]/F");
  ojt->Branch("jet_ph",jet_ph[1],"jet_ph[njet]/F");
  for(int i=0; i<jett[1]->GetEntries(); ++i)
    {
      jett[1]->GetEntry(i);
      for(int j=0; j<njetj; ++j)
	{
	  acetahist->Fill(aceta[j]);
	}
      for(int j=0; j<allcomp; ++j)
	{
	  alcethist->Fill(alcet[j]);
	}
    }
  
  jetfile->WriteObject(acetahist,acetahist->GetName());
  jetfile->WriteObject(alcethist,alcethist->GetName());
  int cancount = 0;
  TCanvas* c = new TCanvas("","",800,480);
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
      event_disrt[i]->GetZaxis()->SetTitleOffset(0.75);
      event_disrt[i]->GetYaxis()->SetTitleOffset(1);
      event_disrt[i]->GetZaxis()->SetTitle("EM Scale Tower Energy [GeV]");
      event_disrt[i]->GetZaxis()->SetRangeUser(0.1,5);
      event_disrt[i]->GetXaxis()->SetNdivisions(6,kFALSE);
      event_disrt[i]->GetXaxis()->SetTitleSize(0.08);
      event_disrt[i]->GetYaxis()->SetTitleSize(0.08);
      event_disrt[i]->GetZaxis()->SetTitleSize(0.08);
      event_disrt[i]->GetXaxis()->SetTitleOffset(0.6);
      event_disrt[i]->GetXaxis()->SetLabelSize(0.08);
      event_disrt[i]->GetYaxis()->SetLabelSize(0.08);
      event_disrt[i]->GetZaxis()->SetLabelSize(0.08);
      event_disrt[i]->GetXaxis()->SetLabelOffset(-0.02);
    }
  event_sum->GetXaxis()->SetTitle("Tower Sum #eta Bin");
  event_sum->GetYaxis()->SetTitle("Tower Sum #phi Bin");
  event_sum->GetZaxis()->SetTitle("EM Scale Tower Energy [GeV]");
  event_sum->GetYaxis()->SetTitleOffset(1);
  event_sum->GetZaxis()->SetTitleOffset(0.75);
  event_sum->GetZaxis()->SetRangeUser(0.1,5);
  event_sum->GetXaxis()->SetNdivisions(6,kFALSE);
  event_sum->GetXaxis()->SetTitleSize(0.08);
  event_sum->GetZaxis()->SetTitleSize(0.08);
  event_sum->GetXaxis()->SetTitleOffset(0.6);
  event_sum->GetYaxis()->SetTitleSize(0.08);
  event_sum->GetXaxis()->SetLabelSize(0.08);
  event_sum->GetYaxis()->SetLabelSize(0.08);
  event_sum->GetZaxis()->SetLabelSize(0.08);
  event_sum->GetXaxis()->SetLabelOffset(-0.02);
  string subdet[4] = {"#bf{EMCal}","#bf{IHCal}","#bf{OHCal}","#bf{Sum}"};
  int ncircle = 64;
  int highejet = 0;
  float fakejets = 0;
  int dispcount = 1;
  int shighe = 0;
  int sshighe = 0;
  int maxeh = 5;
  int eventbase[2] = {0};
  int prevt = 0;
  TH1D* jetTrigE[4];
  for(int i=0; i<4; i++)
    {
      //for(int j=0; j<4; ++j)
      //{
      jetTrigE[i] = new TH1D(("jetTrigE"+to_string(i)).c_str(),"",400,0,40);
	  //}
    }
  TH1D* h1_zdist = new TH1D("h1_zdist","",200,-100,100);
  TH1D* h1_mbdq = new TH1D("h1_mbdq","",1000,0,10);
  for(int h=nosim; h<2; ++h)
    {
      cancount = 0;
  for(int i=0; i<tree[h]->GetEntries(); ++i)
    {
      tree[h]->GetEntry(i);
      //highejet = 0;
      //shighe = 0;
      //sshighe = 0;
      //int passcut = 1;
      //if(i % 1000 == 0) cout << i << endl;

      if(vtx[2] > 100)
	{
	  if(ismb[h]) outnmb--;
	  continue;
	}

      if(h==1) 
	{
      if((trigvec >> 16) & 1 || (trigvec >> 17) & 1 || (trigvec >> 18) & 1 || (trigvec >> 19) & 1)
	{
	  njmb++;
	  for(int j=0; j<njet[h]; ++j)
	    {
	      float ET = jet_e[h][j]/cosh(jet_et[h][j]);
	      if(ET<4) continue;
	      if(abs(jet_et[h][j]) > 1.1) continue;
	      if((trigvec >> 16) & 1)
		{
		  nJetTrig[0]++;
		  jetTrigE[0]->Fill(ET);
		}
	      if((trigvec>>17) &1)
		{
		  nJetTrig[1]++;
		  jetTrigE[1]->Fill(ET);
		}
	      if((trigvec>>18) &1)
		{
		  nJetTrig[2]++;
		  jetTrigE[2]->Fill(ET);
		}
	      if((trigvec>>19) &1)
		{
		  nJetTrig[3]++;
		  jetTrigE[3]->Fill(ET);
		}
	    }
	}
	}
      ojt->Fill();
      int trigvecproxy = 0;
      if(filename != simliststr) trigvecproxy = (trigvec>>10) & 1;
      if(!ismb[h])// && !trigvecproxy && filename != simliststr)
	{
	  //cout << filename << " " << simliststr << " " << to_string(ismb[h]) << endl;
	  continue;
	}
      h1_zdist->Fill(vtx[2]);
      for(int j=0; j<128; ++j)
	{
	  h1_mbdq->Fill(mbdq[j]);
	}

      for(int j=0; j<cluster_n; ++j)
	{
	  h1_cluster_E->Fill(clusterE[j]);
	  if(abs(clusterphi[j]) > 0.25) h1_cluster_eta->Fill(clustereta[j]);
	  h2_cluster_eta_E->Fill(clustereta[j],clusterE[j]);
	  h2_cluster_eta_phi->Fill(clustereta[j],clusterphi[j]);
	}
      if(filename == simliststr)
	{
	  //	  cout << "madeit" << endl;
	  for(int j=0; j<ntj; ++j)
	    {
	      float ET = tjet_e[j]/cosh(tjet_eta[j]);
	      if(ET < 4) continue;
	      if(abs(tjet_eta[j]) > 1.1) continue;
	      tjetE->Fill(ET);
	      h2_tjet_eta_phi->Fill(tjet_eta[j],tjet_phi[j]);
	    }
	}
      for(int j=0; j<3; ++j)
	{
	  float calo_E = 0;
	  float occ = 0;
	  
	  for(int k=0; k<sectorrt[h][j]; ++k)
	    {
	      h1_tower_E[j]->Fill(enrt[h][j][k]);
	      calo_E += enrt[h][j][k];
	      if(enrt[h][j][k] > 0.03) occ+=1;
	    }
	  
	  calo_E = caloE[j];
	  h1_calo_E[j]->Fill(calo_E);
	  h1_calo_occ[j]->Fill(occ/(j==0?24576:1536));
	  h2_occ_E[j]->Fill(occ/(j==0?24576:1536),calo_E);
	}
      if(evtnum[h] < prevt)
	{
	  eventbase[h] += 10000;
	}
      //if(!njet[h]) continue;
      //if(njet[h] > 4) continue;
      //event_display->Reset();
      //int breakvar = 0;
      float maxjet[4] = {0};
      //float maxjetseedD[4] = {0};
      //float maxjeteh[4] = {0};
      //if(trigvec >> 10 & 1)
      //{
      for(int j=0; j<njet[h]; ++j)
	{

	  for(int k=0; k<4; ++k)
	    {
	      if(jet_e[h][j] > maxjet[k])
		{
		  if(k==0)
		    {
		      maxjet[k] = jet_e[h][j];
		      //maxjetseedD[k] = seedD[j];
		      //maxjeteh[k] = ehjet[j];
		    }
		  if(k==1 && ehjet[h][j] < maxeh)
		    {
		      maxjet[k] = jet_e[h][j];
		      //maxjetseedD[k] = seedD[j];
		      //maxjeteh[k] = ehjet[j];
		    }
		  if(k==2 && seedD[h][j] < 0.65)
		    {
		      maxjet[k] = jet_e[h][j];
		      //maxjetseedD[k] = seedD[j];
		      //maxjeteh[k] = ehjet[j];
		    }
		  if(k==3 && seedD[h][j] <0.65 && ehjet[h][j] < maxeh)
		    {
		      maxjet[k] = jet_e[h][j];
		      //maxjetseedD[k] = seedD[j];
		      //maxjeteh[k] = ehjet[j];
		    }
		}
	    }
	}
      //}
      for(int j=0; j<4; ++j)
	{
	  if(maxjet[j] < 4)
	    {
	      maxjet[j] = 0;
	    }
	}
      for(int j=0; j<4; ++j)
	{
	  int it = 40;
	  while(it < 10*maxjet[j])
	    {
	      h1_rej[h][j]->Fill(h1_rej[h][j]->GetBinCenter(it+1));
	      ++it;
	    }
	}
      
      for(int j=0; j<3; ++j)
	{
	  for(int k=0; k<sectorrt[h][j]; ++k)
	    {
	      h2_cal_eta_phi[h][j]->Fill(etart[h][j][k],phirt[h][j][k],enrt[h][j][k]);
	    }
	}
      
      int njo[4] = {0};	    
      float leadjetET = 0;
      float leadjetEta = -2;
      for(int k=0; k<njet[h]; ++k)
	{
	  float ET = jet_e[h][k]/cosh(jet_et[h][k]);
	  if(abs(jet_et[h][k]) > 1.1) continue;
	  if(ET < 4) continue;
	  if(ET > 4) njo[0]++;
	  if(ET > 5) njo[1]++;
	  if(ET > 7) njo[2]++;
	  if(ET > 10) njo[3]++;
	  if(ET > leadjetET)
	    {
	      leadjetET = ET;
	      leadjetEta = jet_et[h][k];
	    }
	  /*
	  if(get_phi(jet_ph[k]) > 29.5 && get_phi(jet_ph[k]) < 37.5)
	    {
	      continue;
	    }
	  */
	  if(jet_e[h][k] < 7)
	    {
	      jetfrac[h][0]->Fill(1.-1./ehjet[h][k]);
	    }
	  else if(jet_e[h][k] < 10)
	    {
	      jetfrac[h][1]->Fill(1.-1./ehjet[h][k]);
	    }
	  else
	    {
	      jetfrac[h][2]->Fill(1.-1./ehjet[h][k]);
	    }
	  //if(seedD[k] > 0.65) continue;
	  
	  hejet[h]->Fill(ehjet[h][k]);
	  jetE[h][0]->Fill(ET);
	  int whichjet = 0;
	  if(ET < 6) whichjet = 0;
	  else if(ET < 10) whichjet = 1;
	  else if(ET < 15) whichjet = 2;
	  else if(ET < 20) whichjet = 3;
	  else whichjet = 4;
	  h1_eta[h][0][whichjet]->Fill(jet_et[h][k]);
	  h1_phi[h][0][whichjet]->Fill(jet_ph[h][k]);
	  h2_jet_eta_phi[h][whichjet]->Fill(jet_et[h][k],jet_ph[h][k]);
	  h2_jet_eta_e[h][0]->Fill(jet_et[h][k],jet_e[h][k]);
	  if(ehjet[h][k] < maxeh && ehjet[h][k] > 0)
	    {
	      h2_jet_eta_e[h][1]->Fill(jet_et[h][k],jet_e[h][k]);
	      jetE[h][1]->Fill(ET);
	      h1_eta[h][1][whichjet]->Fill(jet_et[h][k]);
	      h1_phi[h][1][whichjet]->Fill(jet_ph[h][k]);
	    }
	  if(seedD[h][k] < 0.65)
	    {
	      h2_jet_eta_e[h][2]->Fill(jet_et[h][k],jet_e[h][k]);
	      jetE[h][2]->Fill(ET);
	      h1_eta[h][2][whichjet]->Fill(jet_et[h][k]);
	      h1_phi[h][2][whichjet]->Fill(jet_ph[h][k]);
	    }
	  if(seedD[h][k] < 0.65 && ehjet[h][k] < maxeh && ehjet[h][k] > 0)
	    {
	      h2_jet_eta_e[h][3]->Fill(jet_et[h][k],jet_e[h][k]);
	      h1_eta[h][3][whichjet]->Fill(jet_et[h][k]);
	      h1_phi[h][3][whichjet]->Fill(jet_ph[h][k]);
	      jetE[h][3]->Fill(ET);
	    }

	  for(int l = 0; l<njet[h]; ++l)
	    {
	      if(l==k) continue;
	      if(jet_e[h][k] < jet_e[h][l]) continue;
	      float ET1 = jet_e[h][k]/cosh(jet_et[h][k]);
	      float ET2 = jet_e[h][l]/cosh(jet_et[h][l]);
	      float AJ = (ET1-ET2)/(ET1+ET2);
	      float odph = jet_ph[h][k] - jet_ph[h][l];
	      if(odph < 0) odph += 2*M_PI;
	      if(odph > 3*M_PI/4 && odph < 5*M_PI/4)
		{
		  h1_AJ[h][0]->Fill(AJ);
		  if(seedD[h][k] < 0.65) h1_AJ[h][2]->Fill(AJ);
		  if(ehjet[h][k] < maxeh) h1_AJ[h][1]->Fill(AJ);
		  if(ehjet[h][k] < maxeh && seedD[h][k] < 0.65) h1_AJ[h][3]->Fill(AJ);
		
		  if(jet_e[h][k] > 7 && jet_e[h][l] > 4)
		    {
		      h1_AJ[h][8]->Fill(AJ);
		      if(seedD[h][k] < 0.65) h1_AJ[h][10]->Fill(AJ);
		      if(ehjet[h][k] < maxeh) h1_AJ[h][9]->Fill(AJ);
		      if(ehjet[h][k] < maxeh && seedD[h][k] < 0.65) h1_AJ[h][11]->Fill(AJ);
		    }
		}
	      if(jet_e[h][k] > 10 && jet_e[h][l] > 5)
		{
		  if(odph > 3*M_PI/4 && odph < 5*M_PI/4)
		    {
		      h1_AJ[h][4]->Fill(AJ);
		      if(seedD[h][k] < 0.65) h1_AJ[h][5]->Fill(AJ);
		      if(ehjet[h][k] < maxeh) h1_AJ[h][6]->Fill(AJ);
		      if(ehjet[h][k] < maxeh && seedD[h][k] < 0.65) h1_AJ[h][7]->Fill(AJ);
		    }
		  if(abs(jet_et[h][k]) > 0.7 || abs(jet_et[h][l]) > 0.7) continue;
		  float dphi = jet_ph[h][k] - jet_ph[h][l];
		  float deta = jet_et[h][k] - jet_et[h][l];
		  //if(sqrt(dphi*dphi+deta*deta) < 0.4 && seedD[h][k] < 0.65 && seedD[h][l] < 0.65 && ehjet[h][k] < maxeh && ehjet[h][k] > 0 && ehjet[h][l] > 0 && ehjet[h][l] < maxeh) dodraw = 1;
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
	  /*
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
	  e_stream << std::fixed << std::setprecision(1) << jet_e[h][k];
	  std::string e_string = e_stream.str();
	  drawText((e_string+" GeV").c_str(),get_eta(jet_et[h][k]),get_phi(jet_ph[h][k])+(get_phi(jet_ph[h][k])>55?-6:4.5),(get_eta(jet_et[h][k])>12?1:0),kBlack,0.07,42,false);

	  //d->cd();
	  //event_sum->Draw("LEGO");
	  */
	}
      if(leadjetET != 0) ljetE->Fill(leadjetET);
      if(leadjetEta != -2) ljetEta->Fill(leadjetEta);
      for(int k=0; k<4; ++k)
	{
	  h1_mlt[h][k]->Fill(njo[k]);
	}
      //if(seedD[h][k] < 0.65) h1_mlt[h][2]->Fill(njo5);
      //if(ehjet[h][k] < maxeh) h1_mlt[h][1]->Fill(njo5);
      //if(ehjet[h][k] < maxeh && seedD[h][k] < 0.65) h1_mlt[h][3]->Fill(njo5);
      /*
      c->cd(1);
      drawText(subdet[0].c_str(),0.75,0.9,1,kBlack,0.1);
      drawText(("#it{p}+#it{p} 200 GeV Run "+runnum).c_str(),0,0.97,0,kBlack,0.09);
      c->cd(2);
      drawText(subdet[1].c_str(),0.75,0.9,1,kBlack,0.1);
      drawText(("Event "+to_string(eventbase[h]+evtnum[h])).c_str(),0,0.97,0,kBlack,0.09);
      c->cd(3);
      drawText("2024/06/12",0.4,0.97,0,kBlack,0.1);
      drawText(subdet[2].c_str(),0.75,0.9,1,kBlack,0.09);
      c->cd(4);
      drawText(subdet[3].c_str(),0.75,0.9,1,kBlack,0.09);
      sphenixtext(1,0.97,1,0.09);
      //if(dispcount % 18 == 17 && highejet)
      if(dodraw)
	{
	  //cout << "Saving display " << dispcount/18 << endl;
	  //c->SaveAs(("./output/img/candidate_"+runnum+"_"+to_string(dispcount/18)+".pdf").c_str());
	  c->SaveAs(("./output/img/candidate_"+runnum+"_"+to_string(h)+"_"+to_string(eventbase[h]+evtnum[h])+".pdf").c_str());
	  //cout << "Saved!" << endl;
	  //c->Clear("D");
	}
      ++cancount;
      if(shighe)
	{
	  c->SaveAs(("./output/hmg/candidate_"+runnum+"_superhighE_"+to_string(h)+"_"+to_string(eventbase[h]+evtnum[h])+".pdf").c_str());
	  //d->SaveAs(("./output/hmg/candidate_"+runnum+"_superhighElego_"+to_string(eventbase[h]+evtnum[h])+".pdf").c_str());
	}
      if(sshighe) 
	{
	  
	  c->SaveAs(("./output/smg/candidate_"+runnum+"_supersuperhighE_"+to_string(h)+"_"+to_string(eventbase[h]+evtnum[h])+".pdf").c_str());
	  //d->SaveAs(("./output/smg/candidate_"+runnum+"_supersuperhighElego_"+to_string(eventbase[h]+evtnum[h])+".pdf").c_str());
	}
      if(highejet)
	{
	  //cout << dispcount << endl;
	  ++dispcount;
	}
      //cout << prevt << " " << evtnum[h] << endl;
      */
      prevt = evtnum[h];
      
    }
    }
  //if(dispcount % 18 != 0)
  //  {
  //    c->SaveAs(("./output/img/candidate_"+runnum+"_last.pdf").c_str());
  //  }
  d->cd();
  float erejf[9] = {14.24,42.76,111.95,281.76,681.85,1466.6,3033.03,6272.51,11086.03};
  float x[9] = {4,5,6,7,8,9,10,11,12};

  for(int h=nosim; h<2; ++h)
    {
      gPad->SetLogy(0);
      for(int i=0; i<3; ++i)
	{
	  //jetfrac[h][i]->Scale(1./(nevents[h]));
	}
      //h1_dphi->SetMarkerColor(kMagenta+3);
      for(int i=0; i<12; ++i)
	{
	  if(i<4)
	    {
	  h1_dphi[h][i]->GetYaxis()->SetLabelSize(0.04);
	  h1_dphi[h][i]->GetYaxis()->SetTitleSize(0.05);
	  h1_dphi[h][i]->GetYaxis()->SetTitleOffset(1.3);
	  h1_dphi[h][i]->GetXaxis()->SetLabelSize(0.04);
	  h1_dphi[h][i]->GetXaxis()->SetTitleSize(0.05);
	  h1_dphi[h][i]->GetXaxis()->SetTitleOffset(1);
	  //h1_dphi[h][i]->Scale(1./nevents[h]);
	  h1_dphi[h][i]->GetYaxis()->SetTitle("Counts");
	  h1_dphi[h][i]->GetXaxis()->SetTitle("#Delta#phi [rad]");
	  h1_dphi[h][i]->SetFillColorAlpha(kMagenta,0.4);
	  h1_dphi[h][i]->SetMarkerColor(kMagenta+1);
	  h1_dphi[h][i]->SetMarkerStyle(43);
	  h1_dphi[h][i]->SetMarkerSize(2.5);
	  /*
	  h1_dphi[h][i]->Draw("PE");
	  sphenixtext(0.96,0.96,1);
	  drawText("2024/06/12",0.1,0.97,0,kBlack,0.04);
	  drawText("#it{E}_{jet,1} > 10 GeV", 0.93, 0.9, 1, kBlack, 0.04);
	  drawText("#it{E}_{jet,2} > 5 GeV",0.93,0.85,1,kBlack,0.04);
	  drawText("Anti-#it{k}_{T} #it{R}=0.4",0.93,0.8,1,kBlack,0.04);
	  drawText("Calorimeter towers",0.93,0.75,1,kBlack,0.04);
	  drawText("at EM scale",0.93,0.7,1,kBlack,0.04);
	  drawText("|#it{#eta}^{jet}| < 0.7",0.93,0.65,1,kBlack,0.04);
	  if(i>1) drawText("E_{lead tower} < 0.65 E_{jet}", 0.93, 0.6, 1, kBlack, 0.04);
	  if(i==1 || i==3) drawText(("E_{EM} < "+to_string(maxeh)+" E_{HAD}").c_str(), 0.93, (i==1?0.6:0.55), 1, kBlack, 0.04);
	  d->SaveAs(("output/gmg/"+runnum+"_"+to_string(h)+"_"+to_string(cancount)+"_"+to_string(i)+"_dphi_1d.pdf").c_str());
	  */
	    }
	  h1_AJ[h][i]->GetYaxis()->SetLabelSize(0.02);
	  h1_AJ[h][i]->GetYaxis()->SetTitleSize(0.03);
	  h1_AJ[h][i]->GetYaxis()->SetTitleOffset(1.7);
	  h1_AJ[h][i]->GetXaxis()->SetLabelSize(0.02);
	  h1_AJ[h][i]->GetXaxis()->SetTitleSize(0.03);
	  //h1_AJ[h][i]->Scale(1./nevents[h]);
	  h1_AJ[h][i]->GetYaxis()->SetTitle("Counts");
	  h1_AJ[h][i]->GetXaxis()->SetTitle("A_{J}");
	  h1_AJ[h][i]->SetFillColorAlpha(kMagenta,0.4);
	  /*
	  h1_AJ[h][i]->Draw("HIST");
	  if(i>3 && i<7)drawText("Leading Jet E > 10 GeV", 0.5, 0.85, 0, kBlack, 0.02);
	  if(i>3 && i<7)drawText("Subleading Jet E > 5 GeV", 0.5, 0.825, 0, kBlack, 0.02);
	  if(i>7)drawText("Leading Jet E > 7 GeV", 0.5, 0.85, 0, kBlack, 0.02);
	  if(i>7)drawText("Subleading Jet E > 4 GeV", 0.5, 0.825, 0, kBlack, 0.02);
	  if(i%4>1) drawText("E_{lead tower} < 0.65 E_{jet}", 0.5, 0.8, 0, kBlack, 0.02);
	  if(i%4==1 || i%4==3) drawText(("E_{EM} < "+to_string(maxeh)+" E_{HAD}").c_str(), 0.5, (i==1?0.8:0.775), 0, kBlack, 0.02);
	  d->SaveAs(("output/gmg/"+runnum+"_"+to_string(h)+"_"+to_string(cancount)+"_"+to_string(i)+"_AJ_1d.pdf").c_str());
	  */
	}
      //h2_dphi_deta[h]->SetMarkerColor(kMagenta+3);
      h2_dphi_deta[h]->GetYaxis()->SetLabelSize(0.02);
      h2_dphi_deta[h]->GetYaxis()->SetTitleSize(0.03);
      h2_dphi_deta[h]->GetYaxis()->SetTitleOffset(1.7);
      h2_dphi_deta[h]->GetXaxis()->SetLabelSize(0.02);
      h2_dphi_deta[h]->GetXaxis()->SetTitleSize(0.03);
      //h2_dphi_deta[h]->Scale(1./nmb[h]);
      h2_dphi_deta[h]->GetYaxis()->SetTitle("#Delta#phi [rad]");
      h2_dphi_deta[h]->GetXaxis()->SetTitle("#Delta#eta [rad]");
      h2_dphi_deta[h]->GetZaxis()->SetTitle("Counts");
      //h2_dphi_deta[h]->Draw("COLZ");
      //d->SaveAs(("output/gmg/"+runnum+"_"+to_string(h)+"_"+to_string(cancount)+"_dphi_deta_2d.pdf").c_str());
      //jetfrac[h][0]->Draw();
      //d->SaveAs(("output/gmg/"+runnum+"_"+to_string(h)+"_"+to_string(cancount)+"_jetfrac0.pdf").c_str());
      //jetfrac[h][1]->Draw();
      //d->SaveAs(("output/gmg/"+runnum+"_"+to_string(h)+"_"+to_string(cancount)+"_jetfrac1.pdf").c_str());
      //jetfrac[h][2]->Draw();
      //d->SaveAs(("output/gmg/"+runnum+"_"+to_string(h)+"_"+to_string(cancount)+"_jetfrac2.pdf").c_str());
    }
  d->SetLogy();
  //float goodfrac = 1 - fakejets/tree->GetEntries();
  TGraph* erejg = new TGraph(9,x,erejf);
  TLegend* jetEleg = new TLegend(0.7,0.7,0.9,0.9);
  for(int i=0; i<4; ++i)
    {
      jetTrigE[i]->SetMarkerColor(kRed+i);
      jetTrigE[i]->SetMarkerSize(3);
      jetTrigE[i]->SetMarkerStyle(39+i);
    }
  jetEleg->AddEntry(jetE[0][0],"Minbias Pythia","p");
  jetEleg->AddEntry(jetE[1][0],"Minbias Data","p");
  jetEleg->AddEntry(jetTrigE[0],"Jet 4 GeV Trigger","p");
  jetEleg->AddEntry(jetTrigE[1],"Jet 6 GeV Trigger","p");
  jetEleg->AddEntry(jetTrigE[2],"Jet 8 GeV Trigger","p");
  jetEleg->AddEntry(jetTrigE[3],"Jet 10 GeV Trigger","p");
  for(int i=0; i<4; ++i)
    {
      
      //jetE[0][i]->Scale(1./(nmb[0]));
      //jetE[1][i]->Scale(1./(nmb[1]));
      //jetE[0][i]->SetMarkerColor(kGreen);
      //jetE[0][i]->GetYaxis()->SetRangeUser(0.00000003,1);
      //jetE[0][i]->Draw("P");
      /*
      jetE[1][i]->Draw("P");
      if(i==0)
	{
	  for(int j=0; j<4; ++j)
	    {
	      //jetTrigE[j]->Scale(1./nJetTrig[j]);
	      jetTrigE[j]->Draw("SAME P");
	    }
	}
      d->SaveAs(("output/gmg/"+runnum+"_"+to_string(cancount)+"_jetE"+to_string(i)+".pdf").c_str());
      */
      
      for(int h=nosim; h<2; ++h)
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
      //h1_rej[0][0]->Draw("P");
      //h1_rej[1][0]->Draw("P");
      /*
      for(int i=1; i<4; ++i)
	{
	  for(int h=nosim; h<2; ++h)
	    {
	      h1_rej[h][i]->Draw("SAME P");
	    }
	}
      */
      erejg->SetLineWidth(2);
      //erejg->Draw("SAME");
    }
  /*
  TLegend* leg = new TLegend(0.12,0.6,0.4,0.85);
  //leg->SetTextSize(1);
  for(int h=nosim; h<2; ++h)
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
  d->SaveAs(("output/gmg/"+runnum+"_"+to_string(cancount)+"_rej"+".pdf").c_str());
  */
  //hejet[h]->Scale(1./(nevents[h]));
  /*
  for(int h=nosim; h<2; ++h)
    {
      hejet[h]->Draw();
      d->SaveAs(("output/gmg/"+runnum+"_"+to_string(h)+"_"+to_string(cancount)+"_hejet.pdf").c_str());
    }
  */
  cout << "evtnum/ntree/nevents/nmb: " << cancount << " " << tree[1]->GetEntries() << " " << nevents[1]  << " " << nmb[1] << endl;
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
  for(int i=0; i<3; ++i)
    {
      outfile->WriteObject(h1_calo_occ[i],h1_calo_occ[i]->GetName());
      outfile->WriteObject(h1_calo_E[i],h1_calo_E[i]->GetName());
      outfile->WriteObject(h1_tower_E[i],h1_tower_E[i]->GetName());
      outfile->WriteObject(h2_occ_E[i],h2_occ_E[i]->GetName());
    }
  outfile->WriteObject(h2_tjet_eta_phi,h2_tjet_eta_phi->GetName());
  outfile->WriteObject(h1_cluster_E,h1_cluster_E->GetName());
  outfile->WriteObject(h1_cluster_eta,h1_cluster_eta->GetName());
  outfile->WriteObject(h2_cluster_eta_E,h2_cluster_eta_E->GetName());
  outfile->WriteObject(h2_cluster_eta_phi,h2_cluster_eta_phi->GetName());
  outfile->WriteObject(tjetE,tjetE->GetName());
  jetfile->WriteObject(ojt,ojt->GetName());
  outt->Fill();
  outfile->WriteObject(outt,outt->GetName());
  outfile->WriteObject(h1_zdist,h1_zdist->GetName());
  outfile->WriteObject(h1_mbdq,h1_mbdq->GetName());
  outfile->WriteObject(ljetE,ljetE->GetName());
  outfile->WriteObject(ljetEta,ljetEta->GetName());
  //outfile->WriteObject(jett[1],jett[1]->GetName());
  for(int h=nosim; h<2; ++h)
    {
      for(int i=0; i<12; ++i)
	{
	  if(i<3) outfile->WriteObject(jetfrac[h][i],jetfrac[h][i]->GetName());
	  if(i<3) outfile->WriteObject(h2_cal_eta_phi[h][i],h2_cal_eta_phi[h][i]->GetName());
	  if(i<5) outfile->WriteObject(h2_jet_eta_phi[h][i],h2_jet_eta_phi[h][i]->GetName());
	  if(i < 4)
	    {

	      outfile->WriteObject(h2_jet_eta_e[h][i],h2_jet_eta_e[h][i]->GetName());
	      outfile->WriteObject(jetE[h][i],jetE[h][i]->GetName());
	      outfile->WriteObject(h1_dphi[h][i],h1_dphi[h][i]->GetName());
	      outfile->WriteObject(h1_rej[h][i],h1_rej[h][i]->GetName());
	      if(h==1) outfile->WriteObject(jetTrigE[i],jetTrigE[i]->GetName());
	      for(int j=0; j<5; ++j)
		{
		  outfile->WriteObject(h1_eta[h][i][j],h1_eta[h][i][j]->GetName());
		  outfile->WriteObject(h1_phi[h][i][j],h1_phi[h][i][j]->GetName());
		}
	      outfile->WriteObject(h1_mlt[h][i],h1_mlt[h][i]->GetName());
	    }
	  outfile->WriteObject(h1_AJ[h][i],h1_AJ[h][i]->GetName());
	}
    }
  outfile->Write();
  return 0;
}
