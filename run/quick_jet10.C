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
#include <RooUnfold.h>
#include <RooUnfoldBayes.h>
#include <RooUnfoldResponse.h>
static const float radius_EM = 93.5;
static const float minz_EM = -130.23;
static const float maxz_EM = 130.23;

static const float radius_IH = 127.503;
static const float minz_IH = -170.299;
static const float maxz_IH = 170.299;

static const float radius_OH = 225.87;
static const float minz_OH = -301.683;
static const float maxz_OH = 301.683;

void jetMatch(int nt, int nr, float truthJet[][2], float recoJet[][2], int *whichMatch)
{
  const int nMaxJet = 100;
  float dRMat[nMaxJet][nMaxJet];
  for(int i=0; i<nMaxJet; ++i)
    {
      whichMatch[i] = -1;
    }
  for(int i=0; i<nt; ++i)
    {
      for(int j=0; j<nr; ++j)
	{
	  dRMat[i][j] = sqrt(pow(truthJet[i][0]-recoJet[j][0],2)+pow(truthJet[i][1]-recoJet[j][1],2));
	}
    }

  for(int i=0; i<nt; ++i)
    {
      float drMin = 10;
      int ri = -1;
      int ti = -1;
      for(int j=0; j<nt; ++j)
	{
	  if(whichMatch[j] != -1) continue;
	  for(int k=0; k<nr; ++k)
	    {
	      if(dRMat[j][k] < drMin && dRMat[j][k] < 0.3)
		{
		  ti = j;
		  ri = k;
		  drMin = dRMat[j][k];
		}
	    }
	}
      if(ti != -1) whichMatch[ti] = ri;
    }
}
	  
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

int quick_jet10(string filebase="", string samplestring="jet10", int njob=0, int dotow = 0)
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
  int ntj;
  float tjet_et[100];
  float tjet_phi[100];
  float tjet_eta[100];
  float frcoh[100];
  float dPhiLayer[100];
  float dPhi2pc[1000];
  float dEta2pc[1000];
  int n2pc;
  tree->SetBranchAddress("bbfqavec",&bbfqavec);
  tree->SetBranchAddress("vtx",vtx);
  tree->SetBranchAddress("njet",&njet);
  tree->SetBranchAddress("frcem",frcem);
  tree->SetBranchAddress("jet_ph",jet_ph);
  tree->SetBranchAddress("jet_e",jet_e);
  tree->SetBranchAddress("jet_et",jet_et);
  tree->SetBranchAddress("frcoh",frcoh);
  tree->SetBranchAddress("n2pc",&n2pc);

  tree->SetBranchAddress("ntj",&ntj);
  tree->SetBranchAddress("tjet_et",tjet_et);
  tree->SetBranchAddress("tjet_eta",tjet_eta);
  tree->SetBranchAddress("tjet_phi",tjet_phi);
  
  tree->SetBranchAddress("dPhiLayer",dPhiLayer);
  tree->SetBranchAddress("dPhi2pcd",dPhi2pc);
  tree->SetBranchAddress("dEta2pcd",dEta2pc);
  cout << "end getting branches" << endl;

  TFile* jetfile = TFile::Open(("output/simhists/run_"+samplestring+"_"+idstr+"_"+to_string(njob)+"_simhists.root").c_str(),"RECREATE");
  

  int eventbase = {0};
  int prevt = 0;
  const int nh2 = 6;
  TH2F* hists2[nh2];
  const int nz = 6;
  TH1F* h1_zdist[nz];
  const int maxJetToUse = 100;
  for(int i=0; i<nz; ++i)
    {
      h1_zdist[i] = new TH1F(("h1_zdist"+to_string(i)).c_str(),"",200,-150,150);
    }
  TH1F* h1_forrat[5];
  for(int i=0; i<5; ++i)
    {
      h1_forrat[i] = new TH1F(("h1_forrat_tosee"+to_string(i)).c_str(),"",1000,0,i<3?100:1);
    }

  float xlo[nh2] = {-0.2,-0.2,-0.2,-0.2,-0.2,-0.2};
  float xhi[nh2] = {1.2,1.2,1.2,1.2,1.2,1.2};
  float ylo[nh2] = {0,0,0,0,-0.2,-0.2};
  float yhi[nh2] = {100,100,100,100,1.2,1.2};
  
  for(int i=0; i<nh2; ++i)
    {
      hists2[i] = new TH2F(("hists2"+to_string(i)).c_str(),"",100,xlo[i],xhi[i],100,ylo[i],yhi[i]);
    }

  TH2F* h2_n2pc[2];
  TH2F* h2_dPhiLayer[2];

  for(int i=0; i<2; ++i)
    {
      h2_n2pc[i] = new TH2F(("dancheck_n2pc"+to_string(i)).c_str(),"",101,-1.2,1.2,257,-M_PI,M_PI);
      h2_dPhiLayer[i] = new TH2F(("dancheck_dPhiLayer"+to_string(i)).c_str(),"",80,-0.4,0.4,120,-0.1,1.1);
    }

  const int numTh2f = 36;
  TH2F* lJetEtaPhi[numTh2f];
  TH2F* lJetFrcemET[numTh2f];
  TH2F* lJetFrcemEta[numTh2f];
  TH2F* lJetFrcemPhi[numTh2f];
  TH2F* lJetEtaET[numTh2f];
  TH2F* lJetPhiET[numTh2f];

  for(int i=0; i<numTh2f; ++i)
    {
      lJetEtaPhi[i] = new TH2F(("lJetEtaPhi"+to_string(i)).c_str(),"",150,-1.5,1.5,128,-M_PI,M_PI);
      lJetFrcemET[i] = new TH2F(("lJetFrcemET"+to_string(i)).c_str(),"",150,-0.25,1.25,100,0,100);
      lJetFrcemEta[i] = new TH2F(("lJetFrcemEta"+to_string(i)).c_str(),"",150,-0.25,1.25,150,-1.5,1.5);
      lJetFrcemPhi[i] = new TH2F(("lJetFrcemPhi"+to_string(i)).c_str(),"",150,-0.25,1.25,128,-M_PI,M_PI);
      lJetEtaET[i] = new TH2F(("lJetEtaET"+to_string(i)).c_str(),"",150,-1.5,1.5,100,0,100);
      lJetPhiET[i] = new TH2F(("lJetPhiET"+to_string(i)).c_str(),"",128,-M_PI,M_PI,100,0,100);
    }
      
  
  TH1F* h1_g20_dijet = new TH1F("anotherdancheck","",60,-0.1,1.1);
  

  const int nbinx = 11;
  float binsx[nbinx+1] = {10,12,14,17,20,24,28,33,38,44,50,70};
  const int nbiny = 10;
  float binsy[nbiny+1] = {12,14,17,20,24,28,33,38,44,50,70};
  /*
    for(int i=-2; i<nbin-2; ++i)
    {
      bins[i+2] = 15*pow(6,((float)i)/(nbin-4));
    }
  */
  TH1F* h1_ucspec = new TH1F("h1_ucspec","ucspec",nbiny,binsy);
  TH1F* h1_cspec = new TH1F("h1_cspec","cspec",nbiny,binsy);
  TH1F* h1_tjetspec = new TH1F("h1_tjetspec","tjetspec",nbinx,binsx);
  TH2D* h2_resp = new TH2D("response_matrix","respmat",nbiny,binsy,nbinx,binsx);

  TH1F* h1_miss = new TH1F("h1_miss","miss",nbinx,binsx);
  TH1F* h1_fake = new TH1F("h1_fake","fake",nbiny,binsy);
  RooUnfoldResponse response(h1_cspec,h1_tjetspec,h2_resp);
  TH1F* xJ[2];
  xJ[0]=new TH1F("xJ0_sim","",100,0,1);
  xJ[1]=new TH1F("xJ1_sim","",100,0,1);
  int whichhist[3];
  float thresh[numTh2f/6] = {8,15,20,25,35,40};
  cout << "start processing: " << endl;
  for(int i=0; i<tree->GetEntries(); ++i)
    {

      tree->GetEntry(i);   
      

      if(abs(vtx[2]) > 30)
	{
	  continue;
	}
      float ljetET = 0;
      float subjetET = 0;
      float ljetph = 0;
      float subjetph = 0;
      float subjeteta = 0;
      float ljetfrcem = -1;
      int isdijet = 0;
      float ljeteta = 0;
      float ljetfrcoh = -1;
      float closejetdphi = M_PI;
      for(int j=0; j<njet; ++j)
	{
	  if(jet_e[j] > ljetET)
	    {
	      subjetET = ljetET;
	      subjetph = ljetph;
	      subjeteta = ljeteta;
	      ljetET = jet_e[j];
	      ljetph = jet_ph[j];
	      ljetfrcem = frcem[j];
	      ljeteta = jet_et[j];
	      ljetfrcoh = frcoh[j];
	    }
   	}
      for(int j=0; j<njet; ++j)
	{
	  float testdphi = ljetph - jet_ph[j];
	  if(testdphi < 0.05) continue;
	  if(check_bad_jet_eta(jet_et[j],vtx[2],0.4)) continue;
	  if(testdphi > M_PI) testdphi = 2*M_PI - testdphi;
	  if(testdphi < closejetdphi) closejetdphi = testdphi;
	}
      //cout << "got cut params" << endl;
      bool ljetHighEta = check_bad_jet_eta(ljeteta, vtx[2], 0.4);
      bool subjetHighEta = check_bad_jet_eta(subjeteta, vtx[2], 0.4);
      float truthJet[100][2];
      float recoJet[100][2];
      if(ljetET < 8 || ljetHighEta) continue;
      for(int j=0; j<njet; ++j)
	{
	  if(check_bad_jet_eta(jet_et[j],vtx[2],0.4)) continue;
	  h1_ucspec->Fill(jet_e[j]);
	  h2_dPhiLayer[0]->Fill(dPhiLayer[j],frcem[j]);
	}
      //cout << "ntj: " << ntj << endl;
      for(int j=0; j<ntj; ++j)
	{
	  if(check_bad_jet_eta(tjet_eta[j],vtx[2],0.4)) continue;
	  h1_tjetspec->Fill(tjet_et[j]);
	  truthJet[j][0] = tjet_eta[j];
	  truthJet[j][1] = tjet_phi[j];
	}



      
      if(subjetET > 8) isdijet = 1;
      if(subjetHighEta) isdijet = 0;
      float dphi = abs(ljetph - subjetph);
      if(dphi > M_PI) dphi = 2*M_PI - dphi;

      bool hdPhiCut = dphi < 3*M_PI/4 && isdijet;
      bool bbCut = (bbfqavec >> 5) & 1;
      //bool sdPhiCut = (ljetfrcem < 0.4 && dphi < 0.15) && isdijet;
      bool lETCut = ljetfrcem < 0.1 && (ljetET > (50*ljetfrcem+20)) && (hdPhiCut || !isdijet);
      bool hETCut = ljetfrcem > 0.9 && (ljetET > (-50*ljetfrcem+75)) && (hdPhiCut || !isdijet);
      bool ihCut = ljetfrcoh + ljetfrcem < 0.65;
      bool fullcut = bbCut || lETCut || hETCut || ihCut;
      h1_forrat[2]->Fill(ljetET);


      

      //cout << "before th2f filling" << endl;
      for(int j=0; j<numTh2f/6; ++j)
	{
	  if(ljetET < thresh[j]) continue;
	  if(fullcut) continue;
	  if(isdijet) whichhist[0] = j + (dphi>3*M_PI/4?numTh2f/6:0);
	  else whichhist[0] = -1;
	  if(isdijet) whichhist[1] = j + numTh2f/3 + (closejetdphi < M_PI/4? numTh2f/6:0);
	  else whichhist[1] = -1;
	  whichhist[2] = j + 2*numTh2f/3;
	  for(int k=0; k<3; ++k)
	    {
	      if(whichhist[k] < 0 || whichhist[k] > numTh2f - 1) continue;
	      lJetFrcemET[whichhist[k]]->Fill(ljetfrcem,ljetET);
	      lJetFrcemEta[whichhist[k]]->Fill(ljetfrcem,ljeteta);
	      lJetFrcemPhi[whichhist[k]]->Fill(ljetfrcem,ljetph);
	      lJetEtaET[whichhist[k]]->Fill(ljeteta,ljetET);
	      lJetPhiET[whichhist[k]]->Fill(ljetph,ljetET);
	      lJetEtaPhi[whichhist[k]]->Fill(ljeteta,ljetph);
	    }
	}
      //cout << "after th2f filling" << endl;
      if(ljetET > 20 && dphi > 7*M_PI/8)
	{
	  h1_g20_dijet->Fill(ljetfrcem);
	}

      if(isdijet && ljetET > 20 && subjetET > 10)
	{
	  h1_forrat[3]->Fill((ljetET-subjetET)/(ljetET+subjetET));
	  xJ[0]->Fill(subjetET/ljetET);
	}
      if(bbCut) h1_zdist[1]->Fill(vtx[2]);
      if(ihCut) h1_zdist[2]->Fill(vtx[2]);
      if(lETCut) h1_zdist[3]->Fill(vtx[2]);
      if(hETCut) h1_zdist[4]->Fill(vtx[2]);
      
      hists2[0]->Fill(ljetfrcem,ljetET);
      hists2[2]->Fill(ljetfrcoh,ljetET);
      hists2[4]->Fill(ljetfrcoh,ljetfrcem);
      //cout << "filled hists2" << endl;
      for(int j=0; j<n2pc; ++j)
	{
	  h2_n2pc[0]->Fill(dEta2pc[j],dPhi2pc[j]);
	}

      if(!fullcut)
	{
	  for(int j=0; j<n2pc; ++j)
	    {
	      h2_n2pc[1]->Fill(dEta2pc[j],dPhi2pc[j]);
	    }
	  hists2[1]->Fill(ljetfrcem,ljetET);
	  hists2[3]->Fill(ljetfrcoh,ljetET);
	  hists2[5]->Fill(ljetfrcoh,ljetfrcem);
	  h1_zdist[5]->Fill(vtx[2]);
	  h1_forrat[0]->Fill(ljetET);
	  //h1_cspec->Fill(ljetET);
	  for(int j=0; j<njet; ++j)
	    {
	      if(check_bad_jet_eta(jet_et[j],vtx[2],0.4)) continue;
	      if(jet_e[j] >8)
		{
		  recoJet[j][0] = jet_et[j];
		  recoJet[j][1] = jet_ph[j];
		  h1_cspec->Fill(jet_e[j]);
		  h2_dPhiLayer[1]->Fill(dPhiLayer[j],frcem[j]);
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
      int whichMatch[100];
      jetMatch(ntj,njet,truthJet,recoJet,whichMatch);
      for(int j=0; j<ntj; ++j)
	{
	  if(whichMatch[j] != -1)
	    {
	      response.Fill(jet_e[whichMatch[j]],tjet_et[j]);
	    }
	  else
	    {
	      response.Miss(tjet_et[j]);
	      h1_miss->Fill(tjet_et[j]);
	    }
	}
      bool isFake[100];
      for(int j=0; j<njet; ++j)
	{
	  isFake[j] = true;
	}
      for(int j=0; j<ntj; ++j)
	{
	  if(whichMatch[j] != -1)
	    {
	      isFake[whichMatch[j]] = false;
	    }
	}
      for(int j=0; j<njet; ++j)
	{
	  if(isFake[j]) response.Fake(jet_e[j]);
	  h1_fake->Fill(jet_e[j]);
	}
      h1_zdist[0]->Fill(vtx[2]); 
    }

  cout << "finished filling" << endl;
  h1_zdist[0]->Write();
  //h1_ucspec->Scale(4e-5/2.8e6);
  h1_ucspec->Write();
  //h1_cspec->Scale(4e-5/2.8e6);
  h1_cspec->Write();
  xJ[0]->Write();
  xJ[1]->Write();
  cout << "wrote spectra" << endl;
  for(int i=0; i<5; ++i) h1_forrat[i]->Write();
  cout << "writing hists2" << endl;
  for(int i=0; i<6; ++i)
    {
      hists2[i]->Write();
      if(i>0 && i<5)
	{
	  h1_zdist[i]->Write();
	}
    }
  cout << "writing truth spectra" << endl;
  //h1_tjetspec->Scale(4e-5/2.8e6);
  h1_tjetspec->Write();
  cout << "writing n2pc" << endl;
  for(int i=0; i<2; ++i)
    {
      h2_n2pc[i]->Write();
      h2_dPhiLayer[i]->Write();
    }
  cout << "writing th2f" << endl;
  for(int i=0; i<numTh2f; ++i)
    {
      lJetFrcemET[i]->Write();
      lJetFrcemEta[i]->Write();
      lJetFrcemPhi[i]->Write();
      lJetEtaET[i]->Write();
      lJetPhiET[i]->Write();
      lJetEtaPhi[i]->Write();
    }
  response.Hresponse()->Write("hresponse");
  jetfile->Write();
  return 0;
}
