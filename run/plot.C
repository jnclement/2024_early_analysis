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

void get_scaledowns(int runnumber, int scaledowns[])
{

  TSQLServer *db = TSQLServer::Connect("odbcn://daq","","");

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
      //cout << nrows << " " << nfields << endl;
      for (int i = 0; i < nrows; i++) {
	row = res->Next();
	for (int j = 0; j < nfields; j++) {
	  scaledowns[is] = stoi(row->GetField(j));
	  //cout << is << endl;
	  if(is == 10 || is==17 || is==18 || is==19) cout  << is << ":" << scaledowns[is] << " ";
	}
	delete row;
      }


      delete res;
    }
  cout << endl;
  delete db;
}


float min(float a, float b)
{
  if(a<b) return a;
  return b;
}
void std_text(TCanvas* thecan, string* texts, int* dotext, int ntext, float textsize, float textx, float texty, int rightalign, int nevtsim, int nevtdat, int fitgdat, float gmd, float ged, int fitgsim, float gms, float ges, int run = -1)
{
  std::stringstream stream[6];
  if(fitgdat)
    {
      stream[0] << std::fixed << std::setprecision(4) << gmd;
      stream[1] << std::fixed << std::setprecision(4) << ged;
    }
  if(fitgsim)
    {
      stream[2] << std::fixed << std::setprecision(4) << gms;
      stream[3] << std::fixed << std::setprecision(4) << ges;
    }
  stream[4] << std::setprecision(1) << std::scientific << (float)nevtsim;
  stream[5] << std::setprecision(1) << std::scientific << (float)nevtdat;
  thecan->cd();
  sphenixtext(0.9,0.91,1);
  drawText(("n_{MB evt,PYTHIA} = "+stream[4].str()+"  n_{MB evt,data} = "+stream[5].str()).c_str(),0.1,0.96,0,kBlack,0.03);
  //if(run == -1) drawText("Good runs 47002-47060",0.2,0.91,0,kBlack,0.03);
  //else drawText(("Run "+to_string(run)).c_str(),0.2,0.91,0,kBlack,0.03);
  float drawy = texty;
  for(int i=0; i<ntext; ++i)
    {
      if(dotext[i])
	{
	  drawText(texts[i].c_str(),textx,drawy,rightalign,kBlack,textsize);
	  drawy -= 5*textsize/4;
	}
    }
  if(fitgdat)
    {
      drawText(("#mu_{data} = "+stream[0].str()+" #sigma_{data} = "+stream[1].str()).c_str(),textx,drawy,rightalign,kBlack,textsize);
      drawy -= 5*textsize/4;
    }
  if(fitgsim)
    {
      drawText(("#mu_{sim} = "+stream[2].str()+" #sigma_{sim} = "+stream[3].str()).c_str(),textx,drawy,rightalign,kBlack,textsize);
    }
}

void std_hist(TH1D* dathist, TH1D* simhist)
{
  simhist->GetYaxis()->SetTitle("Event Normalized Counts");
  simhist->SetMarkerSize(2);
  dathist->SetMarkerSize(2);
  simhist->SetMarkerStyle(25);
  simhist->SetMarkerColor(kGreen+2);
  dathist->SetMarkerStyle(20);
  simhist->SetLineWidth(2);
  dathist->SetLineWidth(2);
  dathist->SetMarkerColor(kMagenta+1);
  dathist->GetXaxis()->SetTitleSize(0.03);
  dathist->GetXaxis()->SetLabelSize(0.02);
  dathist->GetYaxis()->SetTitleSize(0.03);
  dathist->GetYaxis()->SetLabelSize(0.02);
  dathist->GetYaxis()->SetTitleOffset(1.7);
  simhist->GetXaxis()->SetTitleSize(0.03);
  simhist->GetXaxis()->SetLabelSize(0.02);
  simhist->GetYaxis()->SetTitleSize(0.03);
  simhist->GetYaxis()->SetLabelSize(0.02);
  simhist->GetYaxis()->SetTitleOffset(1.7);
}

float max(float a, float b)
{
  if(a>b) return a;
  return b;
}

int plot(string filebase = "summed_dat.root", string sfilebase="summed_sim.root", string filelist="lists/sumdatlist.list",int therun = -1)
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
  const int nruns = 286;//33;//72//10;//94
  const int nsd = 64;
  int sds[nruns][nsd] = {0};
  /*
  sds[0][1] = 1;
  sds[1][1] = 1;
  sds[2][1] = 1;
  sds[0][0] = 100;
  sds[1][0] = 100;
  sds[2][0] = 100;
  sds[3][0] = 100;
  sds[4][0] = 100;
  sds[5][0] = 100;
  sds[6][0] = 100;
  sds[7][0] = 100;
  sds[8][0] = 25;
  sds[9][0] = 25;
  sds[10][0] = 25;
  sds[11][0] = 25;
  sds[12][0] = 25;
  sds[13][0] = 50;
  sds[14][0] = 50;
  sds[14][0] = 30;
  sds[15][0] = 22;
  sds[16][0] = 17;
  sds[17][0] = 17;
  sds[21][0] = 70;
  sds[22][0] = 70;
  sds[23][0] = 70;
  sds[24][0] = 27;
  sds[25][0] = 27;
  sds[26][0] = 27;

  sds[0][0] = 20554;
  sds[1][0] = 4162;
  sds[2][0] = 590;
  sds[3][0] = 438;
  sds[4][0] = 5062;
  sds[5][0] = 606;
  sds[7][0] = 454;
  sds[8][0] = 280;
  sds[9][0] = 250;
  */
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0);
  TFile* rootfile[3];
  rootfile[0] = TFile::Open(sfilebase.c_str());
  TTree* simt = (TTree*)rootfile[0]->Get("st");
  int nsimmb = 0;
  int simmb = 0;
  int simevt = 0;
  simt->SetBranchAddress("outnmb",&simmb);
  simt->SetBranchAddress("outevt",&simevt);
  for(int i=0; i<simt->GetEntries(); i++)
    {
      simt->GetEntry(i);
      nsimmb += simmb;
    }
  int runnumbers[nruns] = {47002,47003,47004,47005,47006,47007,47009,47010,47013,47014,47015,47016,47017,47018,47019,47020,47021,47022,47032,47033,47034,47035,47036,47037,47038,47039,47040,47041,47042,47043,47050,47051,47052,47053,47054,47055,47056,47057,47058,47059,47060,47061,47064,47066,47068,47079,47080,47081,47082,47086,47089,47090,47098,47099,47100,47101,47102,47104,47105,47106,47114,47115,47116,47117,47118,47119,47120,47121,47122,47123,47124,47125,47126,47127,47128,47129,47130,47131,47132,47133,47135,47136,47137,47138,47139,47140,47141,47142,47143,47153,47154,47155,47156,47157,47158,47160,47161,47162,47182,47187,47197,47198,47201,47202,47203,47204,47208,47211,47216,47217,47218,47219,47220,47222,47229,47230,47286,47287,47288,47289,47291,47293,47296,47297,47298,47300,47301,47303,47305,47306,47308,47309,47310,47311,47312,47313,47314,47315,47316,47323,47325,47329,47330,47332,47333,47334,47357,47359,47360,47361,47362,47363,47375,47376,47377,47378,47379,47380,47381,47382,47389,47391,47392,47393,47394,47395,47396,47397,47398,47399,47403,47406,47439,47440,47441,47442,47443,47444,47445,47450,47451,47452,47453,47454,47455,47456,47457,47458,47459,47460,47461,47462,47463,47464,47474,47475,47476,47477,47478,47479,47480,47481,47482,47483,47484,47485,47486,47487,47488,47489,47490,47491,47492,47494,47495,47496,47497,47501,47502,47503,47505,47506,47507,47508,47509,47513,47514,47515,47516,47517,47519,47520,47521,47522,47523,47524,47525,47531,47532,47538,47539,47540,47547,47548,47549,47550,47551,47552,47557,47558,47559,47566,47567,47568,47570,47611,47617,47634,47635,47636,47637,47638,47647,47649,47650,47651,47652,47653,47654,47656,47657,47658,47659,47660,47661,47662,47663,47664,47665,47666,47667,47698,47715,47722,47724,47725}; //{47474,47476,47480,47481,47484,47485,47491,47492,47494,47495,47503,47505,47506,47507,47513,47514,47516,47522,47524,47525,47540,47548,47552,47557,47568,47634,47636,47638,47657,47658,47659,47661,47662}; //{47002,47006,47013,47017,47021,47034,47038,47042,47052,47056,47060,47068,47082,47098,47102,47114,47118,47122,47126,47130,47135,47139,47143,47156,47161,47197,47203,47216,47220,47286,47291,47298,47305,47310,47314,47325,47333,47360,47375,47379,47389,47394,47398,47439,47443,47451,47455,47459,47463,47476,47480,47484,47488,47492,47497,47505,47509,47516,47521,47525,47539,47549,47557,47567,47617,47637,47650,47654,47659,47663,47667,47724}; //{47002,47005,47007,47009,47017,47019,47032,47037,47056,47060};//{47289,47293,47297,47298,47303,47305,47306,47308,47309,47310,47315,47316,47323,47325,47330,47334,47360,47377,47378,47379,47381,47382,47391,47393,47394,47395,47396,47397,47398,47399,47443,47445,47451,47455,47456,47457,47458,47459,47460,47464,47474,47480,47481,47484,47485,47491,47492,47494,47495,47497,47502,47503,47505,47506,47507,47508,47513,47514,47516,47522,47524,47525,47540,47548,47552,47557,47558,47568,47634,47636,47657,47658,47659,47661,47662,47666,47667,/*47698,*/47715,47716,47717,47718,47720,47722,47723,47724,47725,47726,47727,47728,47729,47730,47731,47732,47733};//,47741,47745,47746,47747,47748,47749,47751,47752,47753,47757,47758};
  long long unsigned int nevents = 0;
  int evtct = 0;
  int mbevt = 0;
  int njmb = 0;
  int jmb = 0;
  int maxeh = 5;
  long long unsigned int nmb = 0;
  int nJetTrig[4] = {0};
  float effevt[3] = {0};
  int whichrun = 0;
  int runmb = 0;
  string line;
  ifstream file;
  file.open(filelist);
  int nfile = 0;
  for(int i=0; i<nruns; ++i)
    {
      get_scaledowns(runnumbers[i],sds[i]);
    }
  while(getline(file,line))
    {
      if(therun != -1 && nfile != therun)
	{
	  ++whichrun;
	  ++nfile;
	  continue;
	}

      rootfile[2] = TFile::Open(line.c_str());
      cout << "filename: " << line << " address: " << rootfile[2] << endl;
      TTree* tree = (TTree*)rootfile[2]->Get("outt");
      tree->SetBranchAddress("outevt",&evtct);
      tree->SetBranchAddress("outnmb",&mbevt);
      tree->SetBranchAddress("njmb",&njmb);
      tree->SetBranchAddress("nJetTrig",nJetTrig);

      for(int i=0; i<tree->GetEntries(); ++i)
	{
	  tree->GetEntry(i);
	  nevents += evtct;
	  nmb += (sds[whichrun][10]==-1?0:1)*mbevt;
	  runmb += (sds[whichrun][10]==-1?0:1)*mbevt;
	  jmb += njmb;
	}
      cout << "runmb: "<< runmb << endl;
      if(sds[whichrun][17] != -1) effevt[0] += (1.*(1+sds[whichrun][10]))/(1.*(1+sds[whichrun][17]))*runmb;//(nJetTrig[1]+mbevt);
      if(sds[whichrun][18] != -1) effevt[1] += (1.*(1+sds[whichrun][10]))/(1.*(1+sds[whichrun][18]))*runmb;//(nJetTrig[2]+mbevt);
      if(sds[whichrun][19] != -1) effevt[2] += (1.*(1+sds[whichrun][10]))/(1.*(1+sds[whichrun][19]))*runmb;//(nJetTrig[3]+mbevt);      
      runmb = 0;
      whichrun++;
      ++nfile;
    }
  if(therun == -1)
    {
      cout << "thegrombo" << endl;
      rootfile[1] = TFile::Open(filebase.c_str());
    }
  else
    {
      rootfile[1] = rootfile[2];
      cout << "Address: " << rootfile[1] << endl;
    }
  //cout << whichrun << endl;
  //effevt[0];
  //effevt[1];
  //effevt[2];
  cout << nmb  << " " << nevents << " " << nsimmb << " " << effevt[0] << " " << nJetTrig[1] << endl;
  
  TH1D* h1_dphi[2][4];
  TH1D* h1_rej[2][4];
  TF1* ajgaus[2][12];
  TF1* dphigaus[2][4];
  TF1* zgaus[2];
  TH1D* h1_mlt[2][4];
  TH1D* h1_phi[2][4][5];
  TH1D* h1_eta[2][4][5];
  TH1D* jetfrac[2][4];
  TH1D* jetTrigE[4];
  TH2D* h2_jet_eta_phi[2][5];
  TH2D* h2_tjet_eta_phi;
  TH2D* h2_jet_eta_e[2][4];
  TH2D* h2_cal_eta_phi[2][3];
  TH1D* h1_zdist[2];
  TH1D* h1_mbdq[2];
  TH1D* ljetE[2];
  TH1D* ljetEta[2];
  TH1D* tjetE;
  TH1D* h1_tower_E[2][3];
  TH1D* h1_calo_E[2][3];
  TH1D* h1_calo_occ[2][3];
  TH2D* h2_occ_E[2][3];
  TH1D* h1_cluster_E[2];
  TH1D* h1_cluster_eta[2];
  TH2D* h2_cluster_eta_E[2];
  TH2D* h2_cluster_eta_phi[2];
  tjetE = (TH1D*)rootfile[0]->Get("tjetE");
  tjetE->Scale(1./nsimmb);
  for(int h=0; h<2; ++h)
    {
      h1_cluster_E[h] = (TH1D*)rootfile[h]->Get("h1_cluster_E");
      h1_cluster_eta[h] = (TH1D*)rootfile[h]->Get("h1_cluster_eta");
      h2_cluster_eta_E[h] = (TH2D*)rootfile[h]->Get("h2_cluster_eta_E");
      h2_cluster_eta_phi[h] = (TH2D*)rootfile[h]->Get("h2_cluster_eta_phi");
      h1_cluster_E[h]->Scale(1./(h==0?nsimmb:nmb));
      h1_cluster_eta[h]->Scale(1./(h==0?nsimmb:nmb));
      h2_cluster_eta_E[h]->Scale(1./(h==0?nsimmb:nmb));
      h2_cluster_eta_phi[h]->Scale(1./(h==0?nsimmb:nmb));
      h1_zdist[h] = (TH1D*)rootfile[h]->Get("h1_zdist");
      h1_zdist[h]->Scale(1./(h==0?nsimmb:nmb));
      cout << h1_zdist[h]->GetBinContent(0) << " " << h1_zdist[h]->GetBinContent(201) << " " << h1_zdist[h]->GetBinContent(2) << endl;
      h1_zdist[h]->Fit("gaus","QS");
      zgaus[h] = h1_zdist[h]->GetFunction("gaus");
      h1_mbdq[h] = (TH1D*)rootfile[h]->Get("h1_mbdq");
      ljetE[h] = (TH1D*)rootfile[h]->Get("ljetE");
      ljetEta[h] = (TH1D*)rootfile[h]->Get("ljetEta");
      for(int i=0; i<4; ++i)
	{
	  if(i<3)
	    {
	      h1_tower_E[h][i] = (TH1D*)rootfile[h]->Get(("h1_tower_E"+to_string(i)).c_str());
	      h1_tower_E[h][i]->Scale(1./(h==0?nsimmb:nmb));
	      h1_calo_E[h][i] = (TH1D*)rootfile[h]->Get(("h1_calo_E"+to_string(i)).c_str());
	      h1_calo_E[h][i]->Scale(1./(h==0?nsimmb:nmb));
	      h1_calo_occ[h][i] = (TH1D*)rootfile[h]->Get(("h1_calo_occ"+to_string(i)).c_str());
	      h1_calo_occ[h][i]->Scale(1./(h==0?nsimmb:nmb));
	      h2_occ_E[h][i] = (TH2D*)rootfile[h]->Get(("h2_occ_E"+to_string(i)).c_str());
	      h2_occ_E[h][i]->Scale(1./(h==0?nsimmb:nmb));
	    }
	  if(h==1 && i > 0)
	    {
	      jetTrigE[i] = (TH1D*)rootfile[h]->Get(("jetTrigE"+to_string(i)).c_str());
	      //cout << jetTrigE[i]->GetBinContent(100) << endl;
	      cout << "Jet trig E integral/effevt: " <<jetTrigE[i]->Integral(130,400) << " " << effevt[i-1] << endl;
	      jetTrigE[i]->Scale(1./effevt[i-1]);

	    }
	  
	  h1_dphi[h][i] = (TH1D*)rootfile[h]->Get(("h1_dphi"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_dphi[h][i]->Scale(1./(h==0?nsimmb:nmb));
	  h1_dphi[h][i]->Fit("gaus","QS","",1.,5.);
	  dphigaus[h][i] = h1_dphi[h][i]->GetFunction("gaus");
	  h1_rej[h][i] = (TH1D*)rootfile[h]->Get(("h1_rej"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_rej[h][i]->Rebin(5);
	  h1_mlt[h][i] = (TH1D*)rootfile[h]->Get(("h1_mlt"+to_string(1)+"_"+to_string(i)).c_str());
	  
	  h2_jet_eta_e[h][i] = (TH2D*)rootfile[h]->Get(("h2_jet_eta_e"+to_string(1)+to_string(i)).c_str());
	  h2_jet_eta_e[h][i]->Scale(1./(h==0?nsimmb:nmb));
	  if(h==0) h2_tjet_eta_phi = (TH2D*)rootfile[h]->Get("h2_tjet_eta_phi");
	  if(i<3)
	    {
	      h2_cal_eta_phi[h][i] = (TH2D*)rootfile[h]->Get(("h2_cal_eta_phi"+to_string(1)+to_string(i)).c_str());
	      //h2_cal_eta_phi[h][i]->Rebin2D(4,8);
	      h2_cal_eta_phi[h][i]->Scale(1./(h==0?nsimmb:nmb));
	    }
	  for(int j=0; j<5; ++j)
	    {
	      if(i==0) 
		{
		  h2_jet_eta_phi[h][j] = (TH2D*)rootfile[h]->Get(("h2_jet_eta_phi"+to_string(1)+to_string(j)).c_str());
		  //h2_jet_eta_phi[h][j]->Rebin2D(2,4);
		  h2_jet_eta_phi[h][j]->Scale(1./(h==0?nsimmb:nmb));
		}
	      h1_phi[h][i][j] = (TH1D*)rootfile[h]->Get(("h1_phi"+to_string(1)+"_"+to_string(i)+"_"+to_string(j)).c_str());
	      h1_eta[h][i][j] = (TH1D*)rootfile[h]->Get(("h1_eta"+to_string(1)+"_"+to_string(i)+"_"+to_string(j)).c_str());
	    }
	  //cout <<endl << h1_phi[h][i][0]->GetEntries() <<" " << h1_eta[h][i][0]->GetEntries() << endl << endl;
	  jetfrac[h][i] = (TH1D*)rootfile[h]->Get(("jetfrac"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_mlt[h][i]->Scale(1./(h==0?(nsimmb):nmb));
	  for(int j=0; j<5; ++j)
	    {
	      h1_phi[h][i][j]->Scale(1./(h1_phi[h][i][j]->Integral()));
	      h1_eta[h][i][j]->Scale(1./(h1_eta[h][i][j]->Integral()));
	    }
	  if(i<3) jetfrac[h][i]->Scale(1./(h==0?(nsimmb):nmb));
	}
    }
  TH1D* h1_AJ[2][12];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  h1_AJ[h][i] = (TH1D*)rootfile[h]->Get(("h1_AJ"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_AJ[h][i]->Scale(1./h1_AJ[h][i]->Integral());//1./(h==0?(tree->GetEntries()*nsimmb):nmb));
	  h1_AJ[h][i+4] = (TH1D*)rootfile[h]->Get(("h1_AJ"+to_string(1)+"_"+to_string(i+4)).c_str());
	  //cout << h1_AJ[h][i+4] << endl;
	  h1_AJ[h][i+4]->Scale(1./h1_AJ[h][i+4]->Integral());//1./(h==0?(tree->GetEntries()*nsimmb):nmb));
	  h1_AJ[h][i]->Fit("gaus","QS");
	  ajgaus[h][i] = h1_AJ[h][i]->GetFunction("gaus");
	  h1_AJ[h][i+4]->Fit("gaus","QS");
	  ajgaus[h][i+4] = h1_AJ[h][i+4]->GetFunction("gaus");
	  h1_AJ[h][i+8] = (TH1D*)rootfile[h]->Get(("h1_AJ"+to_string(1)+"_"+to_string(i+8)).c_str());
	  h1_AJ[h][i+8]->Scale(1./h1_AJ[h][i+8]->Integral());//1./(h==0?(tree->GetEntries()*nsimmb):nmb));
	  h1_AJ[h][i+8]->Fit("gaus","QS");
	  ajgaus[h][i+8] = h1_AJ[h][i+8]->GetFunction("gaus");
	}
    }
  TH1D* jetE[2][4];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  jetE[h][i] = (TH1D*)rootfile[h]->Get(("jetE"+to_string(1)+"_"+to_string(i)).c_str());
	  //cout << "jetE" << h << " " << i << " nentries: " << jetE[h][i]->GetEntries() << endl;
	  cout << "jetE integral/effevt "<< jetE[h][i]->Integral(130,400) << " " << (h==0?(nsimmb):nmb) << endl;
	  jetE[h][i]->Scale(1./(h==0?(nsimmb):nmb));
	  jetE[h][i]->GetXaxis()->SetTitle("E_{T,jet} [GeV]");
	}
      ljetE[h]->Scale(1./(h==0?nsimmb:nmb));
      ljetEta[h]->Scale(1./ljetEta[h]->Integral());
    }
  TCanvas* d = new TCanvas("","",1000,1000);
  //  TH1D* jetTrigE[4];
  /*
  for(int i=0; i<4; i++)
    {
      jetTrigE[i] = (TH1D*)rootfile[h]->Get(("jetTrigE"+to_string(i)).c_str());
    }
  */
  cout << "Passed getting" << endl;
  const int ntext = 14;
  string texts[ntext] =
    {
      "Anti-#it{k}_{T} #it{R}=0.4",
      "#it{E}_{T,jet,1} > 10 GeV",
      "#it{E}_{T,jet,2} > 5 GeV",
      "|#it{#eta}^{jet}| < 0.7",
      "#it{E}_{lead tower} < 0.65 #it{E}_{jet}",
      "#it{E}_{EM} < "+to_string(maxeh)+" #it{E}_{HAD}",
      "#it{E}_{T,jet,2} > 7 GeV",
      "#it{E}_{T,jet,2} > 4 GeV",
      "#it{E}_{T,jet} > 4 GeV",
      "4 GeV < #it{E}_{T,jet} < 6 GeV",
      "6 GeV < #it{E}_{T,jet} < 10 GeV",
      "10 GeV < #it{E}_{T,jet} < 15 GeV",
      "15 GeV < #it{E}_{T,jet} < 20 GeV",
      "#it{E}_{T,jet} > 20 GeV"
    };

  int dotext[ntext] = {0};
  float stdsize = 0.025;
  float stdx = 0.87;
  float stdy = 0.88-stdsize;
  int stdright = 1;
  int nevtsim = nsimmb;
  int nevtdat = nmb;
  int fitgdat = 1;
  int fitgsim = 0;
  float gmd, ged, gms, ges;
      
  for(int h=0; h<2; ++h)
    {
      gPad->SetLogy(0);
      for(int i=0; i<4; ++i)
	{
	  h1_dphi[h][i]->SetMarkerSize(3);
	  h1_dphi[h][i]->Draw("PE");
	  if(dphigaus[h][i]) dphigaus[h][i]->SetLineColor(kMagenta+1);
	  if(i>1) dotext[4] = 1;
	  else dotext[4] = 0;
	  if(i==1 || i==3) dotext[5] = 1;
	  else dotext[5] = 0;
	  dotext[3] = 1;
	  if(dphigaus[h][i])
	    {
	      gmd = dphigaus[h][i]->GetParameter(1);
	      ged = dphigaus[h][i]->GetParameter(2);
	    }
	  std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
	  d->SaveAs(("output/rmg/summed_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(h)+"_"+to_string(i)+"_dphi_1d.png").c_str());
	}
    }
  TLegend* ajleg = new TLegend(0,0,0.4,0.08);
  ajleg->SetFillStyle(0);
  ajleg->SetFillColor(0);
  ajleg->SetTextFont(42);
  ajleg->SetBorderSize(0);
  h1_AJ[0][0]->SetMarkerSize(2);
  h1_AJ[1][0]->SetMarkerSize(2);
  h1_AJ[0][0]->SetMarkerStyle(25);
  h1_AJ[0][0]->SetMarkerColor(kGreen+2);
  h1_AJ[1][0]->SetMarkerStyle(20);
  h1_AJ[1][0]->SetMarkerColor(kMagenta+1);
  ajleg->AddEntry(h1_AJ[0][0],"MBD N/S>=1 PYTHIA","p");
  ajleg->AddEntry(h1_AJ[1][0],"MBD N/S>=1 Data","p");

  for(int i=0; i<3; ++i)
    {
      std_hist(h1_calo_E[1][i],h1_calo_E[0][i]);
      std_hist(h1_calo_occ[1][i],h1_calo_occ[0][i]);
      std_hist(h1_tower_E[1][i],h1_tower_E[0][i]);
      //std_hist(h2_occ_E[1][i],h2_occ_E[0][i]);
      string calo;
      switch (i)
	{
	case 0:
	  calo = "EMCal";
	  break;
	case 1:
	  calo = "IHCal";
	  break;
	case 2:
	  calo = "OHCal";
	  break;
	default:
	  calo = "Error";
	}
      gPad->SetLogy();
      h1_calo_E[0][i]->GetXaxis()->SetTitle((calo+" Total E [GeV]").c_str());
      h1_calo_E[0][i]->GetYaxis()->SetTitle("Event normalized counts");
      h1_calo_E[0][i]->Draw("PE");
      h1_calo_E[1][i]->Draw("SAME P E");
      ajleg->Draw();
      gPad->SaveAs(("output/rmg/totalE"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(i)+".png").c_str());
      //gPad->SetLogy();
      h1_tower_E[0][i]->GetXaxis()->SetTitle((calo+" Tower E [GeV]").c_str());
      h1_tower_E[0][i]->GetYaxis()->SetTitle("Event normalized counts");
      h1_tower_E[0][i]->Draw("PE");
      h1_tower_E[1][i]->Draw("SAME P E");
      ajleg->Draw();
      gPad->SaveAs(("output/rmg/towerE"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(i)+".png").c_str());
      h1_calo_occ[0][i]->GetXaxis()->SetTitle((calo+" Occupancy").c_str());
      h1_calo_occ[0][i]->GetYaxis()->SetTitle("Event normalized counts");
      h1_calo_occ[0][i]->Draw("PE");
      h1_calo_occ[1][i]->Draw("SAME P E");
      ajleg->Draw();
      gPad->SaveAs(("output/rmg/caloocc"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(i)+".png").c_str());
      gPad->SetLogy(0);
      gPad->SetRightMargin(0.2);
      h2_occ_E[0][i]->GetYaxis()->SetTitle((calo+" Total E [GeV]").c_str());
      h2_occ_E[0][i]->GetXaxis()->SetTitle((calo+" Occupancy").c_str());
      h2_occ_E[0][i]->GetZaxis()->SetTitle("Event Normalized Counts");
      h2_occ_E[0][i]->Draw("COLZ");
      gPad->SaveAs(("output/rmg/occE0_"+to_string((therun==-1?-1:runnumbers[therun]))+"_0_"+to_string(i)+".png").c_str());

      h2_occ_E[1][i]->GetYaxis()->SetTitle((calo+" Total E [GeV]").c_str());
      h2_occ_E[1][i]->GetXaxis()->SetTitle((calo+" Occupancy").c_str());
      h2_occ_E[1][i]->GetZaxis()->SetTitle("Event Normalized Counts");
      h2_occ_E[1][i]->Draw("COLZ");
      gPad->SaveAs(("output/rmg/occE0_"+to_string((therun==-1?-1:runnumbers[therun]))+"_1_"+to_string(i)+".png").c_str());
      gPad->SetRightMargin(0.1);
    }
  
  std_hist(h1_cluster_E[1],h1_cluster_E[0]);
  std_hist(h1_cluster_eta[1],h1_cluster_eta[0]);
  h1_cluster_eta[0]->GetYaxis()->SetRangeUser(0,1.1*max(h1_cluster_eta[1]->GetMaximum(),h1_cluster_eta[0]->GetMaximum()));
  h1_cluster_E[0]->GetXaxis()->SetTitle("Cluster E [GeV]");
  h1_cluster_eta[0]->GetXaxis()->SetTitle("Cluster #eta");
  h1_cluster_E[0]->GetYaxis()->SetTitle("Event Normalized Counts");
  h1_cluster_eta[0]->GetYaxis()->SetTitle("Event Normalized Counts");
  h1_cluster_E[0]->GetYaxis()->SetRangeUser(0.9*min(h1_cluster_E[1]->GetBinContent(h1_cluster_E[1]->FindLastBinAbove()),h1_cluster_E[0]->GetBinContent(h1_cluster_E[0]->FindLastBinAbove())),1.1*max(h1_cluster_E[1]->GetMaximum(),h1_cluster_E[0]->GetMaximum()));
  gPad->SetLogy();
  h1_cluster_E[0]->Draw("PE");
  h1_cluster_E[1]->Draw("SAME PE");
  ajleg->Draw();
  gPad->SaveAs(("output/rmg/clusterE"+to_string((therun==-1?-1:runnumbers[therun]))+".png").c_str());
  gPad->SetLogy(0);
  h1_cluster_eta[0]->Draw("PE");
  h1_cluster_eta[1]->Draw("SAME PE");
  ajleg->Draw();
  gPad->SaveAs(("output/rmg/clustereta"+to_string((therun==-1?-1:runnumbers[therun]))+".png").c_str());
  gPad->SetRightMargin(0.2);
  h2_cluster_eta_E[0]->GetXaxis()->SetTitle("Cluster #eta");
  h2_cluster_eta_E[0]->GetYaxis()->SetTitle("Cluster E [GeV]");
  h2_cluster_eta_E[0]->GetZaxis()->SetTitle("Event Normalized Counts");
  h2_cluster_eta_phi[0]->GetXaxis()->SetTitle("Cluster #eta");
  h2_cluster_eta_phi[0]->GetYaxis()->SetTitle("Cluster #phi [rad]");
  h2_cluster_eta_phi[0]->GetZaxis()->SetTitle("Event Normalized Counts");
  h2_cluster_eta_E[0]->Draw("COLZ");
  gPad->SaveAs(("output/rmg/clusteretaE"+to_string((therun==-1?-1:runnumbers[therun]))+"_0.png").c_str());
  h2_cluster_eta_phi[0]->Draw("COLZ");
  gPad->SaveAs(("output/rmg/clusteretaphi"+to_string((therun==-1?-1:runnumbers[therun]))+"_0.png").c_str());
  h2_cluster_eta_E[1]->GetXaxis()->SetTitle("Cluster #eta");
  h2_cluster_eta_E[1]->GetYaxis()->SetTitle("Cluster E [GeV]");
  h2_cluster_eta_E[1]->GetZaxis()->SetTitle("Event Normalized Counts");
  h2_cluster_eta_phi[1]->GetXaxis()->SetTitle("Cluster #eta");
  h2_cluster_eta_phi[1]->GetYaxis()->SetTitle("Cluster #phi [rad]");
  h2_cluster_eta_phi[1]->GetZaxis()->SetTitle("Event Normalized Counts");
  h2_cluster_eta_E[1]->Draw("COLZ");
  gPad->SaveAs(("output/rmg/clusteretaE"+to_string((therun==-1?-1:runnumbers[therun]))+"_1.png").c_str());
  h2_cluster_eta_phi[1]->Draw("COLZ");
  gPad->SaveAs(("output/rmg/clusteretaphi"+to_string((therun==-1?-1:runnumbers[therun]))+"_1.png").c_str());
  gPad->SetRightMargin(0.1);
  for(int i=0; i<12; ++i)
    {
      h1_AJ[0][i]->GetYaxis()->SetTitle("Integral Normalized Counts");
      h1_AJ[0][i]->SetMarkerSize(2);
      h1_AJ[1][i]->SetMarkerSize(2);
      h1_AJ[0][i]->SetMarkerStyle(25);
      h1_AJ[0][i]->SetMarkerColor(kGreen+2);
      h1_AJ[1][i]->SetMarkerStyle(20);
      h1_AJ[1][i]->SetMarkerColor(kMagenta+1);
      h1_AJ[0][i]->Draw("PE");
      h1_AJ[1][i]->Draw("SAME P E");
      sphenixtext(0.9,0.91,1);
      fitgdat = 1;
      fitgsim = 1;
      if(ajgaus[1][i])
	{
	  gmd = ajgaus[1][i]->GetParameter(1);
	  ged = ajgaus[1][i]->GetParameter(2);
	}
      if(ajgaus[0][i])
	{
	  gms = ajgaus[0][i]->GetParameter(1);
	  ges = ajgaus[0][i]->GetParameter(2);
	}
      for(int j=0; j<ntext; ++j)
	{
	  dotext[j] = 0;
	}
      dotext[0] = 1;
      if(i>3&&i<8)
	{
	  dotext[1] = 1;
	  dotext[2] = 1;
	}
      if(i>7)
	{
	  dotext[6] = 1;
	  dotext[7] = 1;
	}
      if(i%4>1)
	{
	  dotext[4] = 1;
	}
      if(i%4==1 || i%4==3)
	{
	  dotext[5] = 1;
	}
      //dotext[3] = 1;
      if(ajgaus[0][i] && ajgaus[1][i])
	{
	  ajgaus[0][i]->SetLineColor(kGreen+2);
	  ajgaus[0][i]->Draw("SAME");
	  ajgaus[1][i]->SetLineColor(kMagenta+1);
	  ajgaus[1][i]->Draw("SAME");
	}
      ajleg->Draw();
      std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
      d->SaveAs(("output/rmg/summed_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(i)+"_AJ_1d.png").c_str());
    }
  d->SetLogy();
  TLegend* jetEleg = new TLegend(0.55,0.45,0.9,0.7);
  jetTrigE[1]->SetMarkerStyle(24);
  jetTrigE[2]->SetMarkerStyle(21);
  jetTrigE[3]->SetMarkerStyle(20);
  jetTrigE[1]->SetLineWidth(2);
  jetTrigE[3]->SetLineWidth(2);
  tjetE->SetMarkerStyle(21);
  tjetE->SetMarkerColor(2);
  tjetE->SetMarkerSize(2);
  for(int i=1; i<4; ++i)
    {
      jetTrigE[i]->SetMarkerSize(2);
      jetTrigE[i]->Rebin(10);
    }
  jetTrigE[1]->SetMarkerColor(kBlue);
  jetTrigE[2]->SetMarkerColor(kBlue-7);
  jetTrigE[3]->SetMarkerColor(kBlue+3);
  jetEleg->AddEntry(jetE[0][0],"MBD N/S>=1 PYTHIA/Geant","p");
  jetEleg->AddEntry(jetE[1][0],"MBD N/S>=1 Data","p");
  jetEleg->AddEntry(jetTrigE[1],"MBD N/S>=1 \& Jet 8 GeV","p");
  jetEleg->AddEntry(jetTrigE[2],"MBD N/S>=1 \& Jet 10 GeV","p");
  jetEleg->AddEntry(jetTrigE[3],"MBD N/S>=1 \& Jet 12 GeV","p");
  jetEleg->AddEntry(tjetE,"MBD N/S >= 1 Truth Jets","p");
  jetEleg->SetFillStyle(0);
  jetEleg->SetBorderSize(0);
  jetEleg->SetFillColor(0);

  float jetx[9] = {8.43,10.38,12.79,15.68,19.22,23.76,29.25,36.03,44.16};
  float jety[9] = {0.085,0.0324,0.00793,0.00168,0.000474,0.00011,0.0000194,0.00000297,0.000000368};
  for(int i=0; i<9; ++i)
    {
      jety[i]/=960;
    }
  TGraph* pubjet = new TGraph(9,jetx,jety);
  pubjet->SetMarkerStyle(43);
  pubjet->SetMarkerColor(kCyan+2);
  pubjet->SetMarkerSize(3);
  jetEleg->AddEntry(pubjet,"Datathiefed STAR result","p");
  tjetE->Rebin(10);
  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<ntext; ++j)
	{
	  dotext[j] = 0;
	}
      if(i%4>1)
	{
	  dotext[4] = 1;
	}
      if(i%4==1 || i%4==3)
	{
	  dotext[5] = 1;
	}
      //dotext[3] = 1;
      dotext[0] = 1;
      dotext[8] = 1;
      float minimum = 1;
      jetE[0][i]->Rebin(10);
      jetE[1][i]->Rebin(10);

      for(int i=1; i<4; ++i)
	{
	  //cout << jetTrigE[i]->GetBinContent(100) << endl;
	  //cout << jetTrigE[i]->FindLastBinAbove() << endl;
	  float testval = jetTrigE[i]->GetBinContent(jetTrigE[i]->FindLastBinAbove());
	  //cout << jetTrigE[i]->GetBinContent(100) << endl;
	  minimum = min(minimum, testval);
	}
      //cout << minimum << endl;
      jetE[0][i]->GetYaxis()->SetRangeUser(0.5*minimum, 2*jetE[0][i]->GetMaximum());
      std_hist(jetE[1][i], jetE[0][i]);
      jetE[0][i]->Draw("PE");
      jetE[1][i]->Draw("SAME P E");
      if(i==0) tjetE->Draw("SAME P E");
      cout << "tjetE integral: " << tjetE->Integral() << endl;
      if(i==0)
	{
	  jetTrigE[1]->Draw("SAME P E");
	  jetTrigE[2]->Draw("SAME P E");
	  jetTrigE[3]->Draw("SAME P E");
	}
      pubjet->Draw("SAME P");
      jetEleg->Draw();
      fitgdat = 0;
      fitgsim = 0;
      std_text(d, texts, dotext, ntext, 0.03, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges,(therun==-1?-1:runnumbers[therun]));
      d->SaveAs(("output/rmg/summed_jetE"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(i)+".png").c_str());
    }

  for(int j=0; j<ntext; ++j)
    {
      dotext[j] = 0;
    }

  dotext[0] = 1;
  dotext[8] = 1;
  float minimum = 1;
  ljetE[0]->Rebin(10);
  ljetE[1]->Rebin(10);
  minimum = min(ljetE[0]->GetBinContent(ljetE[0]->FindLastBinAbove()),ljetE[1]->GetBinContent(ljetE[1]->FindLastBinAbove()));
  ljetE[0]->GetYaxis()->SetRangeUser(0.5*minimum, 2*max(ljetE[0]->GetMaximum(),ljetE[1]->GetMaximum()));
  std_hist(ljetE[1], ljetE[0]);
  
  ljetE[0]->Draw("PE");
  ljetE[1]->Draw("SAME P E");
  ajleg->Draw();
  fitgdat = 0;
  fitgsim = 0;
  std_text(d, texts, dotext, ntext, 0.03, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
  d->SaveAs(("output/rmg/summed"+to_string((therun==-1?-1:runnumbers[therun]))+"_ljetE.png").c_str());


  for(int j=0; j<ntext; ++j)
    {
      dotext[j] = 0;
    }
  dotext[0] = 1;
  dotext[8] = 1;
  minimum = 1;
  minimum = min(ljetEta[0]->GetBinContent(ljetEta[0]->FindLastBinAbove()),ljetEta[1]->GetBinContent(ljetEta[1]->FindLastBinAbove()));
  ljetEta[0]->GetYaxis()->SetRangeUser(0, 1.1*max(ljetEta[0]->GetMaximum(),ljetEta[1]->GetMaximum()));
  std_hist(ljetEta[1], ljetEta[0]);
  gPad->SetLogy(0);
  ljetEta[0]->Draw("PE");
  ljetEta[1]->Draw("SAME P E");
  ajleg->Draw();
  fitgdat = 0;
  fitgsim = 0;
  std_text(d, texts, dotext, ntext, 0.03, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
  d->SaveAs(("output/rmg/summed"+to_string((therun==-1?-1:runnumbers[therun]))+"_ljetEta.png").c_str());
  gPad->SetLogy();
  //jetEleg->AddEntry(jetTrigE[0],"Jet 4 GeV Trigger","p");
  //jetEleg->AddEntry(jetTrigE[1],"Jet 6 GeV Trigger","p");
  //jetEleg->AddEntry(jetTrigE[2],"Jet 8 GeV Trigger","p");
  //jetEleg->AddEntry(jetTrigE[3],"Jet 10 GeV Trigger","p");
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  //	  h1_rej[h][i]->Scale(1./tree->GetEntries());
	}
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
  TLegend* leg = new TLegend(0.1,0.57,0.5,0.87);
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
  leg->Draw();  
  d->SaveAs(("output/rmg/summed"+to_string((therun==-1?-1:runnumbers[therun]))+"_rej.png").c_str());


  for(int i=0; i<4; ++i)
    {
      for(int j=0; j<ntext; ++j)
	{
	  dotext[j] = 0;
	}
      if(i%4>1)
	{
	  dotext[4] = 1;
	}
      if(i%4==1 || i%4==3)
	{
	  dotext[5] = 1;
	}
      dotext[0] = 1;
      //dotext[8] = 1;
      //dotext[3] = 1;
      fitgdat = 0;
      fitgsim = 0;
      h1_mlt[0][i]->GetXaxis()->SetTitle("#it{E}_{T,jet} > 5 GeV Multiplicity");
      std_hist(h1_mlt[1][i], h1_mlt[0][i]);
      h1_mlt[0][i]->Draw("PE");
      h1_mlt[1][i]->Draw("SAME P E");
      ajleg->Draw();
      std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
      d->SaveAs(("output/rmg/summed_mult"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(i)+".png").c_str());
    }


  gPad->SetLogy(0);
  for(int j=0; j<ntext; ++j)
    {
      dotext[j] = 0;
    }
  fitgdat = 0;
  fitgsim = 0;
  h1_zdist[0]->GetYaxis()->SetRangeUser(0,1.1*max(h1_zdist[0]->GetMaximum(),h1_zdist[1]->GetBinContent(101)));
  h1_zdist[0]->GetXaxis()->SetTitle("Z-vertex Position (cm)");
  h1_zdist[0]->GetYaxis()->SetTitle("Counts");
  std_hist(h1_zdist[1], h1_zdist[0]);
  h1_zdist[0]->Draw("PE");
  h1_zdist[1]->Draw("SAME P E");
  if(zgaus[0] && zgaus[1])
    {
      gmd = zgaus[1]->GetParameter(1);
      ged = zgaus[1]->GetParameter(2);
      gms = zgaus[0]->GetParameter(1);
      ges = zgaus[0]->GetParameter(2);
    }
  ajleg->Draw();
  fitgsim = 1;
  fitgdat = 1;
  std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
  d->SaveAs(("output/rmg/summed"+to_string((therun==-1?-1:runnumbers[therun]))+"_zdist.png").c_str());
  fitgsim = 0;
  fitgdat = 0;
  //cout << h1_zdist[0]->GetBinContent(0) << " " << h1_zdist[1]->GetBinContent(0) << " " << h1_zdist[0]->GetBinContent(201) << " " << h1_zdist[1]->GetBinContent(201) << endl;
  gPad->SetLogy();
  for(int j=0; j<ntext; ++j)
    {
      dotext[j] = 0;
    }
  fitgdat = 0;
  fitgsim = 0;

  //  h1_mbdq[0]->Scale(1./h1_mbdq[0]->Integral(3,500));
  //  h1_mbdq[1]->Scale(1./h1_mbdq[1]->Integral(3,500));
  h1_mbdq[0]->GetYaxis()->SetRangeUser(0.1,1.1*max(h1_mbdq[1]->GetMaximum(),h1_mbdq[0]->GetMaximum()));
  //  h1_mbdq[0]->GetYaxis()->SetRangeUser(0,0.1);  
  //h1_mbdq[0]->GetXaxis()->SetRangeUser(0,1);
  h1_mbdq[0]->GetXaxis()->SetTitle("MBD PMT Charge [arb]");
  h1_mbdq[0]->GetYaxis()->SetTitle("Counts");
  std_hist(h1_mbdq[1], h1_mbdq[0]);
  h1_mbdq[0]->Draw("PE");
  //h1_mbdq[1]->Draw("SAME P E");
  ajleg->Draw();
  std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
  d->SaveAs(("output/rmg/summed"+to_string((therun==-1?-1:runnumbers[therun]))+"_mbdq.png").c_str());
  
  gPad->SetLogy();
  
  for(int h=0; h<5; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  for(int j=0; j<ntext; ++j)
	    {
	      dotext[j] = 0;
	    }
	  if(i%4>1)
	    {
	      dotext[4] = 1;
	    }
	  if(i%4==1 || i%4==3)
	    {
	      dotext[5] = 1;
	    }
	  dotext[0] = 1;
	  dotext[9+h] = 1;
	  //dotext[3] = 1;
	  fitgdat = 0;
	  fitgsim = 0;
	  std_hist(h1_eta[1][i][h], h1_eta[0][i][h]);
	  gPad->SetLogy(0);
	  h1_eta[0][i][h]->GetYaxis()->SetRangeUser(0,1.1*max(h1_eta[0][i][h]->GetMaximum(),h1_eta[1][i][h]->GetMaximum()));
	  h1_eta[0][i][h]->GetXaxis()->SetTitle("#eta");
	  h1_eta[0][i][h]->GetYaxis()->SetTitle("Integral Normalized Counts");
	  h1_eta[0][i][h]->Draw("PE");
	  h1_eta[1][i][h]->Draw("SAME P E");
	  ajleg->Draw();
	  std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
	  d->SaveAs(("output/rmg/summed_eta"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(i)+"_"+to_string(h)+".png").c_str());
	}
      
      h2_cal_eta_phi[1][0]->GetZaxis()->SetRangeUser(0,h2_cal_eta_phi[1][0]->GetMaximum()/2);
      for(int i=0; i<4; ++i)
	{
	  //cout << "grombo2 " << i << endl;
	  for(int j=0; j<ntext; ++j)
	    {
	      dotext[j] = 0;
	    }
	  fitgdat = 0;
	  fitgsim = 0;
	  gPad->SetRightMargin(0.2);
	  if(i<3)
	    {
	      h2_cal_eta_phi[0][i]->GetZaxis()->SetRangeUser(0,h2_cal_eta_phi[0][i]->GetMaximum());
	      h2_cal_eta_phi[0][i]->GetXaxis()->SetTitle("#eta");
	      h2_cal_eta_phi[0][i]->GetYaxis()->SetTitle("#phi [rad]");
	      h2_cal_eta_phi[0][i]->Draw("COLZ");
	      d->SaveAs(("output/rmg/cal_eta_phi_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(0)+"_"+to_string(i)+".png").c_str());
	      h2_cal_eta_phi[1][i]->GetXaxis()->SetTitle("#eta");
	      h2_cal_eta_phi[1][i]->GetYaxis()->SetTitle("#phi [rad]");
	      h2_cal_eta_phi[1][i]->GetZaxis()->SetRangeUser(0,h2_cal_eta_phi[1][i]->GetMaximum());
	      h2_cal_eta_phi[1][i]->Draw("COLZ");
	      d->SaveAs(("output/rmg/cal_eta_phi_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(1)+"_"+to_string(i)+".png").c_str());
	    }
	  //cout << "grombo0" << endl;
	  gPad->SetRightMargin(0.2);
	  h2_tjet_eta_phi->GetZaxis()->SetRangeUser(0,h2_tjet_eta_phi->GetMaximum());
	  h2_tjet_eta_phi->GetZaxis()->SetTitleOffset(1.7);
	  h2_tjet_eta_phi->GetXaxis()->SetTitle("#eta");
	  h2_tjet_eta_phi->GetYaxis()->SetTitle("#phi [rad]");
	  h2_tjet_eta_phi->Draw("COLZ");
	  h2_tjet_eta_phi->GetZaxis()->SetTitle("N_{jet} Event Normalized");
	  d->SaveAs(("output/rmg/tjet_eta_phi_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(0)+"_"+to_string(i)+".png").c_str());

	  h2_jet_eta_phi[0][i]->GetZaxis()->SetRangeUser(0,h2_jet_eta_phi[0][i]->GetMaximum());
	  h2_jet_eta_phi[0][i]->GetZaxis()->SetTitleOffset(1.7);
	  h2_jet_eta_phi[0][i]->GetXaxis()->SetTitle("#eta");
	  h2_jet_eta_phi[0][i]->GetYaxis()->SetTitle("#phi [rad]");
	  h2_jet_eta_phi[0][i]->Draw("COLZ");
	  h2_jet_eta_phi[0][i]->GetZaxis()->SetTitle("N_{jet} Event Normalized");
	  d->SaveAs(("output/rmg/jet_eta_phi_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(0)+"_"+to_string(i)+".png").c_str());
	  h2_jet_eta_phi[0][i]->Draw("LEGO");
	  d->SaveAs(("output/rmg/jet_eta_phi_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(0)+"_"+to_string(i)+"_lego.png").c_str());
	  h2_jet_eta_phi[1][i]->GetZaxis()->SetTitleOffset(1.7);
	  h2_jet_eta_phi[1][i]->GetXaxis()->SetTitle("#eta");
	  h2_jet_eta_phi[1][i]->GetYaxis()->SetTitle("#phi [rad]");
	  h2_jet_eta_phi[1][i]->GetZaxis()->SetTitle("N_{jet} Event Normalized");
	  h2_jet_eta_phi[1][i]->GetZaxis()->SetRangeUser(0,h2_jet_eta_phi[1][i]->GetMaximum());
	  h2_jet_eta_phi[1][i]->Draw("COLZ");
	  d->SaveAs(("output/rmg/jet_eta_phi_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(1)+"_"+to_string(i)+".png").c_str());
	  h2_jet_eta_phi[1][i]->Draw("LEGO");
	  d->SaveAs(("output/rmg/jet_eta_phi_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(1)+"_"+to_string(i)+"_lego.png").c_str());
	  h2_jet_eta_e[0][i]->GetZaxis()->SetRangeUser(0,h2_jet_eta_e[0][i]->GetMaximum());
	  h2_jet_eta_e[0][i]->GetXaxis()->SetTitle("#eta");
	  h2_jet_eta_e[0][i]->GetYaxis()->SetTitle("#it{E}_{jet} [GeV]");
	  h2_jet_eta_e[0][i]->Draw("COLZ");
	  d->SaveAs(("output/rmg/jet_eta_e_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(0)+"_"+to_string(i)+".png").c_str());
	  //cout << "grombo3" << endl;
	  h2_jet_eta_e[1][i]->GetXaxis()->SetTitle("#eta");
	  h2_jet_eta_e[1][i]->GetYaxis()->SetTitle("#it{E}_{jet} [GeV]");
	  h2_jet_eta_e[1][i]->GetZaxis()->SetRangeUser(0,h2_jet_eta_e[1][i]->GetMaximum());
	  h2_jet_eta_e[1][i]->Draw("COLZ");
	  d->SaveAs(("output/rmg/jet_eta_e_"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(1)+"_"+to_string(i)+".png").c_str());
	  gPad->SetRightMargin(0.1);
	  //cout << "grombo1" << endl;
	  if(i>3) continue;
	  if(i%4>1)
	    {
	      dotext[4] = 1;
	    }
	  if(i%4==1 || i%4==3)
	    {
	      dotext[5] = 1;
	    }
	  dotext[0] = 1;
	  dotext[9+h] = 1;
	  //dotext[3] = 1;
	  fitgdat = 0;
	  fitgsim = 0;
	  std_hist(h1_phi[1][i][h], h1_phi[0][i][h]);
	  h1_phi[0][i][h]->GetYaxis()->SetTitle("Integral Normalized Counts");
	  h1_phi[0][i][h]->GetYaxis()->SetRangeUser(0,1.1*max(h1_phi[0][i][h]->GetMaximum(),h1_phi[1][i][h]->GetMaximum()));
	  h1_phi[0][i][h]->GetXaxis()->SetTitle("EM Fraction");
	  h1_phi[0][i][h]->Draw("PE");
	  h1_phi[1][i][h]->Draw("SAME P E");
	  ajleg->Draw();
	  std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
	  d->SaveAs(("output/rmg/summed_phi"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(i)+"_"+to_string(h)+".png").c_str());
	}
    }
  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<ntext; ++j)
	{
	  dotext[j] = 0;
	}
      dotext[0] = 1;
      dotext[9+i] = 1;
      //dotext[3] = 1;
      fitgdat = 0;
      fitgsim = 0;
      std_hist(jetfrac[1][i], jetfrac[0][i]);
      jetfrac[0][i]->Rebin(4);
      jetfrac[1][i]->Rebin(4);
      jetfrac[0][i]->Scale(1./jetfrac[0][i]->Integral());
      jetfrac[1][i]->Scale(1./jetfrac[1][i]->Integral());
      jetfrac[0][i]->SetLineColor(kGreen+2);
      jetfrac[0][i]->SetFillColorAlpha(kGreen+2,0.3);
      jetfrac[1][i]->SetLineColor(kMagenta+1);
      jetfrac[1][i]->SetFillColorAlpha(kMagenta+1,0.3);
      jetfrac[0][i]->GetYaxis()->SetRangeUser(0,1.1*max(jetfrac[0][i]->GetMaximum(),jetfrac[1][i]->GetMaximum()));
      jetfrac[0][i]->GetXaxis()->SetTitle("EM Fraction");
      jetfrac[0][i]->Draw("HIST");
      jetfrac[1][i]->Draw("SAME HIST");
      ajleg->Draw();
      std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
      d->SaveAs(("output/rmg/summed_emfrac"+to_string((therun==-1?-1:runnumbers[therun]))+"_"+to_string(i)+".png").c_str());
    }
  cout << h1_zdist[0]->Integral() << endl;
  cout << h1_zdist[1]->Integral(3,200) << endl;
  cout << h1_zdist[1]->GetEntries() << endl;
  return 0;
}
