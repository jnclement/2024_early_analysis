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
float min(float a, float b)
{
  if(a<b) return a;
  return b;
}
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

void std_text(TCanvas* thecan, string* texts, int* dotext, int ntext, float textsize, float textx, float texty, int rightalign, int nevtsim, int nevtdat, int fitgdat, float gmd, float ged, int fitgsim, float gms, float ges)
{
  std::stringstream stream[4];
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
  thecan->cd();
  sphenixtext(0.9,0.91,1);
  drawText(("n_{evt,Pythia} = "+to_string(nevtsim)+" n_{evt,data} = "+to_string(nevtdat)).c_str(),0.1,0.96,0,kBlack,0.03);
  drawText("Good runs 45200-45300",0.2,0.91,0,kBlack,0.03);
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
      drawText(("#mu_{data} = "+stream[0].str()+" #pm "+stream[1].str()).c_str(),textx,drawy,rightalign,kBlack,textsize);
      drawy -= 5*textsize/4;
    }
  if(fitgsim)
    {
      drawText(("#mu_{sim} = "+stream[2].str()+" #pm "+stream[3].str()).c_str(),textx,drawy,rightalign,kBlack,textsize);
    }
}

void std_hist(TH1D* dathist, TH1D* simhist)
{
  simhist->GetYaxis()->SetTitle("Event Normalized Counts");
  simhist->SetMarkerSize(3);
  dathist->SetMarkerSize(3);
  simhist->SetMarkerStyle(39);
  simhist->SetMarkerColor(kGreen+2);
  dathist->SetMarkerStyle(43);
  dathist->SetMarkerColor(kMagenta+3);
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

int plot(string filebase = "summed_dat.root", string filelist="sumdatlist.list",string sfilebase="summed_sim.root")
{
  const int nruns = 27;
  const int nsd = 4;
  int sds[nruns][nsd] = {0};
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
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0);
  TFile* rootfile[3];
  rootfile[0] = TFile::Open(sfilebase.c_str());
  rootfile[1] = TFile::Open(filebase.c_str());
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
  while(getline(file,line))
    {
      rootfile[2] = TFile::Open(line.c_str());
      TTree* tree = (TTree*)rootfile[2]->Get("outt");
      tree->SetBranchAddress("outevt",&evtct);
      tree->SetBranchAddress("outnmb",&mbevt);
      tree->SetBranchAddress("njmb",&njmb);
      tree->SetBranchAddress("nJetTrig",nJetTrig);

      for(int i=0; i<tree->GetEntries(); ++i)
	{
	  tree->GetEntry(i);
	  nevents += evtct;
	  nmb += mbevt;
	  runmb += mbevt;
	  jmb += njmb;
	}
      whichrun++;
      effevt[0] += (1.*(1+sds[whichrun][0]))/(1.*(1+sds[whichrun][1]))*runmb;//(nJetTrig[1]+mbevt);
      effevt[1] += (1.*(1+sds[whichrun][0]))/(1.*(1+sds[whichrun][2]))*runmb;//(nJetTrig[2]+mbevt);
      effevt[2] += (1.*(1+sds[whichrun][0]))/(1.*(1+sds[whichrun][3]))*runmb;//(nJetTrig[3]+mbevt);
      runmb = 0;
    }

  cout << whichrun << endl;
  //effevt[0];
  //effevt[1];
  //effevt[2];
  cout << nmb  << " " << nevents << " " << 10000000 << " " << effevt[0] << " " << nJetTrig[1] << endl;
  
  TH1D* h1_dphi[2][4];
  TH1D* h1_rej[2][4];
  TF1* ajgaus[2][12];
  TF1* dphigaus[2][4];
  TH1D* h1_mlt[2][4];
  TH1D* h1_phi[2][4];
  TH1D* h1_eta[2][4];
  TH1D* jetfrac[2][4];
  TH1D* jetTrigE[4];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  if(h==1 && i > 0)
	    {
	      jetTrigE[i] = (TH1D*)rootfile[h]->Get(("jetTrigE"+to_string(i)).c_str());
	      jetTrigE[i]->Scale(1./effevt[i-1]);
	    }
	  
	  h1_dphi[h][i] = (TH1D*)rootfile[h]->Get(("h1_dphi"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_dphi[h][i]->Scale(1./(h==0?10000000:nmb));
	  h1_dphi[h][i]->Fit("gaus","","",1.,5.);
	  dphigaus[h][i] = h1_dphi[h][i]->GetFunction("gaus");
	  h1_rej[h][i] = (TH1D*)rootfile[h]->Get(("h1_rej"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_rej[h][i]->Rebin(5);
	  h1_mlt[h][i] = (TH1D*)rootfile[h]->Get(("h1_mlt"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_phi[h][i] = (TH1D*)rootfile[h]->Get(("h1_phi"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_eta[h][i] = (TH1D*)rootfile[h]->Get(("h1_eta"+to_string(1)+"_"+to_string(i)).c_str());
	  jetfrac[h][i] = (TH1D*)rootfile[h]->Get(("jetfrac"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_mlt[h][i]->Scale(1./(h==0?(10000000):nmb));
	  h1_phi[h][i]->Scale(1./(h==0?(10000000):nmb));
	  h1_eta[h][i]->Scale(1./(h==0?(10000000):nmb));
	  if(i<3) jetfrac[h][i]->Scale(1./(h==0?(10000000):nmb));
	}
    }
  TH1D* h1_AJ[2][12];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  h1_AJ[h][i] = (TH1D*)rootfile[h]->Get(("h1_AJ"+to_string(1)+"_"+to_string(i)).c_str());
	  h1_AJ[h][i]->Scale(1./h1_AJ[h][i]->Integral());//1./(h==0?(tree->GetEntries()*10000000):nmb));
	  h1_AJ[h][i+4] = (TH1D*)rootfile[h]->Get(("h1_AJ"+to_string(1)+"_"+to_string(i+4)).c_str());
	  cout << h1_AJ[h][i+4] << endl;
	  h1_AJ[h][i+4]->Scale(1./h1_AJ[h][i+4]->Integral());//1./(h==0?(tree->GetEntries()*10000000):nmb));
	  h1_AJ[h][i]->Fit("gaus");
	  ajgaus[h][i] = h1_AJ[h][i]->GetFunction("gaus");
	  h1_AJ[h][i+4]->Fit("gaus");
	  ajgaus[h][i+4] = h1_AJ[h][i+4]->GetFunction("gaus");
	  h1_AJ[h][i+8] = (TH1D*)rootfile[h]->Get(("h1_AJ"+to_string(1)+"_"+to_string(i+8)).c_str());
	  h1_AJ[h][i+8]->Scale(1./h1_AJ[h][i+8]->Integral());//1./(h==0?(tree->GetEntries()*10000000):nmb));
	  h1_AJ[h][i+8]->Fit("gaus");
	  ajgaus[h][i+8] = h1_AJ[h][i+8]->GetFunction("gaus");
	}
    }
  TH1D* jetE[2][4];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  jetE[h][i] = (TH1D*)rootfile[h]->Get(("jetE"+to_string(1)+"_"+to_string(i)).c_str());
	  jetE[h][i]->Scale(1./(h==0?(10000000):nmb));
	}
    }
  TCanvas* d = new TCanvas("","",1000,1000);
  //  TH1D* jetTrigE[4];
  /*
  for(int i=0; i<4; i++)
    {
      jetTrigE[i] = (TH1D*)rootfile[h]->Get(("jetTrigE"+to_string(i)).c_str());
    }
  */
  const int ntext = 12;
  string texts[ntext] =
    {
      "Anti-#it{k}_{T} #it{R}=0.4",
      "#it{E}_{jet,1} > 10 GeV",
      "#it{E}_{jet,2} > 5 GeV",
      "|#it{#eta}^{jet}| < 0.7",
      "E_{lead tower} < 0.65 E_{jet}",
      "E_{EM} < "+to_string(maxeh)+" E_{HAD}",
      "Leading Jet E > 7 GeV",
      "Subleading Jet E > 4 GeV",
      "E_{jet} > 4 GeV",
      "E_{jet} < 7 GeV",
      "7 GeV < E_{jet} < 10 GeV",
      "E_{jet} > 10 GeV"
    };

  int dotext[ntext] = {1,1,1,1,1,1,1,1};
  float stdsize = 0.02;
  float stdx = 0.87;
  float stdy = 0.88-stdsize;
  int stdright = 1;
  int nevtsim = 10000000;
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
	  dphigaus[h][i]->SetLineColor(kMagenta);
	  if(i>1) dotext[4] = 1;
	  else dotext[4] = 0;
	  if(i==1 || i==3) dotext[5] = 1;
	  else dotext[5] = 0;
	  gmd = dphigaus[h][i]->GetParameter(1);
	  ged = dphigaus[h][i]->GetParError(1);
	  std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
	  d->SaveAs(("output/rmg/summed_"+to_string(h)+"_"+to_string(i)+"_dphi_1d.pdf").c_str());
	}
    }
  TLegend* ajleg = new TLegend(0,0,0.2,0.09);
  ajleg->SetFillStyle(0);
  ajleg->SetFillColor(0);
  ajleg->SetTextFont(42);
  ajleg->SetBorderSize(0);
  ajleg->AddEntry(h1_AJ[0][0],"Minbias Pythia","p");
  ajleg->AddEntry(h1_AJ[1][0],"Minbias Data","p");
  for(int i=0; i<12; ++i)
    {
      h1_AJ[0][i]->GetYaxis()->SetTitle("Event Normalized Counts");
      h1_AJ[0][i]->SetMarkerSize(3);
      h1_AJ[1][i]->SetMarkerSize(3);
      h1_AJ[0][i]->SetMarkerStyle(39);
      h1_AJ[0][i]->SetMarkerColor(kGreen+2);
      h1_AJ[1][i]->SetMarkerStyle(43);
      h1_AJ[1][i]->SetMarkerColor(kMagenta+3);
      h1_AJ[0][i]->Draw("PE");
      h1_AJ[1][i]->Draw("SAME P E");
      sphenixtext(0.9,0.91,1);
      fitgdat = 1;
      fitgsim = 1;
      gmd = ajgaus[1][i]->GetParameter(1);
      ged = ajgaus[1][i]->GetParError(1);
      gms = ajgaus[0][i]->GetParameter(1);
      ges = ajgaus[0][i]->GetParError(1);
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
      ajgaus[0][i]->SetLineColor(kGreen+2);
      ajgaus[0][i]->Draw("SAME");
      ajgaus[1][i]->SetLineColor(kMagenta+3);
      ajgaus[1][i]->Draw("SAME");
      ajleg->Draw();
      std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
      d->SaveAs(("output/rmg/summed_"+to_string(i)+"_AJ_1d.pdf").c_str());
    }
  d->SetLogy();
  TLegend* jetEleg = new TLegend(0.65,0.5,0.85,0.65);
  for(int i=1; i<4; ++i)
    {
      jetTrigE[i]->SetMarkerStyle(38+i);
      jetTrigE[i]->SetMarkerColor(kRed+i);
      jetTrigE[i]->SetMarkerSize(2);
      jetTrigE[i]->Rebin(10);
    }
  jetEleg->AddEntry(jetE[0][0],"Minbias Pythia","p");
  jetEleg->AddEntry(jetE[1][0],"Minbias Data","p");
  jetEleg->AddEntry(jetTrigE[1],"MBD N/S>1 \& Jet 6 GeV","p");
  jetEleg->AddEntry(jetTrigE[2],"MBD N/S>1 \& Jet 8 GeV","p");
  jetEleg->AddEntry(jetTrigE[3],"MBD N/S>1 \& Jet 10 GeV","p");
  jetEleg->SetFillStyle(0);
  jetEleg->SetBorderSize(0);
  jetEleg->SetFillColor(0);
  

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
      dotext[0] = 1;
      dotext[8] = 1;
      float minimum = 1;
      jetE[0][i]->Rebin(10);
      jetE[1][i]->Rebin(10);
      for(int i=1; i<4; ++i)
	{
	  cout << jetTrigE[i]->GetBinContent(100) << endl;
	  cout << jetTrigE[i]->FindLastBinAbove() << endl;
	  float testval = jetTrigE[i]->GetBinContent(jetTrigE[i]->FindLastBinAbove());
	  //cout << jetTrigE[i]->GetBinContent(100) << endl;
	  minimum = min(minimum, testval);
	}
      cout << minimum << endl;
      jetE[0][i]->GetYaxis()->SetRangeUser(0.5*minimum, 2*jetE[1][i]->GetMaximum());
      std_hist(jetE[1][i], jetE[0][i]);
      jetE[0][i]->Draw("PE");
      jetE[1][i]->Draw("SAME P E");
      if(i==0)
	{
	  jetTrigE[1]->Draw("SAME P E");
	  jetTrigE[2]->Draw("SAME P E");
	  jetTrigE[3]->Draw("SAME P E");
	}
      jetEleg->Draw();
      fitgdat = 0;
      fitgsim = 0;
      std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
      d->SaveAs(("output/rmg/summed_jetE"+to_string(i)+".pdf").c_str());
    }
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
  d->SaveAs("output/rmg/summed_rej.pdf");


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
      dotext[8] = 1;
      fitgdat = 0;
      fitgsim = 0;
      h1_mlt[0][i]->GetXaxis()->SetTitle("Jet > 5 GeV Multiplicity");
      std_hist(h1_mlt[1][i], h1_mlt[0][i]);
      h1_mlt[0][i]->Draw("PE");
      h1_mlt[1][i]->Draw("SAME P E");
      ajleg->Draw();
      std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
      d->SaveAs(("output/rmg/summed_mult"+to_string(i)+".pdf").c_str());
    }

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
      dotext[8] = 1;
      fitgdat = 0;
      fitgsim = 0;
      std_hist(h1_eta[1][i], h1_eta[0][i]);
      gPad->SetLogy(0);
      h1_eta[0][i]->GetYaxis()->SetRangeUser(0,1.1*max(h1_eta[0][i]->GetMaximum(),h1_eta[1][i]->GetMaximum()));
      h1_eta[0][i]->GetXaxis()->SetTitle("#eta");
      h1_eta[0][i]->Draw("PE");
      h1_eta[1][i]->Draw("SAME P E");
      ajleg->Draw();
      std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
      d->SaveAs(("output/rmg/summed_eta"+to_string(i)+".pdf").c_str());
    }

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
      dotext[8] = 1;
      fitgdat = 0;
      fitgsim = 0;
      std_hist(h1_phi[1][i], h1_phi[0][i]);
      h1_phi[0][i]->GetYaxis()->SetRangeUser(0,1.1*max(h1_phi[0][i]->GetMaximum(),h1_phi[1][i]->GetMaximum()));
      h1_phi[0][i]->GetXaxis()->SetTitle("EM Fraction");
      h1_phi[0][i]->Draw("PE");
      h1_phi[1][i]->Draw("SAME P E");
      ajleg->Draw();
      std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
      d->SaveAs(("output/rmg/summed_phi"+to_string(i)+".pdf").c_str());
    }

  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<ntext; ++j)
	{
	  dotext[j] = 0;
	}
      dotext[0] = 1;
      dotext[9+i] = 1;
      fitgdat = 0;
      fitgsim = 0;
      std_hist(jetfrac[1][i], jetfrac[0][i]);
      jetfrac[0][i]->Rebin(4);
      jetfrac[1][i]->Rebin(4);
      jetfrac[0][i]->Scale(1./jetfrac[0][i]->Integral());
      jetfrac[1][i]->Scale(1./jetfrac[1][i]->Integral());
      jetfrac[0][i]->SetLineColor(kGreen+2);
      jetfrac[0][i]->SetFillColorAlpha(kGreen+2,0.3);
      jetfrac[1][i]->SetLineColor(kMagenta+3);
      jetfrac[1][i]->SetFillColorAlpha(kMagenta+3,0.3);
      jetfrac[0][i]->GetYaxis()->SetRangeUser(0,1.1*max(jetfrac[0][i]->GetMaximum(),jetfrac[1][i]->GetMaximum()));
      jetfrac[0][i]->GetXaxis()->SetTitle("EM Fraction");
      jetfrac[0][i]->Draw("HIST");
      jetfrac[1][i]->Draw("SAME HIST");
      ajleg->Draw();
      std_text(d, texts, dotext, ntext, stdsize, stdx, stdy, stdright, nevtsim, nevtdat, fitgdat, gmd, ged, fitgsim, gms, ges);
      d->SaveAs(("output/rmg/summed_emfrac"+to_string(i)+".pdf").c_str());
    }
  return 0;
}
