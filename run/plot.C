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

int plot(string filebase="summed.root")
{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0);
  TFile* rootfile = TFile::Open(filebase.c_str());
  TTree* tree = (TTree*)rootfile->Get("outt");;
  long long unsigned int nevents = 0;
  int evtct = 0;
  int mbevt = 0;
  int njmb = 0;
  int jmb = 0;
  int maxeh = 5;
  long long unsigned int nmb = 0;
  tree->SetBranchAddress("outevt",&evtct);
  tree->SetBranchAddress("outnmb",&mbevt);
  tree->SetBranchAddress("njmb",&njmb);
  int nhists = tree->GetEntries();
  for(int i=0; i<tree->GetEntries(); ++i)
    {
      tree->GetEntry(i);
      nevents += evtct;
      nmb += mbevt;
      jmb += njmb;
    }
  cout << nmb  << " " << nevents << " " << 10000000 << endl;
  TH1D* h1_dphi[2][4];
  TH1D* h1_rej[2][4];
  TF1* ajgaus[2][8];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  h1_dphi[h][i] = (TH1D*)rootfile->Get(("h1_dphi"+to_string(h)+"_"+to_string(i)).c_str());
	  h1_dphi[h][i]->Scale(1./(h==0?(tree->GetEntries()*10000000):nmb));
	  h1_rej[h][i] = (TH1D*)rootfile->Get(("h1_rej"+to_string(h)+"_"+to_string(i)).c_str());
	  h1_rej[h][i]->Rebin(5);
	}
    }
  TH1D* h1_AJ[2][8];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  h1_AJ[h][i] = (TH1D*)rootfile->Get(("h1_AJ"+to_string(h)+"_"+to_string(i)).c_str());
	  h1_AJ[h][i]->Scale(1./(h==0?(tree->GetEntries()*10000000):nmb));
	  h1_AJ[h][i+4] = (TH1D*)rootfile->Get(("h1_AJ"+to_string(h)+"_"+to_string(i+4)).c_str());
	  h1_AJ[h][i+4]->Scale(1./(h==0?(tree->GetEntries()*10000000):nmb));
	  h1_AJ[h][i]->Fit("gaus");
	  ajgaus[h][i] = h1_AJ[h][i]->GetFunction("gaus");
	  h1_AJ[h][i+4]->Fit("gaus");
	  ajgaus[h][i+4] = h1_AJ[h][i+4]->GetFunction("gaus");
	  h1_AJ[h][i+8] = (TH1D*)rootfile->Get(("h1_AJ"+to_string(h)+"_"+to_string(i+8)).c_str());
	  h1_AJ[h][i+8]->Scale(1./(h==0?(tree->GetEntries()*10000000):nmb));
	  h1_AJ[h][i+8]->Fit("gaus");
	  ajgaus[h][i+8] = h1_AJ[h][i+8]->GetFunction("gaus");
	}
    }
  TH1D* jetE[2][4];
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<3; ++i)
	{
	  jetE[h][i] = (TH1D*)rootfile->Get(("jetE"+to_string(h)+"_"+to_string(i)).c_str());
	  jetE[h][i]->Scale(1./(h==0?(tree->GetEntries()*10000000):nmb));
	}
    }
  TCanvas* d = new TCanvas("","",1000,1000);
  TH1D* jetTrigE[4];
  for(int i=0; i<4; i++)
    {
      jetTrigE[i] = (TH1D*)rootfile->Get(("jetTrigE"+to_string(i)).c_str());
    }
  for(int h=0; h<2; ++h)
    {
      gPad->SetLogy(0);
      for(int i=0; i<4; ++i)
	{
	  h1_dphi[h][i]->SetMarkerSize(3);
	  h1_dphi[h][i]->Draw("PE");
	  sphenixtext(0.96,0.96,1);
	  drawText("#it{E}_{jet,1} > 10 GeV", 0.93, 0.9, 1, kBlack, 0.04);
	  drawText("#it{E}_{jet,2} > 5 GeV",0.93,0.85,1,kBlack,0.04);
	  drawText("Anti-#it{k}_{T} #it{R}=0.4",0.93,0.8,1,kBlack,0.04);
	  drawText("|#it{#eta}^{jet}| < 0.7",0.93,0.75,1,kBlack,0.04);
	  if(i>1) drawText("E_{lead tower} < 0.65 E_{jet}", 0.93, 0.7, 1, kBlack, 0.04);
	  if(i==1 || i==3) drawText(("E_{EM} < "+to_string(maxeh)+" E_{HAD}").c_str(), 0.93, (i==1?0.66:0.6), 1, kBlack, 0.04);
	  d->SaveAs(("output/rmg/summed_"+to_string(h)+"_"+to_string(i)+"_dphi_1d.pdf").c_str());
	}
    }
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
      sphenixtext(0.96,0.96,1);
      drawText(("#mu_{data} = "+to_string(ajgaus[1][i]->GetParameter(1))+"#pm"+to_string(ajgaus[1][i]->GetParError(1))).c_str(),0.87,0.65,1,kBlack,0.04);
      drawText(("#mu_{sim} = "+to_string(ajgaus[0][i]->GetParameter(1))+"#pm"+to_string(ajgaus[0][i]->GetParError(1))).c_str(),0.87,0.6,1,kBlack,0.04);
      if(i>3&&i<8)drawText("Leading Jet E > 10 GeV", 0.87, 0.85, 1, kBlack, 0.04);
      if(i>3&&i<8)drawText("Subleading Jet E > 5 GeV", 0.87, 0.8, 1, kBlack, 0.04);
      if(i>7)drawText("Leading Jet E > 7 GeV", 0.87, 0.85, 1, kBlack, 0.04);
      if(i>7)drawText("Subleading Jet E > 4 GeV", 0.87, 0.8, 1, kBlack, 0.04);
      if(i%4>1) drawText("E_{lead tower} < 0.65 E_{jet}", 0.87, 0.75, 1, kBlack, 0.04);
      if(i%4==1 || i%4==3) drawText(("E_{EM} < "+to_string(maxeh)+" E_{HAD}").c_str(), 0.87, (i==1?0.75:0.7), 1, kBlack, 0.04);
      ajgaus[0][i]->SetLineColor(kGreen+2);
      ajgaus[0][i]->Draw("SAME");
      ajgaus[1][i]->SetLineColor(kMagenta+3);
      ajgaus[1][i]->Draw("SAME");
      d->SaveAs(("output/rmg/summed_"+to_string(i)+"_AJ_1d.pdf").c_str());
    }
  d->SetLogy();
  jetE[0][0]->Draw("PE");
  jetE[1][0]->Draw("SAME P E");
  TLegend* jetEleg = new TLegend(0.7,0.7,0.9,0.9);
  jetEleg->AddEntry(jetE[0][0],"Minbias Pythia","p");
  jetEleg->AddEntry(jetE[1][0],"Minbias Data","p");
  jetEleg->SetFillStyle(0);
  jetEleg->SetBorderSize(0);
  jetEleg->SetFillColor(0);
  jetEleg->Draw();
  d->SaveAs("output/rmg/summed_jetE.pdf");
  //jetEleg->AddEntry(jetTrigE[0],"Jet 4 GeV Trigger","p");
  //jetEleg->AddEntry(jetTrigE[1],"Jet 6 GeV Trigger","p");
  //jetEleg->AddEntry(jetTrigE[2],"Jet 8 GeV Trigger","p");
  //jetEleg->AddEntry(jetTrigE[3],"Jet 10 GeV Trigger","p");
  for(int h=0; h<2; ++h)
    {
      for(int i=0; i<4; ++i)
	{
	  h1_rej[h][i]->Scale(1./tree->GetEntries());
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
  TLegend* leg = new TLegend(0.12,0.6,0.55,0.85);
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

  return 0;
}
