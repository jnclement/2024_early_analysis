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
void std_text(TCanvas* thecan, string* texts, int ntext, float textsize, float textx, float texty, int rightalign)
{

  thecan->cd();
  sphenixtext(textx,texty,0,textsize);
  sqrt_s_text(textx,texty-(6*textsize/4),0,textsize);
  float drawy = texty-3*textsize;
  for(int i=0; i<ntext; ++i)
    {
      drawText(texts[i].c_str(),textx,drawy,rightalign,kBlack,textsize);
      drawy -= 6*textsize/4;
    }
}

void std_hist(TH1D* dathist=NULL, TH1D* simhist=NULL, string xtitle="", string ytitle="")
{
  if(dathist)
    {
      dathist->SetMarkerSize(2);
      dathist->SetMarkerStyle(20);
      dathist->SetLineWidth(2);
      dathist->SetMarkerColor(kMagenta+1);
      dathist->SetLineColor(kMagenta+1);
      dathist->GetXaxis()->SetTitleSize(0.03);
      dathist->GetXaxis()->SetLabelSize(0.03);
      dathist->GetYaxis()->SetTitleSize(0.03);
      dathist->GetYaxis()->SetLabelSize(0.03);
      dathist->GetYaxis()->SetTitleOffset(1.85);
      dathist->GetYaxis()->SetTitle(ytitle.c_str());
      dathist->GetXaxis()->SetTitle(xtitle.c_str());
    }
  if(simhist)
    {
      simhist->SetLineWidth(2);
      simhist->SetMarkerStyle(25);
      simhist->SetMarkerColor(kGreen+2);
      simhist->SetLineColor(kGreen+2);
      simhist->GetXaxis()->SetTitle(xtitle.c_str());
      simhist->GetYaxis()->SetTitle(ytitle.c_str());
      simhist->SetMarkerSize(2);
      simhist->GetXaxis()->SetTitleSize(0.03);
      simhist->GetXaxis()->SetLabelSize(0.03);
      simhist->GetYaxis()->SetTitleSize(0.03);
      simhist->GetYaxis()->SetLabelSize(0.03);
      simhist->GetYaxis()->SetTitleOffset(1.85);
    }
}

float max(float a, float b)
{
  if(a>b) return a;
  return b;
}
double mygaus(double *x, double* par)
{
  float meanterm = (x[0]-par[1])/par[2];
  float expterm = -0.5*meanterm*meanterm;
  float fullgaus = par[0]*exp(expterm);
  return (fullgaus+par[3]);
}
int plot_olay(string filebase = "summed_dat.root", string sfilebase="summed_sim.root", string filelist="lists/sumdatlist.list", string rmgdir = "rmg", int therun = -1)
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
  const int nruns = 94;//286;//33;//72//10;//94
  const int nsd = 64;
  int sds[nruns][nsd] = {0};

  gROOT->SetStyle("Plain");
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0);

  int markerstyles[8] = {20,20,20,20,24,24,24,24};//{20,20,21,43,24,24,25,42};
  int markercolors[4] = {kMagenta+1,kGreen+1,kMagenta+1,kGreen+1};
  TFile* olay = new TFile("output/outolay.root");
  TFile* files[2];
  files[0] = new TFile("summed_dat2.root");
  files[1] = new TFile("summed_dat.root");
  TH1D* jets[8];
  TH1D* dphi[2][4];
  TH1D* jfrc[2];
  jets[0] = (TH1D*)olay->Get("jetE10rj");
  jets[4] = (TH1D*)olay->Get("jetE10nr");
  jets[1] = (TH1D*)olay->Get("jetTrigE1rj");
  jets[2] = (TH1D*)olay->Get("jetTrigE2rj");
  jets[3] = (TH1D*)olay->Get("jetTrigE3rj");
  jets[5] = (TH1D*)olay->Get("jetTrigE1nr");
  jets[6] = (TH1D*)olay->Get("jetTrigE2nr");
  jets[7] = (TH1D*)olay->Get("jetTrigE3nr");

  TCanvas* c = new TCanvas("","",1000,1000);
  c->cd();
  for(int h=0; h<2; ++h)
    {
      jfrc[h] = (TH1D*)files[h]->Get(("hejet"+to_string(1)).c_str());
      for(int i=0; i<4; ++i)
	{
	  dphi[h][i] = (TH1D*)files[h]->Get(("h1_dphi"+to_string(1)+"_"+to_string(i)).c_str());
	  
	  //dphi[h][i]->GetXaxis()->SetTitle("#Delta #phi [rad]");
	  //dphi[h][i]->GetYaxis()->SetTitle("Survival Ratio After BG Rejection");
	  //jfrc[h][i]->GetXaxis()->SetTitle("EM Fraction");
	  //jfrc[h][i]->GetYaxis()->SetTitle("Survival Ratio After BG Rejection");
	  std_hist(dphi[h][i],NULL,"#Delta#phi [rad]","Survival Ratio After BG Rejection");
	  std_hist(jfrc[h][i],NULL,"EM Fraction","Survival Ratio After BG Rejection");
	}
    }
  jfrc[1]
  for(int i=0; i<4; ++i)
    {
      cout << i << endl;
      //jets[i]->Rebin(50);
      //jets[i+4]->Rebin(50);
      std_hist(jets[i],NULL,"Raw #it{E}_{T,jet} [GeV]","Survival Fraction after BG Rejection");
      cout << "dividing dphi" << endl;
      dphi[1][i]->Divide(dphi[0][i]);
      cout << "dividing jfrc" << endl;
      if(i<3)jfrc[1][i]->Divide(jfrc[0][i]);
      cout << "dividing jets" << endl;
      jets[i]->Divide(jets[i+4]);
      jets[i]->Draw("PE");
      c->SaveAs(("output/pTsurfrac"+to_string(i)+".png").c_str());
      dphi[1][i]->Draw("PE");
      c->SaveAs(("output/dphisurfrac"+to_string(i)+".png").c_str());
      if(i<3)jfrc[1][i]->Draw("PE");
      c->SaveAs(("output/efrcsurfrac"+to_string(i)+".png").c_str());
    }
  

  //gPad->SetLogy();
  /*
  for(int i=0; i<8; i++)
    {
      jets[i]->GetYaxis()->SetRangeUser(1e-11,1e-2);
      jets[i]->SetMarkerStyle(markerstyles[i]);
      jets[i]->SetMarkerColor(markercolors[i%4]+2*(i/4));
      jets[i]->SetLineColor(markercolors[i%4]+2*(i/4));
      jets[i]->SetLineWidth(2);
      jets[i]->SetMarkerSize(2);
      jets[i]->GetXaxis()->SetTitleSize(0.03);
      jets[i]->GetXaxis()->SetLabelSize(0.03);
      jets[i]->GetYaxis()->SetTitleSize(0.03);
      jets[i]->GetYaxis()->SetLabelSize(0.03);
      jets[i]->GetYaxis()->SetTitleOffset(1.85);
      jets[i]->GetXaxis()->SetRangeUser(10,50);
      jets[i]->GetXaxis()->SetTitle("Raw #it{E}_{T,jet} [GeV]");
      jets[i]->GetYaxis()->SetTitle("#frac{1}{N_{evt}} #frac{dN_{evt}}{d#it{E}_{T,jet}}");
    }
  
  string jettrignum[4] = {"","\& Jet 8 GeV ","\& Jet 10 GeV ","\& Jet 12 GeV "};
  for(int i=0; i<4; ++i)
    {
      jets[i]->Draw("PE");
      jets[i+4]->Draw("SAME PE");
      
      string texts[6] = {
	"\\mathscr{L}_\\text{MB} = 3.6 \\text{ nb}^{-1}",
	"\\mathscr{L}_\\text{jet8} = 714 \\text{ nb}^{-1}",
	"Calorimeter Anti-k_{#kern[-1.0]{T}} R=0.4",
	"|z_{vtx}|<100 cm",
	"Open points without streak rejection",
	"Filled points with rejection"	
      };
      TLegend* jetleg = new TLegend(0.43,0.5,0.8,0.64);
      jetleg->SetTextFont(42);
      jetleg->SetTextSize(0.025);
      jetleg->SetFillStyle(0);
      jetleg->SetBorderSize(0);
      jetleg->SetFillColor(0);
      jetleg->AddEntry(jets[i],("MB "+jettrignum[i]+"Trig. w/ Streak Rej.").c_str(),"p");
      jetleg->AddEntry(jets[i+4],("MB "+jettrignum[i]+"Trig. w/o Streak Rej.").c_str(),"p");
      jetleg->Draw();
      
      std_text(c, texts, 6, 0.025, 0.46, 0.91, 0);

      c->SaveAs(("output/overlayplot"+to_string(i)+".png").c_str());
    }
  */
  return 0;
}
