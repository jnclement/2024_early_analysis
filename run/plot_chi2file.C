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
#include <TDatime.h>

int plot_chi2file()
{
  TCanvas* c = new TCanvas("","",1000,1000);
  gROOT->ProcessLine( "gErrorIgnoreLevel = 1001;");
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptTitle(0);
  c->cd();
  TFile* thefile = new TFile("summed_chi2file.root");
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  const int nh = 3;
  TH2D* hists[nh][6];
  TH2D* h_ecc_E[nh];
  TH2D* h2_g20_ecc_angle[nh];
  TH2D* h2_g20_ecc_frcoh[nh];
  TH2D* h2_g20_ecc_frcem[nh];
  gPad->SetLogz();
  for(int h=0; h<3; ++h)
    {
      h2_g20_ecc_angle[h] = (TH2D*)thefile->Get(("h2_ecc_angle"+to_string(h)).c_str());
      h2_g20_ecc_angle[h]->GetXaxis()->SetTitle("Leading Jet Eccentricity");
      h2_g20_ecc_angle[h]->GetYaxis()->SetTitle("Orientation Angle [rad]");
      h2_g20_ecc_angle[h]->GetYaxis()->SetTitleOffset(1.5);
      h2_g20_ecc_angle[h]->GetZaxis()->SetTitleOffset(1.5);
      h2_g20_ecc_angle[h]->GetZaxis()->SetTitle("Counts");
      h2_g20_ecc_angle[h]->Draw("COLZ");
      sphenixtext();
      c->SaveAs(("output/rmg/h2_g20_ecc_angle"+to_string(h)).c_str());
      
      h2_g20_ecc_frcoh[h] = (TH2D*)thefile->Get(("h2_ecc_frcoh"+to_string(h)).c_str());
      h2_g20_ecc_frcoh[h]->GetXaxis()->SetTitle("Leading Jet Eccentricity");
      h2_g20_ecc_frcoh[h]->GetYaxis()->SetTitle("E_{T,lead jet,OHCal}/E_{T,lead jet}");
      h2_g20_ecc_frcoh[h]->GetYaxis()->SetTitleOffset(1.5);
      h2_g20_ecc_frcoh[h]->GetZaxis()->SetTitleOffset(1.5);
      h2_g20_ecc_frcoh[h]->GetZaxis()->SetTitle("Counts");
      h2_g20_ecc_frcoh[h]->Draw("COLZ");
      sphenixtext();
      c->SaveAs(("output/rmg/h2_g20_ecc_frcoh"+to_string(h)).c_str());
      
      h2_g20_ecc_frcem[h] = (TH2D*)thefile->Get(("h2_ecc_frcem"+to_string(h)).c_str());
      h2_g20_ecc_frcem[h]->GetXaxis()->SetTitle("Leading Jet Eccentricity");
      h2_g20_ecc_frcem[h]->GetYaxis()->SetTitle("E_{T,lead jet,EMCal}/E_{T,lead jet}");
      h2_g20_ecc_frcem[h]->GetYaxis()->SetTitleOffset(1.5);
      h2_g20_ecc_frcem[h]->GetZaxis()->SetTitleOffset(1.5);
      h2_g20_ecc_frcem[h]->GetZaxis()->SetTitle("Counts");
      h2_g20_ecc_frcem[h]->Draw("COLZ");
      sphenixtext();
      c->SaveAs(("output/rmg/h2_g20_ecc_frcem"+to_string(h)).c_str());
      

      h_ecc_E[h] = (TH2D*)thefile->Get(("h2_ecc_E"+to_string(h)).c_str());
      h_ecc_E[h]->GetXaxis()->SetTitle("Leading Jet Eccentricity");
      h_ecc_E[h]->GetXaxis()->SetTitle("E_{T,lead jet}");
      h_ecc_E[h]->GetYaxis()->SetTitleOffset(1.5);
      h_ecc_E[h]->GetZaxis()->SetTitleOffset(1.5);
      h_ecc_E[h]->GetZaxis()->SetTitle("Counts");
      h_ecc_E[h]->Draw("COLZ");
      sphenixtext();
      c->SaveAs(("output/rmg/h_ecc_E.png"+to_string(h)).c_str());
     
      float ETJet[6] = {4,10,20,30,40,50};
      for(int i=0; i<6; ++i)
	{
	  hists[h][i] = (TH2D*)thefile->Get(("h2_ecc_layer"+to_string(i)).c_str());
	  hists[h][i]->GetXaxis()->SetTitle("Leading Jet Eccentricity");
	  hists[h][i]->GetYaxis()->SetTitle("Max Layer E_{T} / E_{T,lead jet}");
	  hists[h][i]->GetYaxis()->SetTitleOffset(1.5);
	  hists[h][i]->GetZaxis()->SetTitleOffset(1.5);
	  hists[h][i]->GetZaxis()->SetTitle("Counts");
	  hists[h][i]->Draw("COLZ");
	  sphenixtext();
	  if(i<5)
	    {
	      
	      stringstream stream[2];
	      stream[0] << std::fixed << std::setprecision(0) << ETJet[i];
	      stream[1] << std::fixed << std::setprecision(0) << ETJet[i+1];
	      string text = stream[0].str() + " < E_{T,jet} < " + stream[1].str();
	      drawText(text.c_str(), 0.15, 0.93);
	    }
	  else
	    {
	  stringstream stream;
	  stream << std::fixed << std::setprecision(0) << ETJet[i];
	  string text = stream.str() + " < E_{T,jet}";
	  drawText(text.c_str(), 0.15, 0.93);
	    } 
	  c->SaveAs(("output/rmg/chi2file"+to_string(h)+"_"+to_string(i)+".png").c_str());
	}
    }
  return 0;
}
