#include <TFile.h>
#include <TH2.h>
#include <TObject.h>
#include <TKey.h>
#include <TList.h>
#include <TString.h>
#include <vector>
#include <iostream>
#include <TCanvas.h>
#include <TH2D.h>
#include <TGaxis.h>
#include "../dlUtility.h"

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

void FormatTH1(TH1* thehist, string xTitle, string yTitle, int histColor, float yMin, float yMax)
{
  thehist->GetYaxis()->SetTitleOffset(2.0);
  thehist->SetMarkerColor(histColor);
  thehist->SetMarkerStyle(20);
  thehist->SetMarkerSize(2);
  thehist->SetLineColor(histColor);
  thehist->Scale(1./thehist->GetBinWidth(1));
  thehist->GetYaxis()->SetRangeUser(yMin, yMax);
  thehist->GetXaxis()->SetTitle(xTitle.c_str());
  thehist->GetYaxis()->SetTitle(yTitle.c_str());
  thehist->GetYaxis()->SetTitleOffset(2.0);
}

void FormatAndDrawHistogram(TCanvas* canvas, TH2* hist, 
                            string fname, const char* xTitle, const char* yTitle, const char* zTitle,
                            float leftMargin, float rightMargin, float topMargin, float bottomMargin,
                            float xTitleOffset, float yTitleOffset, float zTitleOffset,
                            float xTitleSize, float yTitleSize, float zTitleSize,
                            float xLabelSize, float yLabelSize, float zLabelSize,
			    string* texts, TProfile* profx = NULL, TProfile* profy = NULL) {
    
  // Check if the histogram is valid
  if (!hist) {
    std::cerr << "Error: Histogram is null." << std::endl;
    return;
  }

  float zmax = 0;
  float zmin = 100000000;

  for(int i=0; i<hist->GetNbinsX(); ++i)
    {
      for(int j=0; j<hist->GetNbinsY(); ++j)
	{
	  float binc = hist->GetBinContent(i+1,j+1);
	  if(binc > zmax) zmax = binc;
	  if(binc < zmin && binc > 0) zmin = binc;
	}
    }
  hist->GetZaxis()->SetRangeUser(zmin,zmax);
  // Set pad margins
  canvas->SetLeftMargin(leftMargin);
  canvas->SetRightMargin(rightMargin);
  canvas->SetTopMargin(topMargin);
  canvas->SetBottomMargin(bottomMargin);
  gPad->SetLogz();
  // Set axis titles
  hist->SetXTitle(xTitle);
  hist->SetYTitle(yTitle);
  hist->SetZTitle(zTitle);

  // Set title offsets
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
  hist->GetZaxis()->SetTitleOffset(zTitleOffset);

  // Set title sizes
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetZaxis()->SetTitleSize(zTitleSize);

  TGaxis::SetExponentOffset(0,-0.05,"x");

  // Set label sizes
  hist->GetXaxis()->SetLabelSize(xLabelSize);
  hist->GetYaxis()->SetLabelSize(yLabelSize);
  hist->GetZaxis()->SetLabelSize(zLabelSize);

  // Draw the histogram
  hist->Draw("COLZ"); // Use "COLZ" to draw with color palet
  // Update the canvas to display changes
  std_text(canvas,texts,1,0.025,0.25,0.98,0);
  TLine* hiemcut = new TLine(0.9,25,0.9,100);
  TLine* cETcut = new TLine(0.9,25,1.25,25);
  TLine* loemcut = new TLine(0.1,25,0.1,100);
  TLine* dETcut = new TLine(-.25,7.5,0.1,25);
  //hiemcut->Draw();
  //cETcut->Draw();
  //loemcut->Draw();
  //dETcut->Draw();
  canvas->SaveAs((fname + ".png").c_str());
  if(profx)
    {
      profx->SetMarkerStyle(20);
      profx->SetMarkerSize(1);
      profx->SetMarkerColor(kRed);
      profx->SetLineColor(kRed);
      profx->Draw("SAME PS");
      canvas->SaveAs((fname + "_profx.png").c_str());
    }
  if(profy)
    {
      profy->SetMarkerStyle(20);
      profy->SetMarkerSize(1);
      profy->SetMarkerColor(kRed);
      profy->SetLineColor(kRed);
      profy->Draw("PS");
      std_text(canvas,texts,1,0.025,0.25,0.98,0);
      canvas->SaveAs((fname + "_profy.png").c_str());
    }
  canvas->Update();
}


void draw_dancheck(const TString& fileName, const TString& simFileName) {
  // Open the ROOT file
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  const int nhist = 92;
  TFile* file = TFile::Open(fileName);
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << fileName << std::endl;
    return;
  }

  TFile* simfile = TFile::Open(simFileName);
  if (!simfile || simfile->IsZombie()) {
    std::cerr << "Error opening file: " << simFileName << std::endl;
    return;
  }

  // Vector to hold pointers to TH2 histograms
  std::vector<TH2*> histograms;
  TList* keys = file->GetListOfKeys();
  cout << "nKey: " << keys->GetSize() << endl;
  TIter iter(keys);
  TKey* key;
  while ((key = (TKey*)iter())) {
    cout << "key: "<< key << endl;
    TObject* obj = key->ReadObj();
    cout << "obj before: "<<obj << endl;
    if (obj && (obj->IsA()->InheritsFrom("TH2") || obj->IsA()->InheritsFrom("TH1")) ) {
      TString histName = obj->GetName();
      if(histName.BeginsWith("dancheck_"))
	{
	  histograms.push_back((TH2*)obj);
	}
    }
    else cout <<"obj again: "<< obj << endl;
  }


  std::vector<TH2*> hists2_sim;
  TList* simkeys = simfile->GetListOfKeys();
  cout << "nKey: " << keys->GetSize() << endl;
  TIter simiter(simkeys);
  TKey* simkey;
  while ((simkey = (TKey*)simiter())) {
    cout << "key: "<< simkey << endl;
    TObject* obj = simkey->ReadObj();
    cout << "obj before: "<<obj << endl;
    if (obj && (obj->IsA()->InheritsFrom("TH2") || obj->IsA()->InheritsFrom("TH1")) ) {
      TString histName = obj->GetName();
      if(histName.BeginsWith("dancheck_"))
	{
	  hists2_sim.push_back((TH2*)obj);
	}
    }
    else cout <<"obj again: "<< obj << endl;
  }

  // Output the names of the retrieved histograms
  std::cout << "Retrieved 2D histograms:" << std::endl;
  cout << histograms.size() << endl;
  for (auto* hist : histograms) {
    cout << hist << endl;
    std::cout << hist->GetName() << std::endl;
  }


  TCanvas* canvas = new TCanvas("","",1000,1000);

  string ct[1] = {"Cuts applied"};
  string nct[1] = {"No cuts applied"};

  //histograms.at(0)->RebinY(4);
  //histograms.at(2)->RebinY(4);
  //hists2_sim.at(0)->RebinY(4);
  //hists2_sim.at(2)->RebinY(4);

  TProfile* profx[8];
  TProfile* profy[8];

  for(int i=0; i<4; ++i)
    {
      
      if(i%2==0)
	{
	  for(int j=0; j<histograms.at(i)->GetNbinsX(); ++j)
	    {
	      for(int k=0; k<histograms.at(i)->GetNbinsY()/3; ++k)
		{
		  histograms.at(i)->SetBinContent(histograms.at(i)->GetBin(j,k),0);
		}
	      for(int k=2*histograms.at(i)->GetNbinsY()/3; k<histograms.at(i)->GetNbinsY()+1; ++k)
		{
		  histograms.at(i)->SetBinContent(histograms.at(i)->GetBin(j,k),0);
		}
	    }
	}
      
      profx[i] = histograms.at(i)->ProfileX(("pfxdat"+to_string(i)).c_str(),1,-1,"S");//(i%2==0?100:1),(i%2==0?140:-1),"S");
      profy[i] = histograms.at(i)->ProfileY(("pfydat"+to_string(i)).c_str(),1,-1,"S");
      if(i%2==0)histograms.at(i)->GetYaxis()->SetRangeUser(-0.1,0.1);
      //if(i%2==0) histograms.at(i)->GetXaxis()->SetRangeUser(-0,0.5);
    }

  for(int i=0; i<4; ++i)
    {
     
      if(i%2==0)
	{
	  for(int j=0; j<hists2_sim.at(i)->GetNbinsX(); ++j)
	    {
	      for(int k=0; k<hists2_sim.at(i)->GetNbinsY()/3; ++k)
		{
		  hists2_sim.at(i)->SetBinContent(hists2_sim.at(i)->GetBin(j,k),0);
		}
	      for(int k=2*hists2_sim.at(i)->GetNbinsY()/3; k<hists2_sim.at(i)->GetNbinsY()+1; ++k)
		{
		  hists2_sim.at(i)->SetBinContent(hists2_sim.at(i)->GetBin(j,k),0);
		}
	    }
	}
      
      profx[i+4] = hists2_sim.at(i)->ProfileX(("pfxsim"+to_string(i+4)).c_str(),1,-1,"S");//(i%2==0?100:1),(i%2==0?140:-1),"S");
      profy[i+4] = hists2_sim.at(i)->ProfileY(("pfysim"+to_string(i+4)).c_str(),1,-1,"S");
      if(i%2==0)hists2_sim.at(i)->GetYaxis()->SetRangeUser(-0.1,0.1);
      //if(i%2==0) hists2_sim.at(i)->GetXaxis()->SetRangeUser(-0.5,0.5);
    }
  

  FormatAndDrawHistogram(
			 canvas,histograms.at(0),
			 "../output/chi2img/2pc_nocut_data","#Delta#eta","#Delta#phi","N_{Tower}",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 nct, profx[0], profy[0]
			 );

  FormatAndDrawHistogram(
			 canvas,histograms.at(1),
			 "../output/chi2img/dPhiEmfrac_nocut_data","#Delta#phi","E_{T,EM}/E_{T,Jet}","N_{Tower}",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 nct, profx[1], profy[1]
			 );

  FormatAndDrawHistogram(
			 canvas,histograms.at(2),
			 "../output/chi2img/2pc_cut_data","#Delta#eta","#Delta#phi","N_{Tower}",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 ct, profx[2], profy[2]
			 );

  FormatAndDrawHistogram(
			 canvas,histograms.at(3),
			 "../output/chi2img/dPhiEmfrac_cut_data","#Delta#phi","E_{T,EM}/E_{T,Jet}","N_{Tower}",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 ct, profx[3], profy[3]
			 );

  FormatAndDrawHistogram(
			 canvas,hists2_sim.at(0),
			 "../output/chi2img/2pc_nocut_sim","#Delta#eta","#Delta#phi","N_{Tower}",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 nct, profx[4], profy[4]
			 );

  FormatAndDrawHistogram(
			 canvas,hists2_sim.at(1),
			 "../output/chi2img/dPhiEmfrac_nocut_sim","#Delta#phi","E_{T,EM}/E_{T,Jet}","N_{Tower}",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 nct, profx[5], profy[5]
			 );

  FormatAndDrawHistogram(
			 canvas,hists2_sim.at(2),
			 "../output/chi2img/2pc_cut_sim","#Delta#eta","#Delta#phi","N_{Tower}",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 ct, profx[6], profy[6]
			 );

  FormatAndDrawHistogram(
			 canvas,hists2_sim.at(3),
			 "../output/chi2img/dPhiEmfrac_cut_sim","#Delta#phi","E_{T,EM}/E_{T,Jet}","N_{Tower}",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 ct, profx[7], profy[7]
			 );
  // Clean up
  return 0;
}
