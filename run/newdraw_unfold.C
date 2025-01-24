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
#include "dlUtility.h"
#include <RooUnfold.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>
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
                            float xLabelSize, float yLabelSize, float zLabelSize, string* texts) {
    
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
  std_text(canvas,texts,2,0.025,0.25,0.98,0);
  TLine* hiemcut = new TLine(0.9,25,0.9,100);
  TLine* cETcut = new TLine(0.9,25,1.25,25);
  TLine* loemcut = new TLine(0.1,25,0.1,100);
  TLine* dETcut = new TLine(-.25,7.5,0.1,25);
  //hiemcut->Draw();
  //cETcut->Draw();
  //loemcut->Draw();
  //dETcut->Draw();
  canvas->SaveAs((fname + ".png").c_str());
  canvas->Update();
}


void draw_chi2hists(const TString& fileName, const TString& simFileName, const TString& mbFileName, const TString& jet30filename) {
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

  TFile* mbfile = TFile::Open(mbFileName);
  TFile* j30file = TFile::Open(jet30filename);

  // Vector to hold pointers to TH2 histograms
  

  TH1F* datnocut, datallcut, datspeccut, j10allcut, j30allcut, mballcut, j10nocut, j30nocut, mbnocut, j10tp, j30tp, mbtp;
  TH2F* j30r, j10r, mbr;
  fileName->Get("h1_jetSpectra_0",datnocut);
  fileName->Get("h1_jetSpectra_31",datallcut);
  fileName->Get("h1_jetSpectra_30",datspeccut);
  simFileName->Get("h1_cspec",j10allcut);
  simFileName->Get("h1_ucspec",j10nocut);
  simFileName->Get("h1_tjetspec",j10tp);
  mbFileName->Get("h1_cspec",mballcut);
  mbFileName->Get("h1_ucspec",mbnocut);
  mbFileName->Get("h1_tjetspec",mbtp);
  jet30filename->Get("h1_cspec",j30allcut);
  jet30filename->Get("h1_ucspec",j30nocut);
  jet30filename->Get("h1_tjetspec",j30tp);
  simFileName->Get("hresponse",j10r);
  mbFileName->Get("hresponse",mbr);
  jet30filename->Get("hresponse",j30r);
  TF1 fturn = new TF1("fturn","[0]*(TMath::Erf((x-[1])/[2]) + 1)", 4, 30);
  fturn->SetParameter(0,0.478683);
  fturn->SetParameter(1,9.76439);
  fturn->SetParameter(2,4.29056);
  double jet10scale = 3.646e-6 / 4.197e-2;
  double jet30scale = 2.505e-9 / 4.197e-2;

  for(int i=1; i<5; ++i)
    {
      datnocut->SetBinContent(i,datnocut->GetBinContent(i)/fturn->Eval(datnocut->GetBinCenter(i)));
      datallcut->SetBinContent(i,datallcut->GetBinContent(i)/fturn->Eval(datallcut->GetBinCenter(i)));
      datspeccut->SetBinContent(i,datspeccut->GetBinContent(i)/fturn->Eval(datspeccut->GetBinCenter(i)));
    }

  j10allcut->Scale(jet10scale);
  j10nocut->Scale(jet10scale);
  j10tp->Scale(jet10scale);
  j30allcut->Scale(jet30scale);
  j30nocut->Scale(jet30scale);
  j30tp->Scale(jet30scale);
  j10r->Scale(jet10scale);
  j30r->Scale(jet30scale);

  TH1F* comball = mballcut;
  TH1F* combno = mbnocut;
  TH1F* combtp = mbtp;

  TH2F* combr = mbr;

  comball->Add(j10allcut);
  comball->Add(j30allcut);
  combno->Add(j10nocut);
  combno->Add(j30nocut);
  combtp->Add(j10tp);
  combtp->Add(j30tp);
  combr->Add(j10r);
  combr->Add(j30r);

  cout <<"hresp: " << hResp << endl;
  RooUnfoldResponse response(comball,combtp,combr);

  RooUnfoldBayes unfold(&response, datallcut, 2);

  TH1D* hUnfold = (TH1D*) unfold.Hunfold(RooUnfold::kErrors);
  cout <<"HUNFOLD PARAMS:" << endl << hUnfold->GetNbinsX() << endl << hUnfold->GetXaxis()->GetBinWidth(1) << endl;
  for(int i=1; i<11; ++i)
    {
      cout << hUnfold->GetBinContent(i) << endl;
    }

  comball->Scale(1./1e7);
  combno->Scale(1./1e7);
  combtp->Scale(1./1e7);

  FormatTH1(comball, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kViolet+2,1e-12,1e-4);
  FormatTH1(combno, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kRed+2,1e-12,1e-4);
  FormatTH1(datnocut, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kBlack,1e-12,1e-4);
  FormatTH1(combtp, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kRed,1e-12,1e-4);
  FormatTH1(datallcut, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kGreen+2,1e-12,1e-4);
  FormatTH1(datspeccut, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kBlue,1e-12,1e-4);
  FormatTH1(hUnfold, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kGreen,1e-12,1e-4);
  
  TLegend* theleg = new TLegend(0.57,0.55,0.9,0.65);
  theleg->SetFillStyle(0);
  theleg->SetBorderSize(0);
  theleg->AddEntry(comball,"Reco sim all cuts","p");
  theleg->AddEntry(combno,"Reco sim no cuts","p");
  theleg->AddEntry(combtp,"Truth PYTHIA","p");
  theleg->AddEntry(datnocut,"Data no cuts","p");
  theleg->AddEntry(datallcut,"Data all cuts","p");
  theleg->AddEntry(datspeccut,"Data all cuts (no dijet check)","p");
  theleg->AddEntry(hUnfold,"Unfolded","p");
  TCanvas* c1 = new TCanvas("","",1000,1000);
  c1->cd();
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);
  combtp->Draw("PE");
  comball->Draw("SAME PE");
  combno->Draw("SAME PE");
  datnocut->Draw("SAME PE");
  datallcut->Draw("SAME PE");
  datspeccut->Draw("SAME PE");
  hUnfold->Draw("SAME PE");
  theleg->Draw();
  const int ntext = 5;
  string texts[ntext] =
    {
      "Calorimeter Anti-k_{T} R=0.4",
      "|z_{vtx}| < 30 cm",
      "Jet-8 & MBDNS>=1 Triggered Data",
      "PYTHIA8 Jet-10 Sample Sim",
      "\\mathscr{L}_{\\text{data}}=16.17 \\text{pb}^{-1}"
    };
  void std_text(c1, texts, ntext, 0.025, 0.5, 0.8, 0)
  std_text
  c1->SaveAs("output/chi2img/final_unfold.png");

  const int nbinx = 11;
  float binsx[nbinx+1] = {10,12,14,17,20,24,28,33,38,44,50,70};
  TH1F* cdattprat = new TH1F("cdattprat","",nbinx,binsx);
  TH1F* csimtprat = new TH1F("csimtprat","",nbinx,binsx);
  TH1F* cdcsrat = new TH1F("cdcsrat","",nbinx,binsx);

  cdattprat->Divide(datallcut,combtp);
  csimtprat->Divide(comball,combtp);
  for(int i=2; i<nbinx+1; ++i)
    {
      if(comball->GetBinContent(i-1) != 0)
	{
	  cdcsrat->SetBinContent(i,datallcut->GetBinContent(i)/comball->GetBinContent(i-1));
	}
      else
	{
	  cdcsrat->SetBinContent(i,0);
	}
    }
  FormatTH1(cdattprat, "E_{T,jet} [GeV]","Ratio Data/Truth PYTHIA",kBlack,1e-3,1);
  FormatTH1(csimtprat, "E_{T,jet} [GeV]","Ratio Reco Sim/Truth PYTHIA",kBlack,1e-3,1);
  FormatTH1(cdcsrat, "E_{T,jet} [GeV]","Ratio Data/Reco Sim",kBlack,1e-3,1);

  cdattprat->Draw("PE");
  c1->SaveAs("output/chi2img/cdattprat.png");
  csimtprat->Draw("PE");
  c1->SaveAs("output/chi2img/csimtprat.png");
  cdcsrat->Draw("PE");
  c1->SaveAs("output/chi2img/cdcsrat.png")
  return 0;
}
