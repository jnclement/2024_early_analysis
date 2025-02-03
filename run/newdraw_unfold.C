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
  //thehist->Scale(1./thehist->GetBinWidth(1));
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
  std_text(canvas,texts,2,0.025,0.25,0.97,0);
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


void newdraw_unfold(const TString& fileName, const TString& simFileName, const TString& mbFileName, const TString& jet30filename) {
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
  

  TH1F* datnocut;
  TH1F* datallcut;
  TH1F* datspeccut;
  TH1F* j10allcut;
  TH1F* j30allcut;
  TH1F* mballcut;
  TH1F* j10nocut;
  TH1F* j30nocut;
  TH1F* mbnocut;
  TH1F* j10tp;
  TH1F* j30tp; 
  TH1F* mbtp;
  TH2F* j30r;
  TH2F* j10r;
  TH2F* mbr;
  datnocut = static_cast<TH1F*>(file->Get("h1_jetSpectra_0"));
  datallcut = static_cast<TH1F*>(file->Get("h1_jetSpectra_31"));
  datspeccut = static_cast<TH1F*>(file->Get("h1_jetSpectra_30"));
  j10allcut = static_cast<TH1F*>(simfile->Get("h1_cspec"));
  j10nocut = static_cast<TH1F*>(simfile->Get("h1_ucspec"));
  j10tp = static_cast<TH1F*>(simfile->Get("h1_tjetspec"));
  mballcut = static_cast<TH1F*>(mbfile->Get("h1_cspec"));
  mbnocut = static_cast<TH1F*>(mbfile->Get("h1_ucspec"));
  mbtp = static_cast<TH1F*>(mbfile->Get("h1_tjetspec"));
  j30allcut = static_cast<TH1F*>(j30file->Get("h1_cspec"));
  j30nocut = static_cast<TH1F*>(j30file->Get("h1_ucspec"));
  j30tp = static_cast<TH1F*>(j30file->Get("h1_tjetspec"));
  j10r = static_cast<TH2F*>(simfile->Get("hresponse"));
  mbr = static_cast<TH2F*>(mbfile->Get("hresponse"));
  j30r = static_cast<TH2F*>(j30file->Get("hresponse"));



  RooUnfoldResponse* simresp = static_cast<RooUnfoldResponse*>(mbfile->Get("roounfold_response"));
  RooUnfoldResponse* j10resp = static_cast<RooUnfoldResponse*>(simfile->Get("roounfold_response"));
  RooUnfoldResponse* j30resp = static_cast<RooUnfoldResponse*>(j30file->Get("roounfold_response"));
  
  simresp->Add(*j10resp);
  simresp->Add(*j30resp);

  TF1* fturn = new TF1("fturn","[0]*(TMath::Erf((x-[1])/[2]) + 1)", 4, 30);
  fturn->SetParameter(0,0.478683);
  fturn->SetParameter(1,9.76439);
  fturn->SetParameter(2,4.29056);
  double jet10scale = 3.646e-6 / 4.197e-2;
  double jet30scale = 2.505e-9 / 4.197e-2;
  jet30scale /= 2.505e-9 / 3.646e-6;
  for(int i=1; i<5; ++i)
    {
      datnocut->SetBinContent(i,datnocut->GetBinContent(i)/fturn->Eval(datnocut->GetBinCenter(i)));
      datallcut->SetBinContent(i,datallcut->GetBinContent(i)/fturn->Eval(datallcut->GetBinCenter(i)));
      datspeccut->SetBinContent(i,datspeccut->GetBinContent(i)/fturn->Eval(datspeccut->GetBinCenter(i)));

      datnocut->SetBinError(i,datnocut->GetBinError(i)/fturn->Eval(datnocut->GetBinCenter(i)));
      datallcut->SetBinError(i,datallcut->GetBinError(i)/fturn->Eval(datallcut->GetBinCenter(i)));
      datspeccut->SetBinError(i,datspeccut->GetBinError(i)/fturn->Eval(datspeccut->GetBinCenter(i)));

    }
  
  j10allcut->Scale(1./jet10scale);
  j10nocut->Scale(1./jet10scale);
  j10tp->Scale(1./jet10scale);
  j30allcut->Scale(1./jet30scale);
  j30nocut->Scale(1./jet30scale);
  j30tp->Scale(1./jet30scale);
  j10r->Scale(1./jet10scale);
  j30r->Scale(1./jet30scale);
  /*
  float xpoints[24] = {14.8, 17.3, 19.9, 22.2, 24.9, 27.4, 29.8, 32.5, 34.8,37.4,39.9,42.4,44.8,47.5,50,52.3,54.9,57.4,59.9,62.3,67.5,69.9,72.4};
  float ypoints[24] = {0.0000068815953765315700,
		       0.00000262743707983794,
		       0.000001131475401952760,
		       4.87256800571268E-07,
		       2.36668415191766E-07,
		       1.08240027881874E-07,
		       5.25739113935031E-08,
		       2.71198624697299E-08,
		       1.3172550620955E-08,
		       6.79496259165534E-09,
		       3.30042192529477E-09,
		       1.70249818462502E-09,
		       8.78221068172245E-10,
		       4.26566008208939E-10,
		       2.2004091326402E-10,
		       1.00635458971372E-10,
		       4.88802507403217E-11,
		       2.37419189702947E-11,
		       1.08583393742322E-11,
		       4.67601831544398E-12,
		       2.01367322688913E-12,
		       8.16520302177809E-13,
		       3.11752667949471E-13};
  */
  //TGraph* datathief = new TGraph(24,xpoints,ypoints);
  float ppScaleFactor = 28*(2.0*M_PI*2.5*2.0)/(1.e12);
  ifstream tjetfile;
  tjetfile.open("jets_newphenix_sc1.dat");
  TGraph* tjet = new TGraph();
  tjet->SetMarkerStyle(20);
  for (int index=0; index<38; index++) {
    double pt; double yield;
    tjetfile >> pt >> yield;
    // rescale points 
    yield = yield * ppScaleFactor*pt / 2.5;
    tjet->SetPoint(index,pt,yield);
  }
  /*
  TH1F* comball = mballcut;
  TH1F* combno = mbnocut;
  TH1F* combtp = mbtp;
  TH2F* combr = mbr;
  */
  TH1F* comball = j10allcut;
  TH1F* combno = j10nocut;
  TH1F* combtp = j10tp;
  TH2F* combr = j10r;

  comball->Add(j10allcut);
  comball->Add(j30allcut);
  combno->Add(j10nocut);
  combno->Add(j30nocut);
  combtp->Add(j10tp);
  combtp->Add(j30tp);
  combr->Add(j10r);
  combr->Add(j30r);


  const int nbinx = 10;
  const int nbiny = 9;
  float binsx[nbinx+1] = {8,12,16,21,26,32,38,45,52,60,70};
  float binsy[nbiny+1] = {12,16,21,26,32,38,45,52,60,70};
  
  TH1F* tempcomball = comball;
  TH1F* tempcombtp = combtp;
  TH1F* tempdatall = datallcut;
  TH2F* tempr = combr;
  
  
  //TH1F* tempcomball = new TH1F("tempcomball","",nbiny,binsy);
  //TH1F* tempcombtp = new TH1F("tempcombtp","",nbinx,binsx);
  //TH1F* tempdatall = new TH1F("tempdatall","",nbiny,binsy);
  //TH2F* tempr = new TH2F("tempr","",nbiny,binsy,nbinx,binsx);
  /*
  for(int i=3; i<nbinx+3; ++i)
    {
      tempcombtp->SetBinContent(i-2,combtp->GetBinContent(i));
      tempcombtp->SetBinError(i-2,combtp->GetBinError(i));

      for(int j=3; j<nbiny+3; ++j)
	{
	  if(i==3)
	    {
	      tempdatall->SetBinContent(j-2,datallcut->GetBinContent(j+1));
	      tempdatall->SetBinError(j-2,datallcut->GetBinError(j+1));
	      tempcomball->SetBinContent(j-2,comball->GetBinContent(j));
	      tempcomball->SetBinError(j-2,comball->GetBinError(j));
	    }
	  tempr->SetBinContent(j-2,i-2,combr->GetBinContent(j,i));
	  tempr->SetBinError(j-2,i-2,combr->GetBinError(j,i));
	}
    }
  */
  //cout <<"hresp: " << hResp << endl;
  //  RooUnfoldResponse response(comball,combtp,combr);

  //RooUnfoldBayes unfold(&response, datallcut, 2);
  cout << "here1" << endl;
  //RooUnfoldResponse response(tempcomball,tempcombtp,tempr);
  cout << "here2" << endl;
  cout << tempdatall->GetNbinsX() << endl;
  cout << tempdatall->GetBinLowEdge(1) << endl;
  cout << tempr->GetNbinsX() << endl;
  RooUnfoldBayes unfold(simresp, tempcomball, 1, false, true);
  TH1D* hUnfold = (TH1D*) unfold.Hunfold(RooUnfold::kErrors);
  cout <<"HUNFOLD PARAMS:" << endl << hUnfold->GetNbinsX() << endl << hUnfold->GetXaxis()->GetBinWidth(1) << endl;
  for(int i=1; i<nbiny+1; ++i)
    {
      cout << hUnfold->GetBinContent(i) << endl;
    }

  tempcomball->Scale(4e-5/2.8e6);
  combno->Scale(4e-5/2.8e6);
  tempcombtp->Scale(4e-5/2.8e6);
  tempdatall->Scale(1./1.79769e11);
  /*
  FormatTH1(comball, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kViolet+2,1e-12,1e-4);
  FormatTH1(combno, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kRed+2,1e-12,1e-4);
  FormatTH1(datnocut, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kBlack,1e-12,1e-4);
  FormatTH1(combtp, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kRed,1e-12,1e-4);
  FormatTH1(datallcut, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kGreen+2,1e-12,1e-4);
  FormatTH1(datspeccut, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kBlue,1e-12,1e-4);
  FormatTH1(hUnfold, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kGreen,1e-12,1e-4);
  */
  //hUnfold->Scale(1./1.79769e11);
  hUnfold->Scale(4e-5/2.8e6);
  FormatTH1(tempcomball, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kViolet+2,1e-13,1e-5);
  FormatTH1(tempdatall, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kGreen+2,1e-13,1e-5);
  FormatTH1(hUnfold, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kGreen,1e-12,1e-5);
  FormatTH1(tempcombtp, "E_{T,jet} [GeV]","#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",kRed,1e-13,1e-5);
  TLegend* theleg = new TLegend(0.2,0.23,0.5,0.4);
  theleg->SetFillStyle(0);
  theleg->SetBorderSize(0);
  theleg->AddEntry(tempcomball,"Reco sim all cuts","p");
  //theleg->AddEntry(combno,"Reco sim no cuts","p");
  theleg->AddEntry(tempcombtp,"Truth PYTHIA","p");
  //theleg->AddEntry(datnocut,"Data no cuts","p");
  theleg->AddEntry(tempdatall,"Data all cuts","p");
  //theleg->AddEntry(datspeccut,"Data all cuts (no dijet check)","p");
  theleg->AddEntry(hUnfold,"Unfolded","p");
  theleg->AddEntry(tjet,"NLO PQCD Calc.","p");
  TCanvas* c1 = new TCanvas("","",1000,1000);
  c1->cd();
  c1->SetLeftMargin(0.2);
  c1->SetBottomMargin(0.2);
  c1->SetLogy();
  for(int i=1; i<tempcombtp->GetNbinsX()+1; ++i)
    {
      tempcombtp->SetBinContent(i,tempcombtp->GetBinContent(i)/tempcombtp->GetBinWidth(i));
      tempcomball->SetBinContent(i,tempcomball->GetBinContent(i)/tempcomball->GetBinWidth(i));
      tempdatall->SetBinContent(i,tempdatall->GetBinContent(i)/tempdatall->GetBinWidth(i));
      hUnfold->SetBinContent(i,hUnfold->GetBinContent(i)/hUnfold->GetBinWidth(i));

      tempcombtp->SetBinError(i,tempcombtp->GetBinError(i)/tempcombtp->GetBinWidth(i));
      tempcomball->SetBinError(i,tempcomball->GetBinError(i)/tempcomball->GetBinWidth(i));
      tempdatall->SetBinError(i,tempdatall->GetBinError(i)/tempdatall->GetBinWidth(i));
      hUnfold->SetBinError(i,hUnfold->GetBinError(i)/hUnfold->GetBinWidth(i));
    }
  tempcombtp->SetMarkerStyle(21);
  tempcomball->SetMarkerStyle(21);
  tempcombtp->Draw("PE");
  tempcomball->Draw("SAME PE");
  //combno->Draw("SAME PE");
  //datnocut->Draw("SAME PE");
  tempdatall->Draw("SAME PE");
  //datspeccut->Draw("SAME PE");
  
  hUnfold->Draw("SAME PE");
  tjet->Draw("SAME P");
  theleg->Draw();
  const int ntext = 5;
  string texts[ntext] =
    {
      "Calorimeter Anti-k_{T} R=0.4",
      "|z_{vtx}| < 30 cm",
      "Jet-8 & MBDNS>=1 Triggered Data",
      "PYTHIA8 Jet-10&30 Sample Sim",
      "\\mathscr{L}_{\\text{data}}=16.17 \\text{pb}^{-1}"
    };
  std_text(c1, texts, ntext, 0.025, 0.5, 0.87, 0);
  c1->SaveAs("output/chi2img/final_unfold.png");
  c1->SetLogy(0);
  //const int nbinx = 11;
  //float binsx[nbinx+1] = {10,12,14,17,20,24,28,33,38,44,50,70};
  TH1F* cdattprat = new TH1F("cdattprat","",nbinx,binsx);
  TH1F* csimtprat = new TH1F("csimtprat","",nbinx,binsx);
  TH1F* cdcsrat = new TH1F("cdcsrat","",nbinx,binsx);
  TH1F* uftprat = new TH1F("uftprat","",nbinx,binsx);
  
      

  for(int i=2; i<nbinx+1; ++i)
    {
      if(comball->GetBinContent(i-1) != 0)
	{
	  if(tempcomball->GetBinContent(i-1) != 0) cdcsrat->SetBinContent(i,tempdatall->GetBinContent(i-1)/tempcomball->GetBinContent(i-1));
	  float cderr = tempdatall->GetBinError(i-1) / tempdatall->GetBinContent(i-1);
	  float cserr = tempcomball->GetBinError(i-1) / tempcomball->GetBinContent(i-1);
	  float tperr = tempcombtp->GetBinError(i) / tempcombtp->GetBinContent(i);
	  float uferr = hUnfold->GetBinError(i) / hUnfold->GetBinContent(i);
	  cdcsrat->SetBinError(i,cdcsrat->GetBinContent(i)*sqrt(cderr*cderr+cserr*cserr));
	  cdattprat->SetBinContent(i,tempdatall->GetBinContent(i-1)/tempcombtp->GetBinContent(i));
	  cdattprat->SetBinError(i,cdattprat->GetBinContent(i)*sqrt(cderr*cderr+tperr*tperr));
	  csimtprat->SetBinContent(i,tempcomball->GetBinContent(i-1)/tempcombtp->GetBinContent(i));
	  csimtprat->SetBinError(i,csimtprat->GetBinContent(i)*sqrt(tperr*tperr+cserr*cserr));
	  uftprat->SetBinContent(i,hUnfold->GetBinContent(i)/tempcombtp->GetBinContent(i));
	  uftprat->SetBinError(i,uftprat->GetBinContent(i)*sqrt(tperr*tperr+uferr*uferr));
	}
      else
	{
	  cdcsrat->SetBinContent(i,0);
	  cdattprat->SetBinContent(i,0);
	  csimtprat->SetBinContent(i,0);
	}
    }
  FormatTH1(cdattprat, "E_{T,jet} [GeV]","Ratio Data/Truth PYTHIA",kBlack,1e-3,1);
  FormatTH1(csimtprat, "E_{T,jet} [GeV]","Ratio Reco Sim/Truth PYTHIA",kBlack,1e-3,1);
  FormatTH1(cdcsrat, "E_{T,jet} [GeV]","Ratio Data/Reco Sim",kBlack,0,2);

  FormatTH1(uftprat, "E_{T,jet} [GeV]","Ratio Unfolded/Truth PYTHIA",kBlack,0,2.5);
  cdattprat->GetYaxis()->SetRangeUser(0,1);
  csimtprat->GetYaxis()->SetRangeUser(0,1);
  uftprat->GetYaxis()->SetRangeUser(0.9,1.1);
  cdattprat->Draw("PE");
  c1->SaveAs("output/chi2img/cdattprat.png");
  csimtprat->Draw("PE");
  c1->SaveAs("output/chi2img/csimtprat.png");
  //cdcsrat->GetXaxis()->SetRangeUser(0,100);
  cdcsrat->Draw("PE");
  c1->SaveAs("output/chi2img/cdcsrat.png");
  uftprat->Draw("PE");
  c1->SaveAs("output/chi2img/uftprat.png");
  c1->SetLogz();
  c1->SetRightMargin(0.2);
  c1->SetTopMargin(0.2);
  tempr->GetXaxis()->SetTitle("Reco E_{T} [GeV]");
  tempr->GetYaxis()->SetTitle("Truth E_{T} [GeV]");
  tempr->Draw("COLZ");
  std_text(c1, texts, ntext, 0.015, 0.5, 0.98, 0);
  c1->SaveAs("output/chi2img/response.png");
  return 0;
}
