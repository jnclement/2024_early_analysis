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


void draw_chi2hists(const TString& fileName, const TString& simFileName) {
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
  std::vector<TH1*> jetSpectra;
  std::vector<TH1*> xJ;
  // Iterate through all keys in the file
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
      if(histName.BeginsWith("h2_"))
	{
	  histograms.push_back((TH2*)obj);
	}
      if(histName.BeginsWith("h1_"))
	{
	  cout << "test" << endl;
	  jetSpectra.push_back((TH1*)obj);
	}
      if(histName.BeginsWith("x"))
	{
	  cout << "test" << endl;
	  xJ.push_back((TH1*)obj);
	}
    }
    else cout <<"obj again: "<< obj << endl;
  }


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
      if(histName.BeginsWith("h1_"))
	{
	  cout << "test" << endl;
	  jetSpectra.push_back((TH1*)obj);
	}
      if(histName.BeginsWith("x"))
	{
	  cout << "test" << endl;
	  xJ.push_back((TH1*)obj);
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

  string xtitles[nhist] = {
    "#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower",
    "N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet",
    "E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]",
    "E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]",
    "E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]",
    "#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}",
    "Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}",
    "Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in OHCal",
    "Fraction of E_{T,lead jet} in EMCal","Fraction of E_{T,lead jet} in EMCal","Fraction of E_{T,lead jet} in EMCal","Fraction of E_{T,lead jet} in EMCal","Fraction of E_{T,lead jet} in EMCal",
    "#eta","#eta","#eta","#eta",
    "#phi","#phi","#phi",
    "E_{T,lead jet}","E_{T,lead jet}",
    "E_{T,sublead jet}",
    "A_{J}"
  };

  const string names[nhist] = {
    "maxETowChi2_nBadChi2","maxETowChi2_maxTowDiff","maxETowChi2_chi2", "maxETowChi2_frcoh", "maxETowChi2_frcem", "maxETowChi2_eta", "maxETowChi2_phi", "maxETowChi2_jet_ET", "maxETowChi2_dphi", "maxETowChi2_subjet_ET","maxETowChi2_ecc","maxETowChi2_maxTowE","maxETowChi2_subTowE",
    "nBadChi2_maxTowDiff","nBadChi2_chi2", "nBadChi2_frcoh", "nBadChi2_frcem", "nBadChi2_eta", "nBadChi2_phi", "nBadChi2_jet_ET", "nBadChi2_dphi", "nBadChi2_subjet_ET","nBadChi2_ecc","nBadChi2_maxTowE","nBadChi2_subTowE",
    "maxTowDiff_chi2", "maxTowDiff_frcoh", "maxTowDiff_frcem", "maxTowDiff_eta", "maxTowDiff_phi", "maxTowDiff_jet_ET", "maxTowDiff_dphi", "maxTowDiff_subjet_ET","maxTowDiff_ecc","maxTowDiff_maxTowE","maxTowDiff_subTowE",
    "subTowE_chi2", "subTowE_frcoh", "subTowE_frcem", "subTowE_eta", "subTowE_phi", "subTowE_jet_ET", "subTowE_dphi", "subTowE_subjet_ET","subTowE_ecc","subTowE_maxTowE",
    "maxTowE_chi2", "maxTowE_frcoh", "maxTowE_frcem", "maxTowE_eta", "maxTowE_phi", "maxTowE_jet_ET", "maxTowE_dphi", "maxTowE_subjet_ET","maxTowE_ecc",
    "ecc_chi2", "ecc_frcoh", "ecc_frcem", "ecc_eta", "ecc_phi", "ecc_jet_ET", "ecc_dphi", "ecc_subjet_ET",
    "chi2_frcoh", "chi2_frcem", "chi2_eta", "chi2_phi", "chi2_jet_ET", "chi2_dphi", "chi2_subjet_ET",
    "frcoh_frcem", "frcoh_eta", "frcoh_phi", "frcoh_jet_ET", "frcoh_dphi", "frcoh_subjet_ET",
    "frcem_eta", "frcem_phi", "frcem_jet_ET", "frcem_dphi", "frcem_subjet_ET",
    "eta_phi_h1", "eta_jet_ET", "eta_dphi", "eta_subjet_ET",
    "phi_jet_ET", "phi_dphi", "phi_subjet_ET",
    "jet_ET_dphi", "jet_ET_subjet_ET",
    "subjet_ET_dphi",
    "AJ_dphi"
  };


  string ytitles[nhist] = {
    "N Bad #chi^{2}","E_{T,lead tower} - E_{T,sublead tower} [GeV]","Max Tower #chi^{2}","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in EMCal","#eta","#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}","E_{T,lead tower} [GeV]","E_{T,sublead tower} [GeV]",
    "E_{T,lead tower} - E_{T,sublead tower} [GeV]","Max Tower #chi^{2}","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in EMCal","#eta","#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}","E_{T,lead tower} [GeV]","E_{T,sublead tower} [GeV]",
    "Max Tower #chi^{2}","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in EMCal","#eta","#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}","E_{T,lead tower} [GeV]","E_{T,sublead tower} [GeV]",
    "Max Tower #chi^{2}","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in EMCal","#eta","#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}","E_{T,lead tower} [GeV]",
    "Max Tower #chi^{2}","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in EMCal","#eta","#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}",
    "Max Tower #chi^{2}","Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in EMCal","#eta","#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}",
    "Fraction of E_{T,lead jet} in OHCal","Fraction of E_{T,lead jet} in EMCal","#eta","#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}",
    "Fraction of E_{T,lead jet} in EMCal","#eta","#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}",
    "#eta","#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}",
    "#phi","E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}",
    "E_{T,lead jet}","#Delta#phi","E_{T,sublead jet}",
    "#Delta#phi","E_{T,sublead jet}",
    "#Delta#phi",
    "#Delta#phi"
  };

  TLegend* spectraLegend = new TLegend(0.57,0.55,0.9,0.65);
  spectraLegend->SetFillStyle(0);
  spectraLegend->SetFillColor(0);
  spectraLegend->SetTextFont(42);
  spectraLegend->SetBorderSize(0);
  spectraLegend->SetTextSize(0.02);
  TCanvas* canvas = new TCanvas("canvas", "Canvas", 1200, 1200);
  const int nSpectra = 32;
  const int nColors = 7;
  const int nLegend = 7;
  const int nRatio = 9;
  string spectraLegendNames[nLegend] = {"No Cuts","OHCal Strip Cut Only","Add #Delta#phi Cuts","Add E_{T} Cuts","Data All Cuts","Sim No Cuts","Sim All Cuts"};
  int spectraColors[nColors] = {kBlack,kRed+2,kOrange+2,kYellow+2,kGreen+2,kBlue-2,kViolet+2};
  string xSpectraTitle = "E_{T,jet} [GeV]";
  string ySpectraTitle = "#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]";
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.15);
  float stdsize = 0.02;
  float stdx = 0.6;
  float stdy = 0.85;
  int stdright = 0;
  const int ntext = 10;
  const int nSame = 3;
  string calcut = "";

  // 4/0, 4/3, 1/0, 2/0, 3/0
  TH1F* ratios[nRatio];
  gPad->SetLogy(0);
  for(int i=0; i<nRatio; ++i)
    {
      ratios[i] = new TH1F(("ratio_"+to_string(i)).c_str(),"",100,0,100);
    }
  for(int i=33; i<38; ++i)
    {
      jetSpectra.at(i)->Rebin(10);
    }
  //jetSpectra.at(41)->Rebin(10);
  jetSpectra.at(43)->Rebin(10);
  jetSpectra.at(44)->Rebin(10);
  jetSpectra.at(45)->Rebin(10);
  //jetSpectra.at(46)->Rebin(10);
  //jetSpectra.at(47)->Rebin(10);
  string texts[ntext] =
    {
      "Calorimeter Anti-k_{T} R=0.4",
      "|z_{vtx}| < 150 cm",
      "Jet-8 & MBDNS>=1 Triggered Data",
      "PYTHIA8 Jet-10 Sample Sim",
      "\\mathscr{L}_{\\text{data}}=13.37 \\text{pb}^{-1}"
    };
  TLegend* ratLeg = new TLegend(0.57,0.5,0.9,0.6);
  ratLeg->SetTextSize(0.02);
  jetSpectra.at(39)->Scale(1./jetSpectra.at(39)->Integral("WIDTH"));
  jetSpectra.at(47)->Scale(1./jetSpectra.at(47)->Integral("WIDTH"));
  jetSpectra.at(38)->Scale(1./jetSpectra.at(38)->Integral("WIDTH"));
  jetSpectra.at(46)->Scale(1./jetSpectra.at(46)->Integral("WIDTH"));

  ratios[0]->Divide(jetSpectra.at(34),jetSpectra.at(33));
  ratios[1]->Divide(jetSpectra.at(35),jetSpectra.at(33));
  ratios[2]->Divide(jetSpectra.at(36),jetSpectra.at(33));
  ratios[3]->Divide(jetSpectra.at(37),jetSpectra.at(33));
  ratios[4]->Divide(jetSpectra.at(37),jetSpectra.at(36));
  ratios[5]->Divide(jetSpectra.at(44),jetSpectra.at(43));
  ratios[6]->Divide(jetSpectra.at(44),jetSpectra.at(45));
  ratios[7]->Divide(jetSpectra.at(38),jetSpectra.at(46));
  ratios[8]->Divide(jetSpectra.at(39),jetSpectra.at(47));
  gPad->SetLogy();
  for(int i=0; i<48; ++i)
    {
      jetSpectra.at(i)->Draw();
      canvas->SaveAs(("output/chi2img/spectra"+to_string(i)+".pdf").c_str());
    }
  gPad->SetLogy(0);
  string ratLegNames[nRatio] = {"Data 2+ Jets / All Jets","2+ Jets & Hard #Delta#phi / All Jets","All cuts / All Jets","2+ Jets & All Cuts / All Jets","Data 2+ Jets & All Cuts / All Cuts","Sim 2+ Jets All Cuts / All Cuts","Sim 2+ Jets All Cuts / No Cuts","Data No Cuts / Sim No Cuts","Data All Cuts / Sim All Cuts"};
  canvas->SetRightMargin(0.05);
  for(int i=0; i<nRatio; ++i)
    {
      
      ratios[i]->GetYaxis()->SetTitleOffset(1.85);
      ratios[i]->SetMarkerColor(spectraColors[i%nColors]);
      ratios[i]->SetMarkerStyle(20);
      ratios[i]->SetMarkerSize(2);
      ratios[i]->SetLineColor(spectraColors[i%nColors]);
      ratios[i]->Scale(1./ratios[i]->GetBinWidth(1));
      ratios[i]->GetYaxis()->SetRangeUser(0,1.1);
      ratios[i]->GetYaxis()->SetTitle("Ratio");
      if(i!=0 && (i < 4 || i > 6)) continue;
      ratLeg->AddEntry(ratios[i],ratLegNames[i].c_str(),"p");
      ratios[i]->GetXaxis()->SetTitle("E_{T,lead} [GeV]");
      ratios[i]->Draw(i==0?"HIST P":"SAME HIST P");
    }
  stdy = 0.8;
    std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
    stdy = 0.85;
    ratLeg->SetFillStyle(0);
    ratLeg->SetFillColor(0);
    ratLeg->SetTextFont(42);
    ratLeg->SetBorderSize(0);
    ratLeg->SetTextSize(0.02);
    TLegend* ajratleg = new TLegend(0.57,0.5,0.9,0.6);
    ajratleg->SetFillStyle(0);
    ajratleg->SetFillColor(0);
    ajratleg->SetTextFont(42);
    ajratleg->SetBorderSize(0);
    ajratleg->SetTextSize(0.02);
    ratLeg->Draw();
    canvas->SaveAs("output/chi2img/evtRatios_h1.png");
    gPad->SetRightMargin(0.1);
    ratios[7]->GetYaxis()->SetRangeUser(0,2);
    ratios[7]->GetXaxis()->SetTitle("A_{J}");
    ratios[7]->GetYaxis()->SetTitle("Ratio");
    ajratleg->AddEntry(ratios[7],ratLegNames[7].c_str(),"p");
    ajratleg->AddEntry(ratios[8],ratLegNames[8].c_str(),"p");
    ratios[7]->Draw("PE");
    ratios[8]->Draw("SAME PE");
    ajratleg->Draw();
    texts[5] = "E_{T,lead} > 20 GeV";
    texts[6] = "E_{T,sub} > 10 GeV";
    std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
    canvas->SaveAs("output/chi2img/ajratio_h1.png");


    jetSpectra.at(39)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(39)->SetMarkerColor(spectraColors[0]);
    jetSpectra.at(39)->SetMarkerStyle(20);
    jetSpectra.at(39)->SetMarkerSize(2);
    jetSpectra.at(39)->SetLineColor(spectraColors[0]);
    jetSpectra.at(39)->Scale(1./jetSpectra.at(39)->GetBinWidth(1));
    //jetSpectra.at(39)->GetYaxis()->SetRangeUser(0,);
    jetSpectra.at(39)->GetXaxis()->SetTitle("A_{J}");
    jetSpectra.at(39)->GetYaxis()->SetTitle("Normalized Counts");
    jetSpectra.at(39)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(47)->SetMarkerColor(spectraColors[1]);
    jetSpectra.at(47)->SetMarkerStyle(20);
    jetSpectra.at(47)->SetMarkerSize(2);
    jetSpectra.at(47)->SetLineColor(spectraColors[1]);
    jetSpectra.at(47)->Scale(1./jetSpectra.at(47)->GetBinWidth(1));
    //jetSpectra.at(47)->GetYaxis()->SetRangeUser(1e-12,1e-4);
    jetSpectra.at(47)->GetXaxis()->SetTitle("A_{J}");
    jetSpectra.at(47)->GetYaxis()->SetTitle("Normalized Counts");
    TLegend* ajleg = new TLegend(0.57,0.5,0.9,0.6);
    ajleg->SetFillStyle(0);
    ajleg->SetFillColor(0);
    ajleg->SetTextFont(42);
    ajleg->SetBorderSize(0);
    ajleg->SetTextSize(0.02);
    jetSpectra.at(39)->Draw("PE");
    jetSpectra.at(47)->Draw("SAME PE");
    std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
    ajleg->AddEntry(jetSpectra.at(39),"Data All Cuts","p");
    ajleg->AddEntry(jetSpectra.at(47),"Sim All Cuts","p");
    ajleg->Draw();
    canvas->SaveAs("output/chi2img/aj_h1.png");


    jetSpectra.at(38)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(38)->SetMarkerColor(spectraColors[0]);
    jetSpectra.at(38)->SetMarkerStyle(20);
    jetSpectra.at(38)->SetMarkerSize(2);
    jetSpectra.at(38)->SetLineColor(spectraColors[0]);
    jetSpectra.at(38)->Scale(1./jetSpectra.at(38)->GetBinWidth(1));
    //jetSpectra.at(38)->GetYaxis()->SetRangeUser(0,);
    jetSpectra.at(38)->GetXaxis()->SetTitle("A_{J}");
    jetSpectra.at(38)->GetYaxis()->SetTitle("Normalized Counts");
    jetSpectra.at(38)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(46)->SetMarkerColor(spectraColors[1]);
    jetSpectra.at(46)->SetMarkerStyle(20);
    jetSpectra.at(46)->SetMarkerSize(2);
    jetSpectra.at(46)->SetLineColor(spectraColors[1]);
    jetSpectra.at(46)->Scale(1./jetSpectra.at(46)->GetBinWidth(1));
    //jetSpectra.at(46)->GetYaxis()->SetRangeUser(1e-12,1e-4);
    jetSpectra.at(46)->GetXaxis()->SetTitle("A_{J}");
    jetSpectra.at(46)->GetYaxis()->SetTitle("Normalized Counts");
    TLegend* ajleg2 = new TLegend(0.57,0.5,0.9,0.6);
    ajleg2->SetFillStyle(0);
    ajleg2->SetFillColor(0);
    ajleg2->SetTextFont(42);
    ajleg2->SetBorderSize(0);
    ajleg2->SetTextSize(0.02);
    jetSpectra.at(38)->Draw("PE");
    jetSpectra.at(46)->Draw("SAME PE");
    std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
    ajleg2->AddEntry(jetSpectra.at(38),"Data No Cuts","p");
    ajleg2->AddEntry(jetSpectra.at(46),"Sim No Cuts","p");
    ajleg2->Draw();
    canvas->SaveAs("output/chi2img/ajpre_h1.png");
    texts[5] = "";
    texts[6] = "";

    gPad->SetLogy();
    //jetSpectra.at(39)->Scale();
    //jetSpectra.at(40)->Scale(0.001);
    jetSpectra.at(0)->Rebin(10);
    
    jetSpectra.at(0)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(0)->SetMarkerColor(spectraColors[0]);
    jetSpectra.at(0)->SetMarkerStyle(20);
    jetSpectra.at(0)->SetMarkerSize(2);
    jetSpectra.at(0)->SetLineColor(spectraColors[0]);
    jetSpectra.at(0)->Scale(1./jetSpectra.at(0)->GetBinWidth(1));
    jetSpectra.at(0)->GetYaxis()->SetRangeUser(1e-12,1e-4);
    jetSpectra.at(0)->GetXaxis()->SetTitle(xSpectraTitle.c_str());
    jetSpectra.at(0)->GetYaxis()->SetTitle(ySpectraTitle.c_str());

  //jetSpectra.at(0)->GetXaxis()->SetTitleSize(0.02);
  //jetSpectra.at(0)->GetYaxis()->SetTitleSize(0.02);
  //jetSpectra.at(0)->GetXaxis()->SetLabelSize();
  //jetSpectra.at(0)->GetYaxis()->SetTitle();
  
  jetSpectra.at(0)->Draw("PE");
  cout << "drew 0" << endl;
  canvas->SaveAs("testimg.pdf");
  //sleep(30);
  TH1* touse[nLegend-1] = {jetSpectra.at(2),jetSpectra.at(16),jetSpectra.at(28),jetSpectra.at(31),jetSpectra.at(41),jetSpectra.at(42)};
  TH1F* prespcr = new TH1F("prespcr","",100,0,100);
  TH1F* specrat = new TH1F("specrat","",100,0,100);

  spectraLegend->AddEntry(jetSpectra.at(0),spectraLegendNames[0].c_str(),"p");
  jetSpectra.at(0)->Scale(1./280758718700.0);
  TLegend* specL2 = new TLegend(*spectraLegend);
  for(int i=1; i<nLegend; ++i)
    {
      if(i<nLegend-2) touse[i-1]->Scale(1./280758718700.0);
      touse[i-1]->Rebin(10);
      touse[i-1]->GetYaxis()->SetTitleOffset(1.85);
      touse[i-1]->SetMarkerColor(spectraColors[i]);
      touse[i-1]->SetMarkerStyle(20);
      touse[i-1]->SetMarkerSize(2);
      touse[i-1]->SetLineColor(spectraColors[i]);
      touse[i-1]->Scale(1./touse[i-1]->GetBinWidth(1));
      touse[i-1]->GetYaxis()->SetRangeUser(1e-12,1e-4);
      touse[i-1]->GetXaxis()->SetTitle(xSpectraTitle.c_str());
      touse[i-1]->GetYaxis()->SetTitle(ySpectraTitle.c_str());
      touse[i-1]->Draw("SAME PE");
      cout << "drew " << i << endl;
      if(i<4) spectraLegend->AddEntry(touse[i-1],spectraLegendNames[i].c_str(),"p");
      else if(i>3) specL2->AddEntry(touse[i-1],spectraLegendNames[i].c_str(),"p");
      if(i==4)
	{
	  spectraLegend->AddEntry(jetSpectra.at(31),"Add IHCal Frac Cut (All Cuts)","p");
	  spectraLegend->Draw();
	  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
	  canvas->SaveAs("output/chi2img/SpectrumCutShown_h1.png");
	  for(int j=0; j<i; ++j)
	    {
	      spectraLegend->RecursiveRemove(touse[j]);
	      cout << "removed " << j << endl;
	    }
	  canvas->Clear();
	  jetSpectra.at(0)->Draw("PE");
	}
      if(i==nLegend-1) touse[3]->Draw("SAME PE");
    }
  int bitset = 0;
  if(bitset) return 0;
  specL2->Draw();
  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
  canvas->SaveAs(("output/chi2img/jetSpectrumWithCuts-1_h1.png"));
  //for(int i=0; i<nSpectra; i+=3)
  //{
  //  if(i > nSpectra) break;

      /*
      if(i%3==0)
	{
	  calcut = "No Calo Frac Cut";
	}
      else if(i%3==1)
	{
	  calcut = "OH E Frac > 0.9";
	}
      else if(i%3==2)
	{
	  calcut = "EM E Frac > 0.9";
	}
      */
      /*
      jetSpectra.at(0)->Draw("PE");
      for(int j=1; j<nSame+1; ++j)
	{
	  if (i+j >= nSpectra) break;
	  jetSpectra.at(i+j)->Rebin(10);
	  jetSpectra.at(i+j)->GetYaxis()->SetRangeUser(0.1,1e6);
	  jetSpectra.at(i+j)->GetYaxis()->SetTitleOffset(1.85);
	  jetSpectra.at(i+j)->SetMarkerColor(spectraColors[j]);
	  jetSpectra.at(i+j)->SetMarkerStyle(20);
	  jetSpectra.at(i+j)->SetMarkerSize(2);
	  jetSpectra.at(i+j)->SetLineColor(spectraColors[j]);
	  jetSpectra.at(i+j)->Scale(1./jetSpectra.at(i+j)->GetBinWidth(1));
	  jetSpectra.at(i+j)->GetXaxis()->SetTitle(xSpectraTitle.c_str());
	  jetSpectra.at(i+j)->GetYaxis()->SetTitle(ySpectraTitle.c_str());
	  jetSpectra.at(i+j)->Draw("SAME PE");
	  
	}
      */
      //std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
      //canvas->SaveAs(("output/chi2img/jetSpectrumWithCuts"+to_string(i)+".png").c_str());
      //if(i<6)
      //{
      //if(i==0) jetSpectra.at(i)->Draw("PE");
      //else jetSpectra.at(i)->Draw("SAME PE");
      //}
      //if(i<6) spectraLegend->AddEntry(jetSpectra.at(i),spectraLegendNames[i].c_str(),"p");
  //}


  gPad->SetLogy(0);
  specrat->Divide(jetSpectra.at(31),jetSpectra.at(42));

  specrat->GetYaxis()->SetTitleOffset(2.0);
  specrat->SetMarkerColor(spectraColors[0]);
  specrat->SetMarkerStyle(20);
  specrat->SetMarkerSize(2);
  specrat->SetLineColor(spectraColors[0]);
  specrat->Scale(1./specrat->GetBinWidth(1));
  specrat->GetYaxis()->SetRangeUser(0,1);
  specrat->GetXaxis()->SetTitle("E_{T} [GeV]");
  specrat->GetYaxis()->SetTitle("Data/Sim All Cuts Spectrum Ratio");
  specrat->Draw("PE");
  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);

  canvas->SaveAs("output/chi2img/datsimrat_h1.png");


  prespcr->Divide(jetSpectra.at(0),jetSpectra.at(41));
  prespcr->GetYaxis()->SetTitleOffset(2.0);
  prespcr->SetMarkerColor(spectraColors[1]);
  prespcr->SetMarkerStyle(20);
  prespcr->SetMarkerSize(2);

  prespcr->SetLineColor(spectraColors[1]);
  prespcr->Scale(1./prespcr->GetBinWidth(1));
  prespcr->GetXaxis()->SetTitle("E_{T} [GeV]");
  prespcr->GetYaxis()->SetTitle("Data/Sim No Cuts Spectrum Ratio");
  prespcr->GetYaxis()->SetRangeUser(0,3);
  prespcr->Draw("PE");
  //spectraLegend->Draw();
  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);

  canvas->SaveAs("output/chi2img/precutdatsimrat_h1.png");
  gPad->SetBottomMargin(0.1);
  gPad->SetRightMargin(0.1);

  TH1F* xJRat[2];
  xJRat[0] = new TH1F("xjrat0","",100,0,1);
  xJRat[1] = new TH1F("xjrat1","",100,0,1);
  TLegend* xJLeg = new TLegend(0.13,0.4,0.46,0.6);
  int xJColors[4] = {kBlack,kRed+2,kBlue+2,kGreen+2};
  string xJLegNames[4] = {"Data No Cuts","Data With Cuts","Sim No Cuts","Sim With Cuts"};
  for(int i=0; i<4; ++i)
    {
      xJ.at(i)->Scale(1./xJ.at(i)->Integral("WIDTH"));
      xJ.at(i)->SetLineColor(xJColors[i]);
      xJ.at(i)->SetMarkerColor(xJColors[i]);
      xJ.at(i)->SetMarkerStyle(20);
      xJ.at(i)->SetMarkerSize(2);
      xJ.at(i)->GetYaxis()->SetTitleOffset(2.0);
      xJ.at(i)->GetYaxis()->SetTitle("Normalized Counts");
      xJ.at(i)->GetXaxis()->SetTitle("X_{J}");
      if(i==0) xJ.at(i)->Draw("PE");
      else xJ.at(i)->Draw("SAME PE");
      xJLeg->AddEntry(xJ.at(i),xJLegNames[i].c_str(),"p");
    }
  texts[5] = "E_{T,lead} > 20 GeV";
  texts[6] = "E_{T,sub} > 10 GeV";
  stdx = 0.13;
  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
  xJLeg->Draw();
  canvas->SaveAs("output/chi2img/xJ_h1.png");


  int histgroup;
  for(int i=0; i<histograms.size(); ++i)
    {
      if(i>3*nhist-1)
	{
	  histgroup = 3;
	}
      else if(i>2*nhist-1)
	{
	  histgroup = 2;
	}
      else if(i>nhist-1)
	{
	  histgroup = 1;
	}
      else
	{
	  histgroup = 0;
	}
      cout << histograms.at(i)->GetName() << endl;
      //histograms.at(i)->Scale(1/histograms.at(i)->Integral("WIDTH"));
      //cout << "scaled hist" << endl;
      string outdrawname = "output/chi2img/"+names[i%nhist] + "_" + to_string(i/nhist);
      const int ntype = 24;
      string texts[2];
      int threshes[ntype/2] = {8,15,20,25,35,40,45,50,55,60,65,70};
      int slt[ntype/2] = {8,8,10,10,15,15,20,20,25,25,30,30};
      //texts[0] = "E_{T,lead jet} > "+to_string(threshes[(i/nhist)%(ntype/2)])+" GeV";
	if(i<(ntype/2)*nhist) texts[0] = "Events contain only one jet with E_{T,jet} > "+to_string(threshes[(i/nhist)%(ntype/2)])+" GeV";
	else
	{
	  texts[0] = "Events contain at least one jet with E_{T,jet} > "+to_string(threshes[(i/nhist)%(ntype/2)])+" GeV";
	}

	if(i<(ntype/2)*nhist) texts[1] = "Cuts Applied";//"And no subleading jet with E_{T,jet} > "+to_string(slt[(i/nhist)%(ntype/2)])+" GeV";
	else
	{
	  texts[1] = "No Cuts Applied";// "And at least one subleading jet with E_{T,jet} > "+to_string(slt[(i/nhist)%(ntype/2)])+" GeV";
	}
       
      //if(i<(ntype/2)*nhist) histograms.at(i)->Add(histograms.at(i),histograms.at(i+180));
      cout << "prep to draw" << endl;
      //histograms.at(i)->Rebin2D(4,4);
      FormatAndDrawHistogram(
			     canvas,histograms.at(i),
			     outdrawname,xtitles[i%nhist].c_str(),ytitles[i%nhist].c_str(),"N_{jet}",
			     0.15,0.2,0.15,0.15,
			     1.85,1.85,2,
			     0.03,0.03,0.03,
			     0.03,0.03,0.03,
			     texts
			     );
    }
  // Clean up
  file->Close();
  return 0;
}
