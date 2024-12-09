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
  gPad->SetLogz(0);
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


void draw_chi2hists(const TString& fileName) {
  // Open the ROOT file
  gStyle->SetOptTitle(0);
  const int nhist = 92;
  TFile* file = TFile::Open(fileName);
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << fileName << std::endl;
    return;
  }

  // Vector to hold pointers to TH2 histograms
  std::vector<TH2*> histograms;
  std::vector<TH1*> jetSpectra;
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
    "eta_phi", "eta_jet_ET", "eta_dphi", "eta_subjet_ET",
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

  TLegend* spectraLegend = new TLegend(0.45,0.5,0.85,0.7);
  spectraLegend->SetFillStyle(0);
  spectraLegend->SetFillColor(0);
  spectraLegend->SetTextFont(42);
  spectraLegend->SetBorderSize(0);

  TCanvas* canvas = new TCanvas("canvas", "Canvas", 1200, 1200);
  const int nSpectra = 32;
  const int nColors = 4;
  string spectraLegendNames[nSpectra] = {"No Cuts","No Cuts OH E_{T} Frac > 0.9","No Cuts EM E_{T} Frac > 0.9","OH 90% #phi","OH 90% #phi OH E_{T} Frac > 0.9","OH 90% #phi OH E_{T} Frac > 0.9"};
  int spectraColors[nColors] = {kBlack,kRed,kGreen,kBlue};
  string xSpectraTitle = "E_{T,jet} [GeV]";
  string ySpectraTitle = "dN_{jet}/dE_{T,jet}";
  gPad->SetLeftMargin(0.15);
  gPad->SetLogy();
  float stdsize = 0.025;
  float stdx = 0.6;
  float stdy = 0.85;
  int stdright = 0;
  const int ntext = 2;
  const int nSame = 3;
  string calcut = "";
  string texts[ntext] =
    {
      "Calorimeter Anti-k_{T} R=0.4",
      "|z_{vtx}| < 150 cm"
    };
  jetSpectra.at(0)->Rebin(10);
  jetSpectra.at(0)->GetYaxis()->SetTitleOffset(1.85);
  jetSpectra.at(0)->SetMarkerColor(spectraColors[0]);
  jetSpectra.at(0)->SetMarkerStyle(20);
  jetSpectra.at(0)->SetMarkerSize(2);
  jetSpectra.at(0)->SetLineColor(spectraColors[0]);
  jetSpectra.at(0)->Scale(1./jetSpectra.at(0)->GetBinWidth(1));
  jetSpectra.at(0)->GetYaxis()->SetRangeUser(0.5,5e6);
  jetSpectra.at(0)->GetXaxis()->SetTitle(xSpectraTitle.c_str());
  jetSpectra.at(0)->GetYaxis()->SetTitle(ySpectraTitle.c_str());
  jetSpectra.at(0)->Draw("PE");
  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
  canvas->SaveAs(("output/chi2img/jetSpectrumWithCuts-1.png"));
  for(int i=0; i<nSpectra; i+=3)
    {
      if(i > nSpectra) break;

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
      std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
      canvas->SaveAs(("output/chi2img/jetSpectrumWithCuts"+to_string(i)+".png").c_str());
      //if(i<6)
      //{
      //if(i==0) jetSpectra.at(i)->Draw("PE");
      //else jetSpectra.at(i)->Draw("SAME PE");
      //}
      //if(i<6) spectraLegend->AddEntry(jetSpectra.at(i),spectraLegendNames[i].c_str(),"p");
    }

  
  //spectraLegend->Draw();
  
  gPad->SetRightMargin(0.1);

  gPad->SetLogy(0);
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

	if(i<(ntype/2)*nhist) texts[1] = "And no subleading jet with E_{T,jet} > "+to_string(slt[(i/nhist)%(ntype/2)])+" GeV";
	else
	{
	  texts[1] = "And at least one subleading jet with E_{T,jet} > "+to_string(slt[(i/nhist)%(ntype/2)])+" GeV";
	}
	

      //if(i<(ntype/2)*nhist) histograms.at(i)->Add(histograms.at(i),histograms.at(i+180));
      cout << "prep to draw" << endl;
      histograms.at(i)->Rebin2D(4,4);
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
