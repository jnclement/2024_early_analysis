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
//#include <RooUnfold.h>
//#include <RooUnfoldResponse.h>
//#include <RooUnfoldBayes.h>
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
  std_text(canvas,texts,2,0.025,0.15,0.98,0);
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


void draw_chi2hists(const std::string fileName, const std::string simFileName) {
  // Open the ROOT file
  bool isJet30 = simFileName.find("jet30") != std::string::npos?true:false;
  bool isPythia = simFileName.find("herwig") != std::string::npos?false:true;
  float Lintpb = 16.8155;//9.2*0.45;
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  const int nhist = 92;
  TFile* file = TFile::Open(fileName.c_str());
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file: " << fileName << std::endl;
    return;
  }

  TFile* simfile = TFile::Open(simFileName.c_str());
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

  TH2D* hResp;
  std::vector<TH2*> hists2_sim;
  std::vector<TH2*> corPlotsLJet;
  TList* simkeys = simfile->GetListOfKeys();
  cout << "nKey: " << keys->GetSize() << endl;
  TIter simiter(simkeys);
  TKey* simkey;
  bool haveresp = false;
  while ((simkey = (TKey*)simiter())) {
    cout << "key: "<< simkey << endl;
    TObject* obj = simkey->ReadObj();
    cout << "obj before: "<<obj << endl;
    if (obj && (obj->IsA()->InheritsFrom("TH2") || obj->IsA()->InheritsFrom("TH1")) ) {
      TString histName = obj->GetName();
      cout << histName << endl;
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
      if(histName.BeginsWith("hists2"))
	{
	  hists2_sim.push_back((TH2*)obj);
	}
      if(histName.BeginsWith("lJet"))
	{
	  corPlotsLJet.push_back((TH2*)obj);
	}
      if(histName.BeginsWith("response") && !haveresp)
	{
	  cout << "hResp: " << obj << endl;
	  hResp = (TH2D*)obj;
	  haveresp = true;
	}
    }
    else cout <<"obj again: "<< obj << endl;
  }

  TH1F* dijetCheckRat = (TH1F*)file->Get("dijetCheckRat");
  TH1F* dijetCheckRatFull = (TH1F*)file->Get("dijetCheckRatFull");

  TH1F* dijetCheckRatSim = (TH1F*)simfile->Get("dijetCheckRatSim");
  TH1F* dijetCheckRatFullSim = (TH1F*)simfile->Get("dijetCheckRatFullSim");
  //const int nbin = 9;
  //float bins[nbin+1] =
  const int nbin = 12;
  float bins[nbin+1] = {12, 15, 18, 21, 24, 28, 32, 36, 40, 45, 50, 55, 60};//{12,16,21,26,32,38,45,52,60,70};

  hResp->Draw("COLZ");
  // Output the names of the retrieved histograms
  std::cout << "Retrieved 2D histograms:" << std::endl;
  /*
  cout << histograms.size() << endl;
  for (auto* hist : histograms) {
    cout << hist << endl;
    std::cout << hist->GetName() << std::endl;
  }
  */
  cout <<"hresp: " << hResp << endl;
  //RooUnfoldResponse response(jetSpectra.at(41),jetSpectra.at(52),hResp);

  //RooUnfoldBayes unfold(&response, jetSpectra.at(31), 4);

  //TH1D* hUnfold = (TH1D*) unfold.Hunfold(RooUnfold::kErrors);
  //cout <<"HUNFOLD PARAMS:" << endl << hUnfold->GetNbinsX() << endl <<hUnfold->GetXaxis()->GetBinWidth(1) << endl;
  /*
  for(int i=1; i<11; ++i)
    {
      cout << hUnfold->GetBinContent(i) << endl;
    }
  */

  TH1F* xJ2 = static_cast<TH1F*>(file->Get("new_xJ2"));
  TH1F* xJ3 = static_cast<TH1F*>(file->Get("new_xJ3"));
  
  TH1F* specspec = static_cast<TH1F*>(simfile->Get("specspec"));

  
  TH2F* asdich_all = (TH2F*)file->Get("asdich_all");
  TH2F* asdich_fail = (TH2F*)file->Get("asdich_fail");
  TH2F* asdich_all_sim = (TH2F*)simfile->Get("asdich_all_sim");
  TH2F* asdich_fail_sim = (TH2F*)simfile->Get("asdich_fail_sim");
  TH2F* asdich_dphi = (TH2F*)file->Get("asdich_dphi");
  TH2F* asdich_dphi_sim = (TH2F*)simfile->Get("asdich_dphi_sim");
  

  if(isJet30)
    {
      xJ[0] = xJ2;
      xJ[1] = xJ3;
    }

  string xtcp[6] = {"E_{T,EM}/E_{T,jet}","E_{T,EM}/E_{T,jet}","E_{T,EM}/E_{T,jet}","E_{T,jet}","E_{T,jet}","#eta"};

  string ytcp[6] = {"E_{T,jet}","#eta","#phi","#eta","#phi","#phi"};

  string xtitles[nhist] = {
    "#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower","#chi^{2} of Maximum Energy Tower",
    "N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet","N Bad #chi^{2} Towers in Lead Jet",
    "E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]","E_{T,lead tower} - E_{T,sublead tower} [GeV]",
    "E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]","E_{T,sublead tower} [GeV]",
    "E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]","E_{T,lead tower} [GeV]",
    "#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}","#epsilon_{lead jet}",
    "Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}","Max Tower #chi^{2}",
    "Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in OHCal",
    "Fraction of E_{T,jet} in EMCal","Fraction of E_{T,jet} in EMCal","Fraction of E_{T,jet} in EMCal","Fraction of E_{T,jet} in EMCal","Fraction of E_{T,jet} in EMCal",
    "#eta","#eta","#eta","#eta",
    "#phi","#phi","#phi",
    "E_{T,jet}","E_{T,jet}",
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
    "N Bad #chi^{2}","E_{T,lead tower} - E_{T,sublead tower} [GeV]","Max Tower #chi^{2}","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in EMCal","#eta","#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}","E_{T,lead tower} [GeV]","E_{T,sublead tower} [GeV]",
    "E_{T,lead tower} - E_{T,sublead tower} [GeV]","Max Tower #chi^{2}","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in EMCal","#eta","#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}","E_{T,lead tower} [GeV]","E_{T,sublead tower} [GeV]",
    "Max Tower #chi^{2}","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in EMCal","#eta","#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}","E_{T,lead tower} [GeV]","E_{T,sublead tower} [GeV]",
    "Max Tower #chi^{2}","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in EMCal","#eta","#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}","E_{T,lead tower} [GeV]",
    "Max Tower #chi^{2}","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in EMCal","#eta","#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}","#epsilon_{lead jet}",
    "Max Tower #chi^{2}","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in EMCal","#eta","#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}",
    "Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in EMCal","#eta","#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}",
    "Fraction of E_{T,jet} in EMCal","#eta","#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}",
    "#eta","#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}",
    "#phi","E_{T,jet}","#Delta#phi","E_{T,sublead jet}",
    "E_{T,jet}","#Delta#phi","E_{T,sublead jet}",
    "#Delta#phi","E_{T,sublead jet}",
    "#Delta#phi",
    "#Delta#phi"
  };

  TLegend* spectraLegend = new TLegend(0.55,0.55,0.9,0.65);
  spectraLegend->SetFillStyle(0);
  spectraLegend->SetFillColor(0);
  spectraLegend->SetTextFont(42);
  spectraLegend->SetBorderSize(0);
  spectraLegend->SetTextSize(0.02);
  TCanvas* canvas = new TCanvas("canvas", "Canvas", 1200, 1200);

  string adt[2] = {"Data; all entries are dijets, no back to back requirement","All leading jets passing all standard cuts included"};
  gPad->SetLogz();
  FormatAndDrawHistogram(canvas,asdich_all,"output/chi2img/asdich_all","A_{J}","E_{T,lead} [GeV]","Counts",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 adt);

  adt[0] = "Data; all entries are dijets with #Delta#phi > 3#pi/4";
  adt[1] = "All events pass all standard cuts, leading jets only";
  FormatAndDrawHistogram(canvas,asdich_dphi,"output/chi2img/asdich_dphi","A_{J}","E_{T,lead} [GeV]","Counts",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 adt);
  adt[0] = "Data; all entries are dijets with #Delta#phi > 3#pi/4, leading jets only";
  adt[1] = "All events pass standard cuts but FAIL when back-to-back dijets not preserved";
  FormatAndDrawHistogram(canvas,asdich_fail,"output/chi2img/asdich_fail","A_{J}","E_{T,lead} [GeV]","Counts",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 adt);


  string adts[2] = {"Sim; all entries are dijets, no back to back requirement","All leading jets passing all standard cuts included"};
  gPad->SetLogz();
  FormatAndDrawHistogram(canvas,asdich_all_sim,"output/chi2img/asdich_all_sim","A_{J}","E_{T,lead} [GeV]","Counts",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 adts);

  adts[0] = "Sim; all entries are dijets with #Delta#phi > 3#pi/4";
  adts[1] = "All events pass all standard cuts, leading jets only";
  FormatAndDrawHistogram(canvas,asdich_dphi_sim,"output/chi2img/asdich_dphi_sim","A_{J}","E_{T,lead} [GeV]","Counts",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 adts);
  adts[0] = "Sim; all entries are dijets with #Delta#phi > 3#pi/4, leading jets only";
  adts[1] = "All events pass standard cuts but FAIL when back-to-back dijets not preserved";
  FormatAndDrawHistogram(canvas,asdich_fail_sim,"output/chi2img/asdich_fail_sim","A_{J}","E_{T,lead} [GeV]","Counts",
			 0.15,0.2,0.15,0.15,
			 1.85,1.85,2,
			 0.03,0.03,0.03,
			 0.03,0.03,0.03,
			 adts);

  canvas->SetTopMargin(0.1);
  
  const int nSpectra = 32;
  const int nColors = 7;
  const int nLegend = 7;
  const int nRatio = 11;
  string spectraLegendNames[nLegend] = {"Data No Cuts","OHCal Strip Cut Only","Add #Delta#phi Cuts","Add E_{T} Cuts","Data All Cuts","Sim No Cuts","Sim All Cuts"};
  int spectraColors[nColors] = {kBlack,kOrange+2,kBlue-2,kYellow+2,kGreen+2,kRed+2,kViolet+2};
  string xSpectraTitle = "E_{T,jet} [GeV]";
  string ySpectraTitle = "#frac{d^{2}#sigma}{2#pi dE_{T}d#eta} [mb/GeV]";//"#frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]";
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
  int rebinfactor = 40;
  int ratioBins[nRatio] = {1000/rebinfactor,1000/rebinfactor,1000/rebinfactor,1000/rebinfactor,1000/rebinfactor,1000/rebinfactor,1000/rebinfactor,25,25,1000/rebinfactor,1000/rebinfactor};
  for(int i=0; i<nRatio; ++i)
    {
      ratios[i] = new TH1F(("ratio_"+to_string(i)).c_str(),"",ratioBins[i],0,(i<7||i>8)?100:1);
    }
  /*
  for(int i=33; i<38; ++i)
    {
      jetSpectra.at(i)->Rebin(10);
    }
  //jetSpectra.at(41)->Rebin(10);
  jetSpectra.at(43)->Rebin(10);
  jetSpectra.at(44)->Rebin(10);
  jetSpectra.at(45)->Rebin(10);
  jetSpectra.at(46)->Rebin(10);
  jetSpectra.at(47)->Rebin(10);
  */
  string texts[ntext] =
    {
      "Calorimeter Anti-k_{T} R=0.4",
      "|z_{vtx}| < 30 cm",
      "Jet-8 & MBDNS>=1 Triggered Data",
      std::string(isPythia?"PYTHIA":"HERWIG")+" "+(isJet30?"Jet-30":"Jet-10")+" Sample Sim",
      "\\mathscr{L}_{\\text{data}}=16.82 \\text{pb}^{-1}"
    };
  TLegend* ratLeg = new TLegend(0.57,0.4,0.9,0.6,"Dijets / Event:");
  ratLeg->SetTextSize(0.02);
  /*
  jetSpectra.at(38)->Rebin(4);
  jetSpectra.at(46)->Rebin(4);
  jetSpectra.at(39)->Rebin(4);
  jetSpectra.at(47)->Rebin(4);
  */

  TF1* fturn = new TF1("fturn","[0]*(TMath::Erf((x-[1])/[2]) + 1)", 4, 30);
  fturn->SetParameter(0,0.470647);
  fturn->SetParameter(1,11.6307);
  fturn->SetParameter(2,4.36329);
  for(int i=1; i<5; ++i)
    {
      jetSpectra.at(0)->SetBinContent(i,jetSpectra.at(0)->GetBinContent(i)/fturn->Eval(jetSpectra.at(0)->GetBinCenter(i)));
      jetSpectra.at(31)->SetBinContent(i,jetSpectra.at(31)->GetBinContent(i)/fturn->Eval(jetSpectra.at(31)->GetBinCenter(i)));
      jetSpectra.at(30)->SetBinContent(i,jetSpectra.at(30)->GetBinContent(i)/fturn->Eval(jetSpectra.at(30)->GetBinCenter(i)));

      jetSpectra.at(0)->SetBinError(i,jetSpectra.at(0)->GetBinError(i)/fturn->Eval(jetSpectra.at(0)->GetBinCenter(i)));
      jetSpectra.at(31)->SetBinError(i,jetSpectra.at(31)->GetBinError(i)/fturn->Eval(jetSpectra.at(31)->GetBinCenter(i)));
      jetSpectra.at(30)->SetBinError(i,jetSpectra.at(30)->GetBinError(i)/fturn->Eval(jetSpectra.at(30)->GetBinCenter(i)));

    }
  jetSpectra.at(46)->Rebin(40);
  jetSpectra.at(47)->Rebin(40);
  jetSpectra.at(38)->Rebin(40);
  jetSpectra.at(39)->Rebin(40);
  jetSpectra.at(39)->Scale(1./jetSpectra.at(39)->Integral());
  jetSpectra.at(47)->Scale(1./jetSpectra.at(47)->Integral());
  jetSpectra.at(38)->Scale(1./jetSpectra.at(38)->Integral());
  jetSpectra.at(46)->Scale(1./jetSpectra.at(46)->Integral());
  cout << jetSpectra.at(33)->GetNbinsX() << endl;
  cout << jetSpectra.at(34)->GetNbinsX() << endl;
  cout << jetSpectra.at(35)->GetNbinsX() << endl;
  cout << jetSpectra.at(36)->GetNbinsX() << endl;
  cout << jetSpectra.at(37)->GetNbinsX() << endl;
  cout << jetSpectra.at(38)->GetNbinsX() << endl;
  cout << jetSpectra.at(39)->GetNbinsX() << endl;
  cout << jetSpectra.at(40)->GetNbinsX() << endl;
  cout << jetSpectra.at(41)->GetNbinsX() << endl;
  cout << jetSpectra.at(42)->GetNbinsX() << endl;
  cout << jetSpectra.at(43)->GetNbinsX() << endl;
  cout << jetSpectra.at(44)->GetNbinsX() << endl;
  cout << jetSpectra.at(45)->GetNbinsX() << endl;
  cout << jetSpectra.at(46)->GetNbinsX() << endl;
  cout << jetSpectra.at(47)->GetNbinsX() << endl;
  cout << dijetCheckRat->GetNbinsX() << endl;
  cout << dijetCheckRatFull->GetNbinsX() << endl;
  cout << dijetCheckRatSim->GetNbinsX() << endl;
  cout << dijetCheckRatFullSim->GetNbinsX() << endl;
  jetSpectra.at(33)->Rebin(40);
  jetSpectra.at(34)->Rebin(40);
  jetSpectra.at(35)->Rebin(40);
  jetSpectra.at(36)->Rebin(40);
  jetSpectra.at(37)->Rebin(40);
  jetSpectra.at(43)->Rebin(40);
  jetSpectra.at(44)->Rebin(40);
  jetSpectra.at(45)->Rebin(40);
  //jetSpectra.at(46)->Rebin(40);
  //jetSpectra.at(47)->Rebin(40);
  dijetCheckRat->Rebin(40);
  dijetCheckRatFull->Rebin(40);
  cout << "illegal?" << endl;
  dijetCheckRatSim->Rebin(40);
  dijetCheckRatFullSim->Rebin(40);
  ratios[0]->Divide(jetSpectra.at(34),jetSpectra.at(33),1,1,"B");
  cout << "1" << endl;
  ratios[1]->Divide(jetSpectra.at(35),jetSpectra.at(33),1,1,"B");
  cout << "2" << endl;
  ratios[2]->Divide(jetSpectra.at(36),jetSpectra.at(33),1,1,"B");
  cout << "3" << endl;
  ratios[3]->Divide(jetSpectra.at(37),jetSpectra.at(33),1,1,"B");
  cout << "4" << endl;
  ratios[4]->Divide(jetSpectra.at(37),jetSpectra.at(36),1,1,"B");
  cout << "5" << endl;
  ratios[5]->Divide(jetSpectra.at(44),jetSpectra.at(43),1,1,"B");
  cout << "6" << endl;
  ratios[6]->Divide(jetSpectra.at(44),jetSpectra.at(45),1,1,"B");
  cout << "7" << endl;
  ratios[7]->Divide(jetSpectra.at(38),jetSpectra.at(46),1,1,"B");
  cout << "8" << endl;
  ratios[8]->Divide(jetSpectra.at(39),jetSpectra.at(47),1,1,"B");
  cout << "9" << endl;
  ratios[9]->Divide(dijetCheckRat,dijetCheckRatFull,1,1,"B");
  ratios[10]->Divide(dijetCheckRatSim,dijetCheckRatFullSim,1,1,"B");
  gPad->SetLogy();
  gStyle->SetOptTitle();
  for(int i=0; i<54; ++i)
    {
      jetSpectra.at(i)->Draw();
      canvas->SaveAs(("output/chi2img/spectra"+to_string(i)+".pdf").c_str());
    }
  gStyle->SetOptTitle(0);
  gPad->SetLogy(0);
  string ratLegNames[nRatio] = {"Data No Cuts","2+ Jets & Hard #Delta#phi / All Jets","All cuts / All Jets","2+ Jets & All Cuts / All Jets","Data All Cuts","Reco Sim All Cuts","Reco Sim No Cuts","Data No Cuts / Reco Sim No Cuts","Data All Cuts / Reco Sim All Cuts"};
  canvas->SetRightMargin(0.05);
  for(int i=0; i<nRatio; ++i)
    {
      
      ratios[i]->GetYaxis()->SetTitleOffset(1.85);
      ratios[i]->SetMarkerColor(spectraColors[i%nColors]);
      ratios[i]->SetMarkerStyle(20);
      ratios[i]->SetMarkerSize(2);
      ratios[i]->SetLineColor(spectraColors[i%nColors]);
      //ratios[i]->Scale(1./ratios[i]->GetBinWidth(1));
      ratios[i]->GetYaxis()->SetRangeUser(0,1.1);
      ratios[i]->GetYaxis()->SetTitle("Ratio");
      if(i!=0 && (i < 4 || i > 6)) continue;//kif(i!=4 && i!= 5) continue;//
      ratLeg->AddEntry(ratios[i],ratLegNames[i].c_str(),"p");
      ratios[i]->GetXaxis()->SetTitle("E_{T,lead} [GeV]");
      ratios[i]->Draw(i==0?"PE":"SAME PE");
    }
  ratios[9]->SetLineColor(kBlue);
  ratios[9]->SetMarkerColor(kBlue);
  ratios[9]->SetMarkerStyle(20);
  //ratios[9]->Draw("SAME PE");

  ratios[10]->SetLineColor(kOrange);
  ratios[10]->SetMarkerColor(kOrange);
  ratios[10]->SetMarkerStyle(20);
  //ratios[10]->Draw("SAME PE");
  stdy = 0.8;
    std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
    stdy = 0.85;
    ratLeg->SetFillStyle(0);
    ratLeg->SetFillColor(0);
    ratLeg->SetTextFont(42);
    ratLeg->SetBorderSize(0);
    ratLeg->SetTextSize(0.02);
    //ratLeg->AddEntry(ratios[10],"Reco Sim All Cuts (No Dijet Check)","p");
    //ratLeg->AddEntry(ratios[9],"Data All Cuts (No Dijet Check)","p");
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
    ratios[7]->GetXaxis()->SetRangeUser(0,1.1);
    ratios[7]->Draw("PE");
    ratios[8]->Draw("SAME PE");
    ajratleg->Draw();
    texts[5] = "E_{T,lead} > 20 GeV, #Delta#phi > 3#pi/4";
    texts[6] = "E_{T,sub} > 10 GeV";
    std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
    canvas->SaveAs("output/chi2img/ajratio_h1.png");


    jetSpectra.at(39)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(39)->SetMarkerColor(spectraColors[0]);
    jetSpectra.at(39)->SetMarkerStyle(20);
    jetSpectra.at(39)->SetMarkerSize(2);
    jetSpectra.at(39)->SetLineColor(spectraColors[0]);
    jetSpectra.at(39)->Scale(1./jetSpectra.at(39)->GetBinWidth(1));
    jetSpectra.at(39)->GetYaxis()->SetRangeUser(0,5);
    jetSpectra.at(39)->GetXaxis()->SetTitle("A_{J}");
    jetSpectra.at(39)->GetYaxis()->SetTitle("Normalized Counts");
    jetSpectra.at(39)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(47)->SetMarkerColor(spectraColors[1]);
    jetSpectra.at(47)->SetMarkerStyle(20);
    jetSpectra.at(47)->SetMarkerSize(2);
    jetSpectra.at(47)->SetLineColor(spectraColors[1]);
    jetSpectra.at(47)->Scale(1./jetSpectra.at(47)->GetBinWidth(1));
    //jetSpectra.at(47)->GetYaxis()->SetRangeUser(1e-13,1e-4);
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
    ajleg->AddEntry(jetSpectra.at(47),"Reco Sim All Cuts","p");
    ajleg->Draw();
    canvas->SaveAs("output/chi2img/aj_h1.png");


    jetSpectra.at(38)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(38)->SetMarkerColor(spectraColors[0]);
    jetSpectra.at(38)->SetMarkerStyle(20);
    jetSpectra.at(38)->SetMarkerSize(2);
    jetSpectra.at(38)->SetLineColor(spectraColors[0]);
    jetSpectra.at(38)->Scale(1./jetSpectra.at(38)->GetBinWidth(1));
    jetSpectra.at(38)->GetYaxis()->SetRangeUser(0,5);
    jetSpectra.at(38)->GetXaxis()->SetTitle("A_{J}");
    jetSpectra.at(38)->GetYaxis()->SetTitle("Normalized Counts");
    jetSpectra.at(38)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(46)->SetMarkerColor(spectraColors[1]);
    jetSpectra.at(46)->SetMarkerStyle(20);
    jetSpectra.at(46)->SetMarkerSize(2);
    jetSpectra.at(46)->SetLineColor(spectraColors[1]);
    jetSpectra.at(46)->Scale(1./jetSpectra.at(46)->GetBinWidth(1));
    //jetSpectra.at(46)->GetYaxis()->SetRangeUser(1e-13,1e-4);
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
    ajleg2->AddEntry(jetSpectra.at(46),"Reco Sim No Cuts","p");
    ajleg2->Draw();
    canvas->SaveAs("output/chi2img/ajpre_h1.png");
    texts[5] = "";
    texts[6] = "";

    gPad->SetLogy();
    //jetSpectra.at(39)->Scale();
    //jetSpectra.at(40)->Scale(0.001);
    //jetSpectra.at(0)->Rebin(10);
    
    jetSpectra.at(0)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(0)->SetMarkerColor(spectraColors[0]);
    jetSpectra.at(0)->SetMarkerStyle(21);
    jetSpectra.at(0)->SetMarkerSize(2);
    jetSpectra.at(0)->SetLineColor(spectraColors[0]);
    for(int i=1; i<jetSpectra.at(0)->GetNbinsX()+1; ++i)
      {
	jetSpectra.at(0)->SetBinContent(i,jetSpectra.at(0)->GetBinContent(i)/jetSpectra.at(0)->GetBinWidth(i));
	jetSpectra.at(0)->SetBinError(i,jetSpectra.at(0)->GetBinError(i)/jetSpectra.at(0)->GetBinWidth(i));
      }
    jetSpectra.at(0)->Scale(1./(Lintpb*1e9)/1.4/(2*M_PI));//0.45*1.79769e11));
    jetSpectra.at(0)->GetYaxis()->SetRangeUser(1e-13,1e-4);
    jetSpectra.at(0)->GetXaxis()->SetTitle(xSpectraTitle.c_str());
    jetSpectra.at(0)->GetYaxis()->SetTitle(ySpectraTitle.c_str());

    //jetSpectra.at(30)->Rebin(10);
    jetSpectra.at(30)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(30)->SetMarkerColor(kBlue);
    jetSpectra.at(30)->SetMarkerStyle(20);
    jetSpectra.at(30)->SetMarkerSize(2);
    jetSpectra.at(30)->SetLineColor(kBlue);
    for(int i=1; i<jetSpectra.at(30)->GetNbinsX()+1; ++i)
      {
	cout << endl << endl << endl << jetSpectra.at(31)->GetBinContent(i) << endl << endl << endl;	
	jetSpectra.at(30)->SetBinContent(i,jetSpectra.at(30)->GetBinContent(i)/jetSpectra.at(30)->GetBinWidth(i));
	jetSpectra.at(30)->SetBinError(i,jetSpectra.at(30)->GetBinError(i)/jetSpectra.at(30)->GetBinWidth(i));
      }

    jetSpectra.at(30)->GetXaxis()->SetTitle(xSpectraTitle.c_str());
    jetSpectra.at(30)->GetYaxis()->SetTitle(ySpectraTitle.c_str());
    //jetSpectra.at(30)->Scale(24381./19910);
    
    jetSpectra.at(30)->Scale(1./(Lintpb*1e9)/1.4/(2*M_PI));//0.45*1.79769e11));
    jetSpectra.at(30)->GetYaxis()->SetRangeUser(1e-13,1e-4);
  //jetSpectra.at(0)->GetXaxis()->SetTitleSize(0.02);
  //jetSpectra.at(0)->GetYaxis()->SetTitleSize(0.02);
  //jetSpectra.at(0)->GetXaxis()->SetLabelSize();
  //jetSpectra.at(0)->GetYaxis()->SetTitle();
  
    cout << jetSpectra.size() << endl;
    //jetSpectra.at(52)->Rebin(10);
    jetSpectra.at(52)->GetYaxis()->SetTitleOffset(2.0);
    jetSpectra.at(52)->SetMarkerColor(kRed);
    jetSpectra.at(52)->SetMarkerStyle(20);
    jetSpectra.at(52)->SetMarkerSize(2);
    jetSpectra.at(52)->SetLineColor(kRed);
    for(int i=1; i<jetSpectra.at(52)->GetNbinsX()+1; ++i)
      {
	jetSpectra.at(52)->SetBinContent(i,jetSpectra.at(52)->GetBinContent(i)/jetSpectra.at(52)->GetBinWidth(i));
	jetSpectra.at(52)->SetBinError(i,jetSpectra.at(52)->GetBinError(i)/jetSpectra.at(52)->GetBinWidth(i));
      }
    float jet10scalefactor = 41.97/.35/1e7/1.4/(2*M_PI);//*3.646e-3
    jetSpectra.at(52)->GetXaxis()->SetTitle(xSpectraTitle.c_str());
    jetSpectra.at(52)->GetYaxis()->SetTitle("dN_{jet}/dE_{T}"); //ySpectraTitle.c_str());
    jetSpectra.at(52)->Scale(jet10scalefactor);//4e-5/(0.45*2.8e6));//1.11e-11);
    jetSpectra.at(52)->GetYaxis()->SetRangeUser(1e-13,1e-4);
    jetSpectra.at(52)->Draw("PE");
    gPad->SaveAs("output/chi2img/tjetspec.png");

    specspec->GetYaxis()->SetTitleOffset(2.0);
    specspec->SetMarkerColor(kOrange);
    specspec->SetLineColor(kOrange);
    specspec->SetMarkerStyle(20);
    specspec->SetMarkerSize(2);
    specspec->SetLineColor(kOrange);
    for(int i=1; i<specspec->GetNbinsX()+1; ++i)
      {
	specspec->SetBinContent(i,specspec->GetBinContent(i)/specspec->GetBinWidth(i));
	specspec->SetBinError(i,specspec->GetBinError(i)/specspec->GetBinWidth(i));
      }

    specspec->GetXaxis()->SetTitle(xSpectraTitle.c_str());
    //specspec->GetYaxis()->SetRangeUser(1e-13,1e-4);
    specspec->Draw("PE");
    gPad->SaveAs("output/chi2img/specspec.png");
    /*
    hUnfold->GetYaxis()->SetTitleOffset(2.0);
    hUnfold->SetMarkerColor(kRed);
    hUnfold->SetMarkerStyle(20);
    hUnfold->SetMarkerSize(2);
    hUnfold->SetLineColor(kOrange);
    for(int i=1; i<hUnfold->GetNbinsX()+1; ++i)
      {
	hUnfold->SetBinContent(i,hUnfold->GetBinContent(i)/hUnfold->GetBinWidth(i));
      }
    //hUnfold->GetYaxis()->SetRangeUser(1e-13,1e-4);
    hUnfold->GetXaxis()->SetTitle(xSpectraTitle.c_str());
    hUnfold->GetYaxis()->SetTitle(ySpectraTitle.c_str());
    hUnfold->Scale(1./(0.45*1.79769e11));
    hUnfold->Draw("SAME PE");
    gPad->SaveAs("output/chi2img/tjetandunfold.png");
    hUnfold->Draw("PE");
    gPad->SaveAs("output/chi2img/unfold.png");
    */
  jetSpectra.at(0)->Draw("PE");
  cout << "drew 0" << endl;
  canvas->SaveAs("testimg.pdf");
  //sleep(30);
  jetSpectra.at(41)->Scale(jet10scalefactor);//4e-5/(0.45*2.8e6));//1.11e-11);
  specspec->Scale(jet10scalefactor);//4e-5/(0.45*2.8e6));//1.11e-11);
  jetSpectra.at(42)->Scale(jet10scalefactor);//4e-5/(0.45*2.8e6));//1.11e-11);
  TH1* touse[nLegend-1] = {jetSpectra.at(2),jetSpectra.at(16),jetSpectra.at(24),jetSpectra.at(31),jetSpectra.at(41),jetSpectra.at(42)};
  TH1F* prespcr = new TH1F("prespcr","",nbin,bins);
  TH1F* specrat = new TH1F("specrat","",nbin,bins);
  TH1F* djsprat = new TH1F("djsprat","",nbin,bins);

  TH1F* stripcutrat = new TH1F("stripcutrat","",nbin,bins);
  TH1F* etcutrat = new TH1F("etcutrat","",nbin,bins);
  TH1F* allcutrat = new TH1F("allcutrat","",nbin,bins);
  

  //touse[0]->Scale(24381./19910);
  //touse[1]->Scale(24381./19910);
  //touse[2]->Scale(24381./19910);
  //touse[3]->Scale(24381./19910);
  //jetSpectra.at(0)->Scale(24381./19910);
  TLegend* specL2 = new TLegend(*spectraLegend);
  spectraLegend->AddEntry(jetSpectra.at(0),spectraLegendNames[0].c_str(),"p");
  specL2->AddEntry(jetSpectra.at(0),spectraLegendNames[0].c_str(),"p");
  /*
  if(isPythia)
    {
      jetSpectra.at(52)->Scale(4.197e-2/(isJet30?2.505e-9:3.646e-6));
      jetSpectra.at(42)->Scale(4.197e-2/3.646e-6);
      jetSpectra.at(41)->Scale(4.197e-2/3.646e-6);
      //specspec->Scale(4.197e-2/3.646e-6);
    }  
  */
  TLegend* specL3 = new TLegend(0.6,0.55,0.9,0.65);
  specL3->SetBorderSize(0);
  specL3->SetFillStyle(0);

  specL2->SetTextSize(0.015);
  for(int i=1; i<nLegend; ++i)
    {
      if(i==2) continue;

      //if(i==) continue;
      if(i<nLegend-2) touse[i-1]->Scale(1./(Lintpb*1e9)/1.4/(2*M_PI));//0.45*1.79769e11));
      //touse[i-1]->Rebin(10);
      touse[i-1]->GetYaxis()->SetTitleOffset(1.85);
      touse[i-1]->SetMarkerColor(spectraColors[i]);
      touse[i-1]->SetMarkerStyle((i==5)?21:20);
      touse[i-1]->SetMarkerSize(2);
      touse[i-1]->SetLineColor(spectraColors[i]);
      for(int j=1; j<touse[i-1]->GetNbinsX()+1; ++j)
      {
	touse[i-1]->SetBinContent(j,touse[i-1]->GetBinContent(j)/touse[i-1]->GetBinWidth(j));
	touse[i-1]->SetBinError(j,touse[i-1]->GetBinError(j)/touse[i-1]->GetBinWidth(j));
      }
      touse[i-1]->GetYaxis()->SetRangeUser(1e-13,1e-4);
      touse[i-1]->GetXaxis()->SetTitle(xSpectraTitle.c_str());
      touse[i-1]->GetYaxis()->SetTitle(ySpectraTitle.c_str());
      touse[i-1]->Draw("SAME PE");
      cout << "drew " << i << endl;
      if(i<4) spectraLegend->AddEntry(touse[i-1],spectraLegendNames[i].c_str(),"p");
      //if(i==5) continue;
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
  stripcutrat->Divide(jetSpectra.at(2),jetSpectra.at(41));
  etcutrat->Divide(jetSpectra.at(24),jetSpectra.at(41));
  allcutrat->Divide(jetSpectra.at(31),jetSpectra.at(41));
  //specL2->AddEntry(jetSpectra.at(30),"Data All Cuts (No Dijet Check)","p");
  //specL2->AddEntry(specspec,"Reco Sim All Cuts (No Dijet Check)","p");
  specL2->AddEntry(jetSpectra.at(52),(std::string("Truth ")+(isPythia?"PYTHIA":"HERWIG")).c_str(),"p");
  jetSpectra.at(52)->Draw("SAME PE");
  //jetSpectra.at(30)->Draw("SAME PE");
  //specspec->Draw("SAME PE");
  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
  specL2->Draw();
  
  canvas->SaveAs(("output/chi2img/jetSpectrumWithCuts-1_h1.png"));

  jetSpectra.at(42)->Draw("PE");
  jetSpectra.at(41)->Draw("SAME PE");
  jetSpectra.at(52)->Draw("SAME PE");
  specL3->AddEntry(jetSpectra.at(42),"Reco Sim All Cuts","p");
  specL3->AddEntry(jetSpectra.at(41),"Reco Sim No Cuts","p");
  specL3->AddEntry(jetSpectra.at(52),"Truth PYTHIA","p");
  specL3->Draw();
  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
  gPad->SaveAs("output/chi2img/simspeconly.png");

  jetSpectra.at(52)->GetYaxis()->SetRangeUser(0.1,1e6);
  jetSpectra.at(52)->Draw("PE");
  jetSpectra.at(42)->Draw("SAME PE");
  TLegend* tjetleg2 = new TLegend(0.57,0.63,0.9,0.73);
  tjetleg2->SetFillStyle(0);
  tjetleg2->SetFillColor(0);
  tjetleg2->SetTextFont(41);
  tjetleg2->SetBorderSize(0);
  tjetleg2->SetTextSize(0.02);
  tjetleg2->AddEntry(jetSpectra.at(52),(std::string("Truth ")+(isPythia?"PYTHIA":"HERWIG")).c_str(),"p");
  tjetleg2->AddEntry(jetSpectra.at(42),"Reco Sim (With Cuts)","p");
  tjetleg2->Draw();
  texts[0] = "Calorimeter Anti-k_{T} R=0.4";
  texts[1] = "|z_{vtx}| < 30 cm";
  texts[3] = std::string(isPythia?"PYTHIA":"HERWIG")+" "+(isJet30?"Jet-30":"Jet-10")+" Sample Sim";
  //texts[3] = "1-1 Jet Matching Applied";
  /*
  for(int i=3; i<ntext; ++i)
    {
      texts[i] = "";
    }
  */
  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
  canvas->SaveAs((std::string("output/chi2img/")+(isPythia?"PYTHIA":"HERWIG")+"_"+(isJet30?"30":"10")+"_spec.png").c_str());
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
  TH1F* simcutrat = new TH1F("simcutrat","",nbin,bins);
  //jetSpectra.at(42)->Scale(1./(4.197e-2/3.646e-6));
  //jetSpectra.at(42)->Scale(1./(3.646e-3/1e7/1.4));//4e-5/(0.45*2.8e6));//1.11e-11);
  //jetSpectra.at(41)->Scale(1./(4.197e-2/3.646e-6));
  //jetSpectra.at(41)->Scale(1./(3.646e-3/1e7/1.4));//4e-5/(0.45*2.8e6));//1.11e-11);
  //  jetSpectra.at(42)->SetBinContent(6,jetSpectra.at(42)->GetBinContent(6)*jetSpectra.at(42)->GetBinWidth(6)-1);
  //jetSpectra.at(42)->SetBinContent(7,jetSpectra.at(42)->GetBinContent(7)*jetSpectra.at(42)->GetBinWidth(7)-1);
  simcutrat->Divide(jetSpectra.at(42),jetSpectra.at(41),1,1,"B");
  //jetSpectra.at(42)->SetBinContent(6,(jetSpectra.at(42)->GetBinContent(6)+1)/jetSpectra.at(42)->GetBinWidth(6));
  //jetSpectra.at(42)->SetBinContent(7,(jetSpectra.at(42)->GetBinContent(7)+1)/jetSpectra.at(42)->GetBinWidth(7));
  simcutrat->SetMarkerStyle(20);
  simcutrat->SetMarkerSize(2);
  simcutrat->GetYaxis()->SetRangeUser(0.8,1.1);
  simcutrat->GetYaxis()->SetTitle("Reco Sim Cut/No Cut Ratio");
  simcutrat->GetXaxis()->SetTitle("E_{T,jet} [GeV]");
  simcutrat->Draw("PE");
  gPad->SaveAs("output/chi2img/simcutrat.png");

  TH1F* datcutnocutrat = new TH1F("datcutnocutrat","",nbin,bins);
  datcutnocutrat->Divide(jetSpectra.at(31),jetSpectra.at(0),1,1,"B");
  datcutnocutrat->SetMarkerStyle(20);
  datcutnocutrat->SetMarkerSize(2);
  datcutnocutrat->GetYaxis()->SetRangeUser(0,1.1);
  datcutnocutrat->GetYaxis()->SetTitle("Data Cut/No Cut Ratio");
  datcutnocutrat->GetXaxis()->SetTitle("E_{T,jet} [GeV]");
  datcutnocutrat->Draw("PE");
  gPad->SaveAs("output/chi2img/datcutnocutrat.png");

  TH1F* dijet_to_no_dijet_sim = new TH1F("dijet_to_no_dijet_sim","",nbin,bins);
  //specspec->Scale(1./(3.646e-3/1e7/1.4));//4e-5/(0.45*2.8e6));//1.11e-11);
  dijet_to_no_dijet_sim->Divide(specspec,jetSpectra.at(42),1,1,"B");
  //specspec->Scale((3.646e-3/1e7/1.4));//4e-5/(0.45*2.8e6));//1.11e-11);
  dijet_to_no_dijet_sim->SetMarkerStyle(20);
  dijet_to_no_dijet_sim->SetMarkerSize(2);
  dijet_to_no_dijet_sim->GetYaxis()->SetRangeUser(0,1.1);
  dijet_to_no_dijet_sim->GetYaxis()->SetTitle("Reco Sim No Dijet Preservation/Dijet Preservation Ratio");
  dijet_to_no_dijet_sim->GetXaxis()->SetTitle("E_{T,jet} [GeV]");
  dijet_to_no_dijet_sim->Draw("PE");
  gPad->SaveAs("output/chi2img/dijet_to_no_dijet_sim.png");
  //jetSpectra.at(42)->Scale((4.197e-2/3.646e-6));
  //jetSpectra.at(42)->Scale((3.646e-3/1e7/1.4));//4e-5/(0.45*2.8e6));//1.11e-11);
  //jetSpectra.at(41)->Scale((4.197e-2/3.646e-6));
  //jetSpectra.at(41)->Scale((3.646e-3/1e7/1.4));//4e-5/(0.45*2.8e6));//1.11e-11);

  TH1F* dijet_to_no_dijet = new TH1F("dijet_to_no_dijet","",nbin,bins);
  jetSpectra.at(30)->Scale(2*M_PI*Lintpb*1e9*1.4);
  jetSpectra.at(31)->Scale(2*M_PI*Lintpb*1e9*1.4);
  jetSpectra.at(30)->SetBinContent(7,jetSpectra.at(30)->GetBinContent(7)*jetSpectra.at(30)->GetBinWidth(7)-1);
  jetSpectra.at(31)->SetBinContent(7,jetSpectra.at(31)->GetBinContent(7)*jetSpectra.at(31)->GetBinWidth(7)-1);
  //jetSpectra.at(30)->SetBinContent(8,jetSpectra.at(30)->GetBinContent(8)*jetSpectra.at(30)->GetBinWidth(8)-1);
  //jetSpectra.at(30)->SetBinContent(9,jetSpectra.at(30)->GetBinContent(9)*jetSpectra.at(30)->GetBinWidth(9)-1);
  dijet_to_no_dijet->Divide(jetSpectra.at(30),jetSpectra.at(31),1,1,"B");
  jetSpectra.at(30)->SetBinContent(7,(jetSpectra.at(30)->GetBinContent(7)+1)/jetSpectra.at(30)->GetBinWidth(7));
  jetSpectra.at(31)->SetBinContent(7,(jetSpectra.at(31)->GetBinContent(7)+1)/jetSpectra.at(31)->GetBinWidth(7));
  //jetSpectra.at(30)->SetBinContent(8,(jetSpectra.at(30)->GetBinContent(8)+1)/jetSpectra.at(30)->GetBinWidth(8));
  //jetSpectra.at(30)->SetBinContent(9,(jetSpectra.at(30)->GetBinContent(9)+1)/jetSpectra.at(30)->GetBinWidth(9));
  jetSpectra.at(30)->Scale(1./(Lintpb*1e9*1.4)/(2*M_PI));
  jetSpectra.at(31)->Scale(1./(Lintpb*1e9*1.4)/(2*M_PI));
  dijet_to_no_dijet->SetMarkerStyle(20);
  dijet_to_no_dijet->SetMarkerSize(2);
  dijet_to_no_dijet->GetYaxis()->SetRangeUser(0,1.1);
  dijet_to_no_dijet->GetYaxis()->SetTitle("Data No Dijet Preservation/Dijet Preservation Ratio");
  dijet_to_no_dijet->GetXaxis()->SetTitle("E_{T,jet} [GeV]");
  dijet_to_no_dijet->Draw("PE");
  gPad->SaveAs("output/chi2img/dijet_to_no_dijet.png");
  

  //gPad->SetLogy(0);
  specrat->Divide(jetSpectra.at(31),jetSpectra.at(42));
  djsprat->Divide(jetSpectra.at(30),specspec);
  djsprat->SetMarkerColor(kBlue);
  djsprat->SetMarkerStyle(20);
  djsprat->SetLineColor(kBlue);
  djsprat->SetMarkerSize(2);
  specrat->GetYaxis()->SetTitleOffset(2.0);
  specrat->SetMarkerColor(spectraColors[0]);
  specrat->SetMarkerStyle(20);
  specrat->SetMarkerSize(2);
  specrat->SetLineColor(spectraColors[0]);
  specrat->GetYaxis()->SetRangeUser(0,2);
  specrat->GetXaxis()->SetTitle("E_{T,jet} [GeV]");
  specrat->GetYaxis()->SetTitle("Data/Reco Sim All Cuts Spectrum Ratio");
  specrat->Draw("PE");
  //djsprat->Draw("SAME PE");
  stdx = 0.3;
  std_text(canvas, texts, ntext, stdsize, stdx, stdy, stdright);
  stdx = 0.6;
  TLegend* specratleg = new TLegend(0.25,0.25,0.5,0.35);
  specratleg->SetBorderSize(0);
  specratleg->SetFillStyle(0);
  specratleg->AddEntry(specrat,"Dijet Preservation","p");
  specratleg->AddEntry(djsprat,"No Dijet Preservation","p");
  //specratleg->Draw();
  canvas->SaveAs("output/chi2img/datsimrat_h1.png");


  prespcr->Divide(jetSpectra.at(0),jetSpectra.at(41));
  prespcr->GetYaxis()->SetTitleOffset(2.0);
  prespcr->SetMarkerColor(spectraColors[0]);
  prespcr->SetMarkerStyle(21);
  prespcr->SetMarkerSize(2);
  prespcr->SetLineColor(spectraColors[0]);

  stripcutrat->SetMarkerColor(spectraColors[1]);
  stripcutrat->SetMarkerStyle(20);
  stripcutrat->SetMarkerSize(2);
  stripcutrat->SetLineColor(spectraColors[1]);

  etcutrat->SetMarkerColor(spectraColors[3]);
  etcutrat->SetMarkerStyle(20);
  etcutrat->SetMarkerSize(2);
  etcutrat->SetLineColor(spectraColors[3]);

  allcutrat->SetMarkerColor(spectraColors[4]);
  allcutrat->SetMarkerStyle(20);
  allcutrat->SetMarkerSize(2);
  allcutrat->SetLineColor(spectraColors[4]);

  prespcr->GetXaxis()->SetTitle("E_{T,jet} [GeV]");
  prespcr->GetYaxis()->SetTitle("Data/Reco Sim No Cuts Spectrum Ratio");
  prespcr->GetYaxis()->SetRangeUser(0,5);
  prespcr->Draw("PE");
  stripcutrat->Draw("SAME PE");
  etcutrat->Draw("SAME PE");
  allcutrat->Draw("SAME PE");
  spectraLegend->SetX1NDC(0.23);
  spectraLegend->SetX2NDC(0.4);
  spectraLegend->Draw();

  std_text(canvas, texts, ntext, stdsize, 0.23, stdy, stdright);

  canvas->SaveAs("output/chi2img/precutdatsimrat_h1.png");
  gPad->SetBottomMargin(0.1);
  gPad->SetRightMargin(0.1);

  TH1F* xJRat[2];
  xJRat[0] = new TH1F("xjrat0","",100,0,1);
  xJRat[1] = new TH1F("xjrat1","",100,0,1);
  TLegend* xJLeg = new TLegend(0.2,0.4,0.46,0.6);
  xJLeg->SetBorderSize(0);
  xJLeg->SetFillStyle(0);
  int xJColors[5] = {kBlack,kRed+2,kBlue+2,kGreen+2,kOrange+2};
  string xJLegNames[5] = {"Data No Cuts","Data With Cuts","Reco Sim No Cuts","Reco Sim With Cuts",(std::string("Truth ")+(isPythia?"PYTHIA":"HERWIG")).c_str()};
  for(int i=0; i<5; ++i)
    {
      xJ.at(i)->Rebin(5);
      xJ.at(i)->Scale(1./xJ.at(i)->Integral());
      xJ.at(i)->SetLineColor(xJColors[i]);
      xJ.at(i)->SetMarkerColor(xJColors[i]);
      xJ.at(i)->SetMarkerStyle(i%2==1?20:21);
      xJ.at(i)->SetMarkerSize(2);
      xJ.at(i)->GetYaxis()->SetTitleOffset(2.0);
      xJ.at(i)->GetYaxis()->SetTitle("Normalized Counts");
      xJ.at(i)->GetXaxis()->SetTitle("x_{J}");
      xJ.at(i)->GetYaxis()->SetRangeUser(0,0.25);
      if(i==0) xJ.at(i)->Draw("PE");
      else xJ.at(i)->Draw("SAME PE");
      xJLeg->AddEntry(xJ.at(i),xJLegNames[i].c_str(),"p");
    }
  texts[3] = std::string("E_{T,lead} > ")+(isJet30?"30":"14")+" GeV";
  texts[4] = std::string("E_{T,sub} > ")+(isJet30?"20":"10")+" GeV";
  stdx = 0.13;
  std_text(canvas, texts, ntext, stdsize, 0.25, stdy, stdright);
  xJLeg->Draw();
  canvas->SaveAs((std::string("output/chi2img/xJ_h1_")+(isPythia?"PYTHIA_":"HERWIG_")+(isJet30?"30":"10")+".png").c_str());



  
  TH1F* dhCut_to_all = new TH1F("dhCut_to_all","",nbin,bins);
  TH1F* ihCut_to_all = new TH1F("ihCut_to_all","",nbin,bins);
  TH1F* loETCut_to_all = new TH1F("loETCut_to_all","",nbin,bins);
  TH1F* hiETCut_to_all = new TH1F("hiETCut_to_all","",nbin,bins);

  TH1F* dhCut_minus_all = new TH1F("dhCut_minusall","",nbin,bins);
  TH1F* ihCut_minus_all = new TH1F("ihCut_minusall","",nbin,bins);
  TH1F* loETCut_minus_all = new TH1F("loETCut_minusall","",nbin,bins);
  TH1F* hiETCut_minus_all = new TH1F("hiETCut_minusall","",nbin,bins);

  TH1F* dhCutfail_to_all = new TH1F("dhCutfail_to_all","",nbin,bins);
  TH1F* ihCutfail_to_all = new TH1F("ihCutfail_to_all","",nbin,bins);
  TH1F* loETCutfail_to_all = new TH1F("loETCutfail_to_all","",nbin,bins);
  TH1F* hiETCutfail_to_all = new TH1F("hiETCutfail_to_all","",nbin,bins);
  
  //jetSpectra.at(3)->Rebin(10);
  //jetSpectra.at(4)->Rebin(10);
  //jetSpectra.at(5)->Rebin(10);
  //jetSpectra.at(2)->Scale(1./(0.45*1.79769e11));
  jetSpectra.at(3)->Scale(1./(Lintpb*1e9)/(2*M_PI));//0.45*1.79769e11));
  jetSpectra.at(4)->Scale(1./(Lintpb*1e9)/(2*M_PI));//0.45*1.79769e11));
  jetSpectra.at(5)->Scale(1./(Lintpb*1e9)/(2*M_PI));//0.45*1.79769e11));
  dhCut_to_all->Divide(jetSpectra.at(2),jetSpectra.at(0),1,1,"B");
  ihCut_to_all->Divide(jetSpectra.at(3),jetSpectra.at(0),1,1,"B");
  loETCut_to_all->Divide(jetSpectra.at(4),jetSpectra.at(0),1,1,"B");
  hiETCut_to_all->Divide(jetSpectra.at(5),jetSpectra.at(0),1,1,"B");

  
  
  

  dhCut_minus_all->Add(jetSpectra.at(0),jetSpectra.at(2),1,-1);
  ihCut_minus_all->Add(jetSpectra.at(0),jetSpectra.at(3),1,-1);
  loETCut_minus_all->Add(jetSpectra.at(0),jetSpectra.at(4),1,-1);
  hiETCut_minus_all->Add(jetSpectra.at(0),jetSpectra.at(5),1,-1);

  dhCutfail_to_all->Divide(dhCut_minus_all,jetSpectra.at(0));
  ihCutfail_to_all->Divide(ihCut_minus_all,jetSpectra.at(0));
  loETCutfail_to_all->Divide(loETCut_minus_all,jetSpectra.at(0));
  hiETCutfail_to_all->Divide(hiETCut_minus_all,jetSpectra.at(0));
  
  dhCut_minus_all->Draw();
  canvas->SaveAs("output/chi2img/ihcutfail_minus_all.png");



  FormatTH1(dhCut_to_all, "E_{T} [GeV]","Ratio of Passing Spectrum to No Cut Spectrum",spectraColors[0],0,1);
  FormatTH1(ihCut_to_all, "E_{T} [GeV]","Ratio of Passing Spectrum to No Cut Spectrum",spectraColors[1],0,1);
  FormatTH1(loETCut_to_all, "E_{T} [GeV]","Ratio of Passing Spectrum to No Cut Spectrum",spectraColors[2],0,1);
  FormatTH1(hiETCut_to_all, "E_{T} [GeV]","Ratio of Passing Spectrum to No Cut Spectrum",spectraColors[3],0,1);

  FormatTH1(dhCut_minus_all, "E_{T} [GeV]","Failing Cuts #frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",spectraColors[0],1e-13,1e-4);
  FormatTH1(ihCut_minus_all, "E_{T} [GeV]","Failing Cuts #frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",spectraColors[1],1e-11,1e-4);
  FormatTH1(loETCut_minus_all, "E_{T} [GeV]","Failing Cuts #frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",spectraColors[2],1e-11,1e-4);
  FormatTH1(hiETCut_minus_all, "E_{T} [GeV]","Failing Cuts #frac{1}{N_{evt}}#frac{dN_{jet}}{dE_{T,jet}} [GeV^{-1}]",spectraColors[3],1e-11,1e-4);
  FormatTH1(dhCutfail_to_all, "E_{T} [GeV]","Ratio of Failing Spectrum to No Cut Spectrum",spectraColors[0],0,1);
  FormatTH1(ihCutfail_to_all, "E_{T} [GeV]","Ratio of Failing Spectrum to No Cut Spectrum",spectraColors[1],0,1);
  FormatTH1(loETCutfail_to_all, "E_{T} [GeV]","Ratio of Failing Spectrum to No Cut Spectrum",spectraColors[2],0,1);
  FormatTH1(hiETCutfail_to_all, "E_{T} [GeV]","Ratio of Failing Spectrum to No Cut Spectrum",spectraColors[3],0,1);

  TLegend* passing_cut_to_all = new TLegend(0.57,0.55,0.9,0.65);
  TLegend* all_minus_passing = new TLegend(0.57,0.55,0.9,0.65);
  TLegend* failing_cut_to_all = new TLegend(0.57,0.55,0.9,0.65);
  passing_cut_to_all->SetFillStyle(0);
  passing_cut_to_all->SetFillColor(0);
  passing_cut_to_all->SetTextFont(42);
  passing_cut_to_all->SetBorderSize(0);
  passing_cut_to_all->SetTextSize(0.02);

  all_minus_passing->SetFillStyle(0);
  all_minus_passing->SetFillColor(0);
  all_minus_passing->SetTextFont(42);
  all_minus_passing->SetBorderSize(0);
  all_minus_passing->SetTextSize(0.02);

  failing_cut_to_all->SetFillStyle(0);
  failing_cut_to_all->SetFillColor(0);
  failing_cut_to_all->SetTextFont(42);
  failing_cut_to_all->SetBorderSize(0);
  failing_cut_to_all->SetTextSize(0.02);
  gPad->SetLogy(0);
  
  passing_cut_to_all->AddEntry(dhCut_to_all,"Strip Cut","p");
  passing_cut_to_all->AddEntry(ihCut_to_all,"IH Fraction Cut","p");
  passing_cut_to_all->AddEntry(loETCut_to_all,"Low EM Frac. Cut","p");
  passing_cut_to_all->AddEntry(hiETCut_to_all,"High EM Frac. Cut","p");
  
  all_minus_passing->AddEntry(dhCut_minus_all,"Strip Cut","p");
  all_minus_passing->AddEntry(ihCut_minus_all,"IH Fraction Cut","p");
  all_minus_passing->AddEntry(loETCut_minus_all,"Low EM Frac. Cut","p");
  all_minus_passing->AddEntry(hiETCut_minus_all,"High EM Frac. Cut","p");
  gPad->SetLogy(0);
  failing_cut_to_all->AddEntry(dhCutfail_to_all,"Strip Cut","p");
  failing_cut_to_all->AddEntry(ihCutfail_to_all,"IH Fraction Cut","p");
  failing_cut_to_all->AddEntry(loETCutfail_to_all,"Low EM Frac. Cut","p");
  failing_cut_to_all->AddEntry(hiETCutfail_to_all,"High EM Frac. Cut","p");
  

  dhCut_to_all->Draw("P HIST");
  ihCut_to_all->Draw("SAME P HIST");
  loETCut_to_all->Draw("SAME P HIST");
  hiETCut_to_all->Draw("SAME P HIST");
  passing_cut_to_all->Draw();
  canvas->SaveAs("output/chi2img/h1_passing_to_all.png");

  dhCutfail_to_all->Draw("PHIST");
  ihCutfail_to_all->Draw("SAME P HIST");
  loETCutfail_to_all->Draw("SAME P HIST");
  hiETCutfail_to_all->Draw("SAME P HIST");
  failing_cut_to_all->Draw();
  canvas->SaveAs("output/chi2img/h1_failing_to_all.png");

  gPad->SetLogy();
  dhCut_minus_all->Draw("P HIST");
  ihCut_minus_all->Draw("SAME P HIST");
  loETCut_minus_all->Draw("SAME P HIST");
  hiETCut_minus_all->Draw("SAME P HIST");
  all_minus_passing->Draw();
  canvas->SaveAs("output/chi2img/h1_passing_minus_all.png");
  gPad->SetLogy(0);

  const int naxistitles2 = 3;
  std::string xtitles_hists2_sim[naxistitles2] = {"Fraction of E_{T,jet} in EMCal","Fraction of E_{T,jet} in OHCal","Fraction of E_{T,jet} in OHCal"};

  std::string ytitles_hists2_sim[naxistitles2] = {"E_{T,jet}","E_{T,jet}","Fraction of E_{T,jet} in EMCal"};

  std::string filenames2[naxistitles2] = {"output/chi2img/sim_frcem_et","output/chi2img/sim_frcoh_et","output/chi2img/sim_frcoh_frcem"};


  for(int i=0; i<6; ++i)
    {
    std:string cutstring = i%2==0?"_nocuts":"_withcuts";
      std::string texts2[2] = {i%2==0?"No cuts applied":"All cuts applied",""};
      cout << hists2_sim.at(i) << endl;
      FormatAndDrawHistogram(
			     canvas,hists2_sim.at(i),
			     (filenames2[i/2]+cutstring),xtitles_hists2_sim[i/2].c_str(),ytitles_hists2_sim[i/2].c_str(),"N_{jet}",
			     0.15,0.2,0.15,0.15,
			     1.85,1.85,2,
			     0.03,0.03,0.03,
			     0.03,0.03,0.03,
			     texts2
			     );
    }

  int histgroup;
  //return 0;
  for(int i=0; i<(92*3); ++i)//histograms.size(); ++i)
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
      const int ntype = 36;
      string texts[2];
      int threshes[ntype/6] = {8,15,20,25,35,40};
      //int slt[ntype/2] = {8,8,10,10,15,15,20,20,25,25,30,30};
      //texts[0] = "E_{T,jet} > "+to_string(threshes[(i/nhist)%(ntype/2)])+" GeV";
      //texts[0] = "Events contain leading jet with E_{T,jet} > "+to_string(threshes[(i/nhist)%6])+" GeV, "+((i/nhist)%ntype>=24?"inclusive jets":"dijets only");
      if(i<(92))
	{
	  texts[0] = "Leading jets pass nominal event quality cuts";
	  texts[1] = "Leading jets may fail cut without back-to-back preservation";
	}
      else if(i<(92*2))
	{
	  texts[0] = "Leading jets pass nominal event quality cuts";
	  texts[1] = "Leading jets are good even without back-to-back preservation";
	}
      else
	{
	  texts[0] = "Leading jets pass nominal event quality cuts";
	  texts[1] = "Leading jets exclusively FAIL without back-to-back preservation";
	}
      /*
      else
	{
	  texts[0] = "Events contain at least one jet with E_{T,jet} > 8 GeV";//+to_string(threshes[(i/nhist)%(ntype/2)])+" GeV";
	}
      */
	/*
	if(i<(ntype/4)*nhist) texts[1] = "And at least one subleading jet with E_{T,jet} > "+to_string(slt[(i/nhist)%(ntype/2)])+" GeV";
	else
	{
	  texts[1] = "And no subleading jet with E_{T,jet} > "+to_string(slt[(i/nhist)%(ntype/2)])+" GeV";
	}
	*/
      /*
      if((i/nhist)%ntype < ntype/6) texts[1] = "All jet quality cuts applied, #Delta#phi < 3#pi/4, no EM fraction cut";
      else if((i/nhist)%ntype < ntype/3) texts[1] = "All jet quality cuts applied, #Delta#phi > 3#pi/4, no EM fraction cut";
      else if((i/nhist)%ntype < ntype/2) texts[1] = "All jet quality cuts applied, #Delta#phi_{closest non-lead} > #pi/4, no EM fraction cut";
      else if((i/nhist)%ntype < 2*ntype/3) texts[1] = "All jet quality cuts applied,  #Delta#phi_{closest non-lead} < #pi/4, no EM fraction cut";
      else if((i/nhist)%ntype < 5*ntype/6) texts[1] = "All jet quality cuts applied, E_{T,EM}/E_{T,jet} > 0.9";
      else texts[1]="All jet quality cuts applied, E_{T,EM}/E_{T,jet} < 0.9";
      cout << "prep to draw" << endl;
      //histograms.at(i)->Rebin2D(4,4);
      */
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


  for(int i=0; i<36; ++i)//corPlotsLJet.size(); ++i)
    {
      const int ntype = 36;
      cout << corPlotsLJet.at(i)->GetName() << endl;
      //corPlotsLJet.at(i)->Scale(1/corPlotsLJet.at(i)->Integral("WIDTH"));
      //cout << "scaled hist" << endl;
      int threshes[ntype/6] = {8,15,20,25,35,40};      
      string outdrawname = "output/chi2img/"+string(corPlotsLJet.at(i)->GetName()) + "_" + to_string(i/6);
      //int slt[ntype/2] = {8,8,10,10,15,15,20,20,25,25,30,30};
      //texts[0] = "E_{T,jet} > "+to_string(threshes[(i/nhist)%(ntype/2)])+" GeV";
      /*
      texts[0] = "Events contain leading jet with E_{T,jet} > "+to_string(threshes[(i/6)%6])+" GeV, "+((i/6)%ntype>=24?"inclusive jets":"dijets only");
      if((i/6)%ntype < ntype/6) texts[1] = "All jet quality cuts applied, #Delta#phi < 3#pi/4, no EM fraction cut";
      else if((i/6)%ntype < ntype/3) texts[1] = "All jet quality cuts applied, #Delta#phi > 3#pi/4, no EM fraction cut";
      else if((i/6)%ntype < ntype/2) texts[1] = "All jet quality cuts applied,  #Delta#phi_{closest non-lead} > #pi/4, no EM fraction cut";
      else if((i/6)%ntype < 2*ntype/3) texts[1] = "All jet quality cuts applied,  #Delta#phi_{closest non-lead} < #pi/4, no EM fraction cut";
      else if((i/6)%ntype < 5*ntype/6) texts[1] = "All jet quality cuts applied, E_{T,EM}/E_{T,jet} > 0.9";
      else texts[1]="All jet quality cuts applied, E_{T,EM}/E_{T,jet} > 0.9";
      */
      //if(i<(ntype/2)*nhist) corPlotsLJet.at(i)->Add(corPlotsLJet.at(i),corPlotsLJet.at(i+180));
      cout << "prep to draw" << endl;
      //corPlotsLJet.at(i)->Rebin2D(4,4);
      if(i<18)
	{
	  texts[0] = "Reco sim jets; leading jets pass nominal event quality cuts"; 
	}
      else
	{
	  texts[0] = "Truth jets; leading reco jets pass nominal event quality cuts";
	}
      if(i%18<6)
	{
	  texts[1] = "Reco sim jets can either pass or fail without back-to-back preservation";
	}
      else if(i%18<12)
	{
	  texts[1] = "Reco sim jets must pass even without back-to-back preservation";
	}
      else
	{
	  texts[1] = "Reco sim jets must fail without back-to-back preservation";
	}
      FormatAndDrawHistogram(
			     canvas,corPlotsLJet.at(i),
			     outdrawname,xtcp[i%6].c_str(),ytcp[i%6].c_str(),"N_{jet}",
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
