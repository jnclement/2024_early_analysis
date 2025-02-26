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

int make_ratio_hist(string jet10file, string jet30file, string mbfile)
{
  TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  const int nbiny = 12;
  float binsy[nbiny+1] = {12, 15, 18, 21, 24, 28, 32, 36, 40, 45, 50, 55, 60};
  TFile* datf = TFile::Open("chi2file_allhists.root","READ");
  TH1D* dat_cspec = (TH1D*)datf->Get("fullRangeSpectra_30");
  TH1D* dat_ucspec = (TH1D*)datf->Get("fullRangeSpectra_0");
  TH1D* dat_spspec = (TH1D*)datf->Get("fullRangeSpectra_31");
  
  dat_cspec->Rebin(20);
  dat_ucspec->Rebin(20);
  dat_spspec->Rebin(20);

  TFile* jet10tf = TFile::Open(jet10file.c_str(),"READ");
  TFile* jet30tf = TFile::Open(jet30file.c_str(),"READ");
  TFile* mbtf = TFile::Open(mbfile.c_str(),"READ");

  TFile* outfile = TFile::Open("ratiofile.root","RECREATE");
  
  float jet10scale = 3.646e-6/4.197e-2;
  float jet30scale = 2.505e-9/4.197e-2;
  float mbscale = 0.1;

  TH1D* jet10_ucspec = (TH1D*)jet10tf->Get("fullrange_ucspec");
  TH1D* jet10_cspec = (TH1D*)jet10tf->Get("fullrange_cspec");
  TH1D* jet10_specspec = (TH1D*)jet10tf->Get("fullrange_specspec");

  TH1D* jet30_ucspec = (TH1D*)jet30tf->Get("fullrange_ucspec");
  TH1D* jet30_cspec = (TH1D*)jet30tf->Get("fullrange_cspec");
  TH1D* jet30_specspec = (TH1D*)jet30tf->Get("fullrange_specspec");

  TH1D* mb_ucspec = (TH1D*)mbtf->Get("fullrange_ucspec");
  TH1D* mb_cspec = (TH1D*)mbtf->Get("fullrange_cspec");
  TH1D* mb_specspec = (TH1D*)mbtf->Get("fullrange_specspec");
  
  jet10_cspec->Rebin(20);
  jet10_ucspec->Rebin(20);
  jet10_specspec->Rebin(20);
  jet30_cspec->Rebin(20);
  jet30_ucspec->Rebin(20);
  jet30_specspec->Rebin(20);
  mb_cspec->Rebin(20);
  mb_ucspec->Rebin(20);
  mb_specspec->Rebin(20);

  TH1D* full_ucspec = new TH1D("full_ucspec","full_ucspec",50,0,100);//nbiny,binsy);
  full_ucspec->Add(mb_ucspec);
  full_ucspec->Add(jet10_ucspec);
  full_ucspec->Add(jet30_ucspec);

  jet10_ucspec->Scale(1./jet10scale);
  jet10_cspec->Scale(1./jet10scale);
  jet30_ucspec->Scale(1./jet30scale);
  jet30_cspec->Scale(1./jet30scale);
  mb_ucspec->Scale(1./mbscale);
  mb_cspec->Scale(1./mbscale);
  
  TH1D* jet10_cspec_m = new TH1D("jet10_cspec_m","jet10_cspec_m",50,0,100);//nbiny,binsy);
  TH1D* jet10_specspec_m = new TH1D("jet10_specspec_m","jet10_specspec_m",50,0,100);//nbiny,binsy);
  TH1D* jet30_cspec_m = new TH1D("jet30_cspec_m","jet30_cspec_m",50,0,100);//nbiny,binsy);
  TH1D* jet30_specspec_m = new TH1D("jet30_specspec_m","jet30_specspec_m",50,0,100);//nbiny,binsy);
  TH1D* mb_cspec_m = new TH1D("mb_cspec_m","mb_cspec_m",50,0,100);//nbiny,binsy);
  TH1D* mb_specspec_m = new TH1D("mb_specspec_m","mb_specspec_m",50,0,100);//nbiny,binsy);

  for(int i=1; i<100; ++i)
    {
      if(std::abs(jet10_cspec->GetBinContent(i) - jet10_ucspec->GetBinContent(i)) < 1 && jet10_cspec->GetBinContent(i) - 1 > 0)
	{
	  jet10_cspec_m->SetBinContent(i,jet10_cspec->GetBinContent(i)-1);
	  if(i>1 && i<6) jet10_cspec_m->SetBinError(i,sqrt(jet10_cspec_m->GetBinContent(i)));
	}
      else
	{
	  jet10_cspec_m->SetBinContent(i,jet10_cspec->GetBinContent(i));
	  if(i>1 && i<6) jet10_cspec_m->SetBinError(i,sqrt(jet10_cspec_m->GetBinContent(i)));
	}
      if(std::abs(jet10_specspec->GetBinContent(i) - jet10_ucspec->GetBinContent(i)) < 1 && jet10_specspec->GetBinContent(i) - 1 > 0)
	{
	  jet10_specspec_m->SetBinContent(i,jet10_specspec->GetBinContent(i)-1);
	  if(i>1 && i<6) jet10_specspec_m->SetBinError(i,sqrt(jet10_specspec_m->GetBinContent(i)));
	}
      else
	{
	  jet10_specspec_m->SetBinContent(i,jet10_specspec->GetBinContent(i));
	  if(i>1 && i<6) jet10_specspec_m->SetBinError(i,sqrt(jet10_specspec_m->GetBinContent(i)));
	}
      if(std::abs(jet30_cspec->GetBinContent(i) - jet30_ucspec->GetBinContent(i)) < 1 && jet30_cspec->GetBinContent(i) - 1 > 0)
	{
	  jet30_cspec_m->SetBinContent(i,jet30_cspec->GetBinContent(i)-1);
	  if(i>5) jet30_cspec_m->SetBinError(i,sqrt(jet30_cspec_m->GetBinContent(i)));
	}
      else
	{
	  jet30_cspec_m->SetBinContent(i,jet30_cspec->GetBinContent(i));
	  if(i>5) jet30_cspec_m->SetBinError(i,sqrt(jet30_cspec_m->GetBinContent(i)));
	}
      if(std::abs(jet30_specspec->GetBinContent(i) - jet30_ucspec->GetBinContent(i)) < 1 && jet30_specspec->GetBinContent(i) - 1 > 0)
	{
	  jet30_specspec_m->SetBinContent(i,jet30_specspec->GetBinContent(i)-1);
	  if(i>5) jet30_specspec_m->SetBinError(i,sqrt(jet30_specspec_m->GetBinContent(i)));
	}
      else
	{
	  jet30_specspec_m->SetBinContent(i,jet30_specspec->GetBinContent(i));
	  if(i>5) jet30_specspec_m->SetBinError(i,sqrt(jet30_specspec_m->GetBinContent(i)));
	}

      if(std::abs(mb_cspec->GetBinContent(i) - mb_ucspec->GetBinContent(i)) < 1 && mb_cspec->GetBinContent(i) - 1 > 0)
	{
	  mb_cspec_m->SetBinContent(i,mb_cspec->GetBinContent(i)-1);
	  if(i==1) mb_cspec_m->SetBinError(i,sqrt(mb_cspec_m->GetBinContent(i)));
	}
      else
	{
	  mb_cspec_m->SetBinContent(i,mb_cspec->GetBinContent(i));
	  if(i==1) mb_cspec_m->SetBinError(i,sqrt(mb_cspec_m->GetBinContent(i)));
	}
      if(std::abs(mb_specspec->GetBinContent(i) - mb_ucspec->GetBinContent(i)) < 1 && mb_specspec->GetBinContent(i) - 1 > 0)
	{
	  mb_specspec_m->SetBinContent(i,mb_specspec->GetBinContent(i)-1);
	  if(i==1) mb_specspec_m->SetBinError(i,sqrt(mb_specspec_m->GetBinContent(i)));
	}
      else
	{
	  mb_specspec_m->SetBinContent(i,mb_specspec->GetBinContent(i));
	  if(i==1) mb_specspec_m->SetBinError(i,sqrt(mb_specspec_m->GetBinContent(i)));
	}
    }

  jet10_cspec->Scale(jet10scale);
  jet10_cspec_m->Scale(jet10scale);
  jet10_ucspec->Scale(jet10scale);
  //jet10_ucspec_m->Scale(jet10scale);

  mb_cspec->Scale(mbscale);
  mb_cspec_m->Scale(mbscale);
  mb_ucspec->Scale(mbscale);
  //mb_ucspec_m->Scale(mbscale);

  jet30_cspec->Scale(jet30scale);
  jet30_cspec_m->Scale(jet30scale);
  jet30_ucspec->Scale(jet30scale);
  //jet30_ucspec_m->Scale(jet30scale);

  TH1D* full_cspec = new TH1D("full_cspec","full_cspec",50,0,100);//nbiny,binsy);
  TH1D* full_cspec_m = new TH1D("full_cspec_m","full_cspec_m",50,0,100);//nbiny,binsy);
  TH1D* full_specspec = new TH1D("full_specspec","full_specspec",50,0,100);//nbiny,binsy);
  TH1D* full_specspec_m = new TH1D("full_specspec_m","full_specspec_m",50,0,100);//nbiny,binsy);
			       
  full_cspec->Add(mb_cspec);
  full_cspec->Add(jet10_cspec);
  full_cspec->Add(jet30_cspec);

  full_cspec_m->Add(mb_cspec_m);
  full_cspec_m->Add(jet10_cspec_m);
  full_cspec_m->Add(jet30_cspec_m);

  full_specspec->Add(mb_specspec);
  full_specspec->Add(jet10_specspec);
  full_specspec->Add(jet30_specspec);

  full_specspec_m->Add(mb_specspec_m);
  full_specspec_m->Add(jet10_specspec_m);
  full_specspec_m->Add(jet30_specspec_m);
  
  TH1D* ratio_cspec = new TH1D("ratio_cspec","ratio_cspec",50,0,100);//nbiny,binsy);
  TH1D* ratio_cspec_m = new TH1D("ratio_cspec_m","ratio_cspec_m",50,0,100);//nbiny,binsy);
  TH1D* ratio_specspec = new TH1D("ratio_specspec","ratio_specspec",50,0,100);//nbiny,binsy);
  TH1D* ratio_specspec_m = new TH1D("ratio_specspec_m","ratio_specspec_m",50,0,100);//nbiny,binsy);
  
  ratio_cspec->Divide(full_cspec,full_ucspec,1,1,"B");
  ratio_cspec_m->Divide(full_cspec_m,full_ucspec,1,1,"B");
  ratio_specspec->Divide(full_specspec,full_ucspec,1,1,"B");
  ratio_specspec_m->Divide(full_specspec_m,full_ucspec,1,1,"B");
  
  TH1D* ratio_cspec_final = new TH1D("ratio_cspec_final","ratio_cspec_final",50,0,100);//nbiny,binsy);
  TH1D* ratio_specspec_final = new TH1D("ratio_specspec_final","ratio_specspec_final",50,0,100);//nbiny,binsy);

  TH1D* ratio_cspec_specspec = new TH1D("ratio_cspec_specspec","ratiocspecspecspec",50,0,100);

  full_cspec->Draw("PE");
  full_specspec->SetLineColor(kRed);
  full_specspec->Draw("SAME PE");
  c1->SaveAs("output/chi2img/issue.png");
  
  ratio_cspec_specspec->Divide(full_cspec,full_specspec,1,1,"B");
  
  for(int i=1; i<100;++i)//nbiny+1; ++i)
    {
      ratio_cspec_final->SetBinContent(i,ratio_cspec->GetBinContent(i));
      ratio_cspec_final->SetBinError(i,(i==100?ratio_cspec_m:ratio_cspec)->GetBinError(i));
      cout << ratio_cspec_final->GetBinContent(i) << endl;
      ratio_specspec_final->SetBinContent(i,ratio_specspec->GetBinContent(i));
      ratio_specspec_final->SetBinError(i,(i==100?ratio_specspec_m:ratio_specspec)->GetBinError(i));
      //cout << ratio_specspec_final->GetBinContent(i) << endl;
    }

  TGraphAsymmErrors* ratio_cspec_final_asym = new TGraphAsymmErrors(full_cspec,full_ucspec,"b(1,1)");

  ratio_cspec_final->SetMarkerStyle(20);
  ratio_cspec_final->SetMarkerSize(2);
  ratio_cspec_final->GetYaxis()->SetRangeUser(0.9,1.05);
  ratio_cspec_final->GetYaxis()->SetTitle("Reco Sim Cut/No Cut Ratio TH1::Divide()");
  ratio_cspec_final->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_cspec_final->Draw("PE");
  gPad->SaveAs("output/chi2img/ratio_cspec_final.png");
  ratio_cspec_final_asym->SetMarkerStyle(20);
  ratio_cspec_final_asym->SetMarkerSize(2);
  ratio_cspec_final_asym->GetYaxis()->SetRangeUser(0.9,1.05);
  ratio_cspec_final_asym->GetYaxis()->SetTitle("Reco Sim Cut/No Cut Ratio TGraphAsymmErrors::Divide()");
  ratio_cspec_final_asym->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  ratio_cspec_final_asym->Draw("APE");
  gPad->SaveAs("output/chi2img/ratio_cspec_final_asym.png");

  
  
  float Lintpb = 16.82;
  const int ntext = 5;
  string texts[ntext] =
  {
    "Calorimeter Anti-k_{T} R=0.4",
    "|z_{vtx}| < 30 cm",
    "Jet-10 & MBDNS>=1 Triggered Data",
    "\\mathscr{L}_{\\text{data}}=16.82 \\text{pb}^{-1}",
    "Jet Detector Level |#eta|<0.7"
  };

  ratioPanelCanvas(c1, 0.3);

  full_cspec->Scale(1./jet30scale);
  full_ucspec->Scale(1./jet30scale);
  full_specspec->Scale(1./jet30scale);
  
  TF1* myexp = new TF1("myexp","[0]*exp([1]*x)",0,100);
  TF1* smexp = new TF1("simexp","[0]*exp([1]*x)",0,100);
  smexp->SetParameter(1,-0.3);
  dat_cspec->Fit(myexp,"IL0","",20,50);
  
  full_specspec->Fit(smexp,"I WL 0","",20,50);
  myexp->SetLineColor(kBlue);
  smexp->SetLineColor(kBlue);
  TH1D* cspec_ratio = (TH1D*)dat_spspec->Clone();
  TH1D* scspec_ratio = (TH1D*)full_cspec->Clone();
  cspec_ratio->Sumw2();
  //scspec_ratio->Sumw2();
  cspec_ratio->Divide(myexp);
  scspec_ratio->Divide(smexp);
 
  //ratio_cspec_specspec->GetXaxis()->SetTitleSize(0.1);
  //ratio_cspec_specspec->GetXaxis()->SetLabelSize(0.1);
  //ratio_cspec_specspec->GetYaxis()->SetLabelSize(0.05);
  //ratio_cspec_specspec->GetYaxis()->SetTitleSize(0.1);
  c1->cd(1);
  gPad->SetLogy();
  FormatTH1(dat_ucspec, "E_{T}^{jet} [GeV]", "N_{jet}", kBlack, 1e-12, 1e-4);
  FormatTH1(dat_cspec, "E_{T}^{jet} [GeV]", "N_{jet}",kAzure+1,1e-12,1e-4);
  FormatTH1(dat_spspec, "E_{T}^{jet} [GeV]", "N_{jet}",kRed+1,1e-12,1e-4);
  FormatTH1(scspec_ratio,"E_{T}^{jet} [GeV]", "Ratio",kGreen+2,0,2);
  //dat_cspec->Scale(1./(Lintpb*1e9)/1.4/(2*M_PI));
  //dat_ucspec->Scale(1./(Lintpb*1e9)/1.4/(2*M_PI));
  dat_ucspec->GetYaxis()->SetRangeUser(0.5,1e7);
  dat_ucspec->Draw("PE");
  dat_spspec->Draw("SAME PE");
  dat_cspec->Draw("SAME PE");
  myexp->Draw("SAME");
  //ucexp->Draw("SAME");
  texts[3] = "";
  texts[4]= "";
  std_text(c1, texts, ntext, 0.025, 0.45, 0.85, 0);
  TLegend* frcemleg = new TLegend(0.5,0.53,0.9,0.7);
  frcemleg->SetFillStyle(0);
  frcemleg->SetFillColor(0);
  frcemleg->SetBorderSize(0);
  frcemleg->AddEntry(dat_cspec,"Dijet Cut","p");
  frcemleg->AddEntry(dat_ucspec,"No Cuts","p");
  frcemleg->AddEntry(dat_spspec,"Frac. Cut","p");
  frcemleg->AddEntry(myexp,"Fit Dijet Cut","l");

  frcemleg->Draw();

  c1->cd(2);
  FormatTH1(scspec_ratio, "E_{T}^{jet} [GeV]", "Ratio", kGreen+2,0,2);
  FormatTH1(cspec_ratio, "E_{T}^{jet} [GeV]", "Ratio",kRed+1,0,2);
  //ucspec_ratio->GetYaxis()->CenterTitle();
  cspec_ratio->GetXaxis()->SetLabelSize(0.1);
  cspec_ratio->GetYaxis()->SetLabelSize(0.05);
  cspec_ratio->GetXaxis()->SetTitleSize(0.1);
  cspec_ratio->Draw("PE");

  TLine* line = new TLine(0,1,100,1);
  line->SetLineColor(kBlue);
  line->Draw();

  c1->SaveAs("output/chi2img/counts_cuts.png");  
  scspec_ratio->Draw("PE");
  scspec_ratio->GetXaxis()->SetLabelSize(0.1);
  scspec_ratio->GetYaxis()->SetLabelSize(0.05);
  scspec_ratio->GetXaxis()->SetTitleSize(0.1);
  c1->cd(2);

  TLegend* newleg = new TLegend(0.5,0.53,0.9,0.7);
  texts[2] = "Reco PYTHIA Combined Jet10&Jet30&MB";
  FormatTH1(scspec_ratio, "E_{T}^{jet} [GeV]", "Ratio", kGreen+2,0,2);
  FormatTH1(full_cspec, "E_{T}^{jet} [GeV]", "N_{jet}", kGreen+2,0.5,1e10);
  FormatTH1(full_ucspec,"E_{T}^{jet} [GeV]", "N_{jet}",kBlack, 0.5,1e10);
  FormatTH1(full_specspec,"E_{T}^{jet} [GeV]", "N_{jet}", kAzure+1,0.5,1e10);
  newleg->AddEntry(full_cspec,"Frac. Cut","p");
  newleg->AddEntry(full_ucspec,"No Cuts","p");
  newleg->AddEntry(full_specspec,"Dijet Cut","p");
  newleg->AddEntry(smexp,"Fit Dijet Cut","l");
  smexp->SetLineColor(kBlue);


  c1->cd(1);
  gPad->Clear();
  ratioPanelCanvas(c1, 0.3);
  c1->cd(1);
  gPad->SetLogy();
  full_ucspec->Draw("PE");
  full_specspec->Draw("SAME PE");
  full_cspec->Draw("SAME PE");
  smexp->Draw("SAME");
  newleg->SetFillStyle(0);
  newleg->SetFillColor(0);
  newleg->SetBorderSize(0);
  c1->cd(0);
  newleg->Draw();
  c1->cd(1);
  std_text(c1, texts, ntext, 0.025, 0.45, 0.85, 0);

  c1->cd(2);
  scspec_ratio->Draw("PE");
  line->Draw();
  c1->SaveAs("output/chi2img/count_cuts_sim.png");
  
  TFile* jamiefile = new TFile("output/chi2img/forJamie.root","RECREATE");
  jamiefile->cd();
  dat_cspec->Write();
  jamiefile->Write();
  
  outfile->Write();
  outfile->Close();
  jet10tf->Close();
  jet30tf->Close();
  mbtf->Close();
  jamiefile->Close();

  return 0;
}
  
