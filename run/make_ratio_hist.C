int make_ratio_hist(string jet10file, string jet30file, string mbfile)
{
  TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const int nbiny = 12;
  float binsy[nbiny+1] = {12, 15, 18, 21, 24, 28, 32, 36, 40, 45, 50, 55, 60};

  TFile* jet10tf = TFile::Open(jet10file.c_str(),"READ");
  TFile* jet30tf = TFile::Open(jet30file.c_str(),"READ");
  TFile* mbtf = TFile::Open(mbfile.c_str(),"READ");

  TFile* outfile = TFile::Open("ratiofile.root","RECREATE");
  
  float jet10scale = 3.646e-6/4.197e-2;
  float jet30scale = 2.505e-9/4.197e-2;
  float mbscale = 0.1;

  TH1D* jet10_ucspec = (TH1D*)jet10tf->Get("h1_ucspec");
  TH1D* jet10_cspec = (TH1D*)jet10tf->Get("h1_cspec");
  TH1D* jet10_specspec = (TH1D*)jet10tf->Get("specspec");

  TH1D* jet30_ucspec = (TH1D*)jet30tf->Get("h1_ucspec");
  TH1D* jet30_cspec = (TH1D*)jet30tf->Get("h1_cspec");
  TH1D* jet30_specspec = (TH1D*)jet30tf->Get("specspec");

  TH1D* mb_ucspec = (TH1D*)mbtf->Get("h1_ucspec");
  TH1D* mb_cspec = (TH1D*)mbtf->Get("h1_cspec");
  TH1D* mb_specspec = (TH1D*)mbtf->Get("specspec");

  TH1D* full_ucspec = new TH1D("full_ucspec","full_ucspec",nbiny,binsy);
  full_ucspec->Add(mb_ucspec);
  full_ucspec->Add(jet10_ucspec);
  full_ucspec->Add(jet30_ucspec);

  jet10_ucspec->Scale(1./jet10scale);
  jet10_cspec->Scale(1./jet10scale);
  jet30_ucspec->Scale(1./jet30scale);
  jet30_cspec->Scale(1./jet30scale);
  mb_ucspec->Scale(1./mbscale);
  mb_cspec->Scale(1./mbscale);

  TH1D* jet10_cspec_m = new TH1D("jet10_cspec_m","jet10_cspec_m",nbiny,binsy);
  TH1D* jet10_specspec_m = new TH1D("jet10_specspec_m","jet10_specspec_m",nbiny,binsy);
  TH1D* jet30_cspec_m = new TH1D("jet30_cspec_m","jet30_cspec_m",nbiny,binsy);
  TH1D* jet30_specspec_m = new TH1D("jet30_specspec_m","jet30_specspec_m",nbiny,binsy);
  TH1D* mb_cspec_m = new TH1D("mb_cspec_m","mb_cspec_m",nbiny,binsy);
  TH1D* mb_specspec_m = new TH1D("mb_specspec_m","mb_specspec_m",nbiny,binsy);

  for(int i=1; i<nbiny+1; ++i)
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
  jet10_specspec->Scale(jet10scale);
  jet10_specspec_m->Scale(jet10scale);

  mb_cspec->Scale(mbscale);
  mb_cspec_m->Scale(mbscale);
  mb_specspec->Scale(mbscale);
  mb_specspec_m->Scale(mbscale);

  jet30_cspec->Scale(jet30scale);
  jet30_cspec_m->Scale(jet30scale);
  jet30_specspec->Scale(jet30scale);
  jet30_specspec_m->Scale(jet30scale);

  TH1D* full_cspec = new TH1D("full_cspec","full_cspec",nbiny,binsy);
  TH1D* full_cspec_m = new TH1D("full_cspec_m","full_cspec_m",nbiny,binsy);
  TH1D* full_specspec = new TH1D("full_specspec","full_specspec",nbiny,binsy);
  TH1D* full_specspec_m = new TH1D("full_specspec_m","full_specspec_m",nbiny,binsy);
			       
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
  
  TH1D* ratio_cspec = new TH1D("ratio_cspec","ratio_cspec",nbiny,binsy);
  TH1D* ratio_cspec_m = new TH1D("ratio_cspec_m","ratio_cspec_m",nbiny,binsy);
  TH1D* ratio_specspec = new TH1D("ratio_specspec","ratio_specspec",nbiny,binsy);
  TH1D* ratio_specspec_m = new TH1D("ratio_specspec_m","ratio_specspec_m",nbiny,binsy);
  
  ratio_cspec->Divide(full_cspec,full_ucspec,1,1,"B");
  ratio_cspec_m->Divide(full_cspec_m,full_ucspec,1,1,"B");
  ratio_specspec->Divide(full_specspec,full_ucspec,1,1,"B");
  ratio_specspec_m->Divide(full_specspec_m,full_ucspec,1,1,"B");
  
  TH1D* ratio_cspec_final = new TH1D("ratio_cspec_final","ratio_cspec_final",nbiny,binsy);
  TH1D* ratio_specspec_final = new TH1D("ratio_specspec_final","ratio_specspec_final",nbiny,binsy);

  for(int i=1; i<nbiny+1; ++i)
    {
      ratio_cspec_final->SetBinContent(i,ratio_cspec->GetBinContent(i));
      ratio_cspec_final->SetBinError(i,(i==nbiny?ratio_cspec_m:ratio_cspec)->GetBinError(i));
      cout << ratio_cspec_final->GetBinContent(i) << endl;
      ratio_specspec_final->SetBinContent(i,ratio_specspec->GetBinContent(i));
      ratio_specspec_final->SetBinError(i,(i==nbiny?ratio_specspec_m:ratio_specspec)->GetBinError(i));
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


  outfile->Write();
  outfile->Close();
  jet10tf->Close();
  jet30tf->Close();
  mbtf->Close();

  return 0;
}
  
