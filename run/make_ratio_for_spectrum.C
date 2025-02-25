int make_ratio_for_spectrum()

{
  TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const int nbiny = 12;
  double binsy[nbiny+1] = {12, 15, 18, 21, 24, 28, 32, 36, 40, 45, 50, 55, 60};

  TFile* jcf = TFile::Open("chi2file_allhists.root","READ");
  TFile* hpf = TFile::Open("output_data.root","READ");

  TH1D* jcs = (TH1D*) jcf->Get("h1_jetSpectra_31");
  //TH1D* hp1 = (TH1D*) hpf->Get("h_recojet_pt_record");
  TH1D* hp2 = (TH1D*) hpf->Get("h_recojet_pt");

  //hp1 = (TH1D*)hp1->Rebin(nbiny,"",binsy);
  //hp2->Rebin(nbiny,"",binsy);
  //cout << hp2->GetNbinsX() << " " << hp1->GetNbinsX() << endl;
  cout << jcs->GetNbinsX() << endl;
  TH1D* recordrat = new TH1D("recordrat","recordrat",nbiny,binsy);
  TH1D* rat = new TH1D("rat","rat",nbiny,binsy);

  jcs->Draw("PE");
  gPad->SaveAs("output/chi2img/jcs1.png");

  rat->Divide(hp2,jcs);  
  //recordrat->Divide(hp1,jcs);
  for(int i=1; i<nbiny+1; ++i)
    {
      cout << rat->GetBinContent(i) << endl;
      jcs->SetBinContent(i,jcs->GetBinContent(i)/jcs->GetBinWidth(i));
      jcs->SetBinContent(i,jcs->GetBinError(i)/jcs->GetBinWidth(i));
    }


  hp2->Draw("PE");
  gPad->SaveAs("output/chi2img/hp2.png");

  //hp1->Draw("PE");
  //gPad->SaveAs("output/chi2img/hp1.png");
  
  jcs->Draw("PE");
  gPad->SaveAs("output/chi2img/jcs2.png");

  recordrat->SetMarkerStyle(20);
  recordrat->GetXaxis()->SetTitle("p_{T}^jet [GeV]");
  recordrat->GetYaxis()->SetRangeUser(0.75,1.1);
  recordrat->Draw("PE");
  gPad->SaveAs("output/chi2img/recordrat.png");

  rat->SetMarkerStyle(20);
  rat->GetXaxis()->SetTitle("p_{T}^jet [GeV]");
  rat->GetYaxis()->SetRangeUser(0.75,1.1);
  rat->Draw("PE");
  gPad->SaveAs("output/chi2img/rat.png");
  

  return 0;
}
  
