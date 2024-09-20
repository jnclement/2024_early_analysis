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
void get_timestamps(int runnumber, unsigned int* timestamp)
{

  string datestamp;
  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","phnxro","");

  if (db)
    {
      printf("Server info: %s\n", db->ServerInfo());
    }
  else
    {
      printf("bad\n");
    }


  TSQLRow *row;
  TSQLResult *res;
  TString cmd = "";
  char sql[1000];

  cout << runnumber << endl;
  sprintf(sql, "select brtimestamp from run where runnumber = %d;", runnumber);

  res = db->Query(sql);

  int nrows = res->GetRowCount();
  
  int nfields = res->GetFieldCount();
  for (int i = 0; i < nrows; i++) {
    row = res->Next();
    for (int j = 0; j < nfields; j++) {
      datestamp = row->GetField(j);
    }
    delete row;
  }
  delete res;
  cout << datestamp << endl;
  TDatime das = TDatime(datestamp.c_str());
  *timestamp = das.Convert();
}

void get_scaledowns(int runnumber, int scaledowns[])
{

  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","phnxro","");

  if (db)
    {
      printf("Server info: %s\n", db->ServerInfo());
    }
  else
    {
      printf("bad\n");
    }


  TSQLRow *row;
  TSQLResult *res;
  TString cmd = "";
  char sql[1000];

  cout << runnumber << endl;
  for (int is = 0; is < 64; is++)
    {
      sprintf(sql, "select scaledown%02d from gl1_scaledown where runnumber = %d;", is, runnumber);

      res = db->Query(sql);

      int nrows = res->GetRowCount();

      int nfields = res->GetFieldCount();
      //cout << nrows << " " << nfields << endl;   
      for (int i = 0; i < nrows; i++) {
        row = res->Next();
        for (int j = 0; j < nfields; j++) {
          scaledowns[is] = stoi(row->GetField(j));
          //cout << is << endl;       
          if(is == 10 || is==17 || is==18 || is==19) cout  << is << ":" << scaledowns[is] << " ";
        }
        delete row;
      }


      delete res;
    }
  cout << endl;
  return 0;
}

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
  TH2D* hists[6];
  float ETJet[6] = {4,10,20,30,40,50};
  for(int i=0; i<6; ++i)
    {
      hists[i] = (TH2D*)thefile->Get(("h2_ecc_layer"+to_string(i)).c_str());
      hists[i]->GetXaxis()->SetTitle("Jet Eccentricity");
      hists[i]->GetYaxis()->SetTitle("Max Layer E_{T} / E_{T,jet}");
      hists[i]->GetYaxis()->SetTitleOffset(1.5);
      hists[i]->GetZaxis()->SetTitleOffset(1.3);
      hists[i]->GetZaxis()->SetTitle("Counts");
      hists[i]->Draw("COLZ");
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
      c->SaveAs(("output/rmg/chi2file"+to_string(i)+".png").c_str());
    }

  return 0;
}
