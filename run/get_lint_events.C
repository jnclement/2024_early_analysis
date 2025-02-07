float get_effevt(int runnumber, int runevt)
{
  
  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","","");

  TSQLRow *row;
  TSQLResult *res;
  char sql[1000];


  
  sprintf(sql, "select scaledown17 from gl1_scaledown where runnumber = %d;", runnumber);
  res = db->Query(sql);
  if(!res) return 0;
  row = res->Next();
  if(!row) return 0;
  int sd17 = stoi(row->GetField(0));
  
  sprintf(sql, "select scaledown10 from gl1_scaledown where runnumber = %d;",runnumber);// and index = 10;", runnumber);
  res = db->Query(sql);
  if(!res) return 0;
  row = res->Next();
  if(!row) return 0;
  int sd10 = std::stoull(row->GetField(0));
  cout << sd10 << " " << sd17 << endl;

  float effevt = runevt*(1.+sd10)/(1+sd17);
  delete row;
  delete res;
  delete db;
  if(!std::isfinite(sd17) || !std::isfinite(sd10) || sd17 < 0 || !std::isfinite(effevt)) return 0;
  cout << runnumber << " " << effevt << endl;
  return effevt;
}

void get_lint_events(string filename)
{
  
  TFile* mbfile = TFile::Open(filename.c_str());
  gROOT->ProcessLine( "gErrorIgnoreLevel = 2002;");
  int rn;
  int nevtrun;
  int mbevt;
  int prn = 0;
  float effevt = 0;
  TTree* mbtree = (TTree*)mbfile->Get("mbtree");
  mbtree->SetBranchAddress("mbevt",&mbevt);
  mbtree->SetBranchAddress("rn",&rn);
  for(int i=0; i<mbtree->GetEntries(); ++i)
    {
      mbtree->GetEntry(i);
      cout << nevtrun << endl;
      nevtrun += mbevt;
      if(rn == prn) continue;
      prn = rn;
      effevt += get_effevt(rn,nevtrun);
      nevtrun = 0;
    }
  cout << effevt << endl;
  effevt /= ((21e-3)*(1e12));
  cout << effevt << endl;
}
