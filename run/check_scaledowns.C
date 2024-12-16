float get_effevt(int runnumber, int nmb)
{
  
  TSQLServer *db = TSQLServer::Connect("pgsql://sphnxdaqdbreplica:5432/daq","","");

  TSQLRow *row;
  TSQLResult *res;
  char sql[1000];


  
  sprintf(sql, "select scaledown17 from gl1_scaledown where runnumber = %d;", runnumber);
  res = db->Query(sql);
  row = res->Next();
  int sd17 = stoi(row->GetField(0));

  sprintf(sql, "select scaledown10 from gl1_scaledown where runnumber = %d;", runnumber);
  res = db->Query(sql);
  row = res->Next();
  int sd10 = stoi(row->GetField(0));
  /*
  sprintf(sql, "select scaled from gl1_scalers where runnumber = %d and index = 10;", runnumber);
  res = db->Query(sql);
  row = res->Next();
  int nmb = stoi(row->GetField(0));

  sprintf(sql, "select scaled from gl1_scalers where runnumber = %d and index = 17;", runnumber);
  res = db->Query(sql);
  row = res->Next();
  int njt;
  if(row) njt = stoi(row->GetField(0));
  else njt = 0;


  cout << runnumber << " " << sd10 << " " << sd17 << " " << nmb << " " << njt << " " << effevt << endl;
  */
  float effevt = (((1.*(1+sd10))/(1+sd17))*nmb);
  delete row;
  delete res;
  delete db;
  if(!std::isfinite(sd10) || !std::isfinite(sd17) || !std::isfinite(nmb) || sd17 < 0 || !std::isfinite(effevt)) return 0;
  return effevt;
}

void check_scaledowns(string filename)
{
  
  TFile* mbfile = TFile::Open(filename.c_str());
  gROOT->ProcessLine( "gErrorIgnoreLevel = 2002;");
  int nmb, rn;
  float effevt = 0;
  TTree* mbtree = (TTree*)mbfile->Get("mbtree");
  mbtree->SetBranchAddress("mbevt",&nmb);
  mbtree->SetBranchAddress("rn",&rn);
  for(int i=0; i<mbtree->GetEntries(); ++i)
    {
      mbtree->GetEntry(i);
      effevt += get_effevt(rn, nmb);
    }
  cout << effevt << endl;
}
