void check_scaledowns(int runnumber)
{
  gROOT->ProcessLine( "gErrorIgnoreLevel = 2002;");
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
  float effevt = (((1.*(1+sd10))/(1+sd17))*nmb);

  cout << runnumber << " " << sd10 << " " << sd17 << " " << nmb << " " << njt << " " << effevt << endl;
  if(!std::isfinite(sd10) || !std::isfinite(sd17) || !std::isfinite(nmb))cout << "PANIC" << endl;
  delete row;
  delete res;
  delete db;
}
