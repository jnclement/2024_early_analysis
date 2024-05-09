#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
#include <caloreco/CaloTowerStatus.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
//#include <jetbase/FastJetAlgo.h>
//#include <jetbase/JetReco.h>
//#include <jetbase/TowerJetInput.h>
//#include <g4jets/TruthJetInput.h>
#include <fstream>
#include <phool/recoConsts.h>
#include <TSystem.h>
#include <caloreco/CaloTowerCalib.h>
#include <g4mbd/MbdDigitization.h>
#include <mbd/MbdReco.h>
#include <frog/FROG.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <g4centrality/PHG4CentralityReco.h>
#include <centrality/CentralityReco.h>
#include <r24earlytreemaker/R24earlytreemaker.h>
//#include <G4Setup_sPHENIX.C>
using namespace std;
R__LOAD_LIBRARY(libr24earlytreemaker.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libFROG.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libg4vertex.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalo_io.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4mbd.so)
R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbackground.so)
bool file_exists(const char* filename)
{
  std::ifstream infile(filename);
  return infile.good();
}
int run_earlydata(string tag = "", int nproc = 0, int debug = 0, int nevt = 0, int rn = 0, int szs = 0)
{
  cout << "test0" << endl;
  int verbosity = 0;
  string filename = "output/evt/events_"+tag+(tag==""?"":"_");
  filename += (szs?"yszs_":"nszs_")+to_string(szs)+"_"+to_string(rn)+"_";
  filename += to_string(nproc);
  filename += ".root";
  FROG *fr = new FROG();
  cout << "test0.5" << endl;
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libjetbackground.so");
  gSystem->Load("libcalo_io.so");
  gSystem->Load("libg4dst.so");
  cout << "test1" << endl;
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity( verbosity );
  // just if we set some flags somewhere in this macro
  recoConsts *rc =  recoConsts::instance();
  rc->set_uint64Flag("TIMESTAMP",rn);
  ifstream list1;
  string line1;
  list1.open(("./"+to_string(rn)+(szs?"_ysZS":"_noZS")+".list"), ifstream::in);
  if(!list1) exit(1);
  for(int i=0; i<nproc+1; i++)
    {
      getline(list1, line1);
    }
  Fun4AllInputManager *in_1 = new Fun4AllDstInputManager("DSTin1");
  in_1->AddFile(line1);
  se->registerInputManager( in_1 );
  // this points to the global tag in the CDB
  rc->set_StringFlag("CDB_GLOBALTAG","2023p009");//"ProdA_2023");                                     
// The calibrations have a validity range set by the beam clock which is not read out of the prdfs as of now
  CDBInterface::instance()->Verbosity(0);
  int cont = 0;
  /*
  DansSpecialVertex *dsv;
  if(!datormc)
    {
      dsv = new DansSpecialVertex("DansSpecialVertex", "dump.root");
      dsv->SetRunNumber(runnumber);
      dsv->Verbosity(0);
      se->registerSubsystem(dsv);
    }
  */
  //Fun4AllInputManager *intrue2 = new Fun4AllRunNodeInputManager("DST_GEO");

  //CDBInterface *cdb = CDBInterface::instance();
  //std::string geoLocation = cdb->getUrl("calo_geo");
  //intrue2->AddFile(geoLocation);
  //se->registerInputManager(intrue2);
  /*
  CaloTowerStatus *statusEMC = new CaloTowerStatus("CEMCSTATUS");
  statusEMC->set_detector_type(CaloTowerDefs::CEMC);
  statusEMC->set_time_cut(1);
  se->registerSubsystem(statusEMC);

  CaloTowerStatus *statusHCalIn = new CaloTowerStatus("HCALINSTATUS");
  statusHCalIn->set_detector_type(CaloTowerDefs::HCALIN);
  statusHCalIn->set_time_cut(2);
  se->registerSubsystem(statusHCalIn);

  CaloTowerStatus *statusHCALOUT = new CaloTowerStatus("HCALOUTSTATUS");
  statusHCALOUT->set_detector_type(CaloTowerDefs::HCALOUT);
  statusHCALOUT->set_time_cut(2);
  se->registerSubsystem(statusHCALOUT);
  */
  CaloTowerCalib* ctcem = new CaloTowerCalib("EMCALIB");
  ctcem->set_detector_type(CaloTowerDefs::CEMC);
  CaloTowerCalib* ctcih = new CaloTowerCalib("IHCALIB");
  ctcih->set_detector_type(CaloTowerDefs::HCALIN);
  CaloTowerCalib* ctcoh = new CaloTowerCalib("OHCALIB");
  ctcoh->set_detector_type(CaloTowerDefs::HCALOUT);
  se->registerSubsystem(ctcem);
  se->registerSubsystem(ctcih);
  se->registerSubsystem(ctcoh);
  RetowerCEMC *rcemc = new RetowerCEMC();
  rcemc->set_towerinfo(true);
  se->registerSubsystem(rcemc);
  
  JetReco *towerjetreco = new JetReco();
  //towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWER));
  towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO));
  towerjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Tower_HIRecoSeedsRaw_r02");
  towerjetreco->set_algo_node("ANTIKT");
  towerjetreco->set_input_node("TOWER");
  //towerjetreco->Verbosity(verbosity);
  se->registerSubsystem(towerjetreco);
  cout << "test2" << endl;
  R24earlytreemaker *tt = new R24earlytreemaker(filename, debug);
  cout << "test3" << endl;
  se->registerSubsystem( tt );
  cout << "test4" << endl;
  se->Print("NODETREE");
  se->run(nevt);
  cout << "Ran all events" << endl;
  se->End();
  cout << "Ended server" << endl;
  delete se;
  cout << "Deleted server" << endl;
  gSystem->Exit(0);
  cout << "Exited gSystem" << endl;
  return 0;

}
