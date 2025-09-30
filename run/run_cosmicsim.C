#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
//#include <g4jets/TruthJetInput.h>
#include <caloreco/CaloTowerStatus.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/BeamBackgroundFilterAndQA.h>
#include <fstream>
#include <phool/recoConsts.h>
#include <TSystem.h>
#include <caloreco/CaloTowerCalib.h>
#include <frog/FROG.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
//#include <r24earlytreemaker/R24earlytreemaker.h>
#include <chi2checker/Chi2checker.h>
#include <trigzvtxchecker/Trigzvtxchecker.h>
#include <globalvertex/GlobalVertexReco.h>
#include <Calo_Calib.C>
#include <GlobalVertex.h>
#include <TruthJetInput.h>//#include <G4Setup_sPHENIX.C>
//#include <MbdDigitization.h>
#include <MbdReco.h>
using namespace std;

R__LOAD_LIBRARY(libchi2checker.so)
//R__LOAD_LIBRARY(libr24earlytreemaker.so)
R__LOAD_LIBRARY(libg4centrality.so)
//R__LOAD_LIBRARY(libFROG.so)
//R__LOAD_LIBRARY(libg4vertex.so)
//R__LOAD_LIBRARY(libglobalvertex.so);
R__LOAD_LIBRARY(libfun4allraw.so)
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libcalo_io.so)
R__LOAD_LIBRARY(libcalo_reco.so)
R__LOAD_LIBRARY(libg4mbd.so)
R__LOAD_LIBRARY(libmbd_io.so)
R__LOAD_LIBRARY(libmbd.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libg4jets.so)
R__LOAD_LIBRARY(libjetbase.so)
R__LOAD_LIBRARY(libjetbackground.so)
R__LOAD_LIBRARY(libtrigzvtxchecker.so)
R__LOAD_LIBRARY(libg4dst.so)
//gSystem->Load("libg4detectors.so");
//gSystem->Load("libg4detectors.so");

bool file_exists(const char* filename)
{
  std::ifstream infile(filename);
  return infile.good();
}
int run_cosmicsim(string tag = "", int nproc = 0, int debug = 0, int nevt = 0, int rn = 0, int szs = 0, int isdat = 1, int chi2check = 0, int sampletype = -1, string dir = ".", int dowf = 0)
{
  int verbosity = 0;//debug;
  string chi2filename = dir+"/"+to_string(rn)+"_chi2/events_"+tag+"_"+to_string(rn)+"_"+to_string(nproc)+"_"+to_string(nevt)+"_chi2file.root";
  
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc =  recoConsts::instance();
  rc->set_StringFlag("CDB_GLOBALTAG","MDC2");
  rc->set_uint64Flag("TIMESTAMP",28);
  
  se->Verbosity(verbosity);
  // just if we set some flags somewhere in this macro


  Fun4AllInputManager *in_1 = new Fun4AllDstInputManager("DSTin1");

  string line1;
  line1 = "./dsts/"+to_string(nproc)+"/cosmics_"+to_string(nproc)+".root";
  in_1->AddFile(line1);
  se->registerInputManager( in_1 );
  
  CDBInterface::instance()->Verbosity(0);

  Process_Calo_Calib();






  
  
  NullFilter::Config cfg_null {
    .verbosity = verbosity,
    .debug = false
  };

  StreakSidebandFilter::Config cfg_sideband {
    .verbosity = verbosity,
    .debug = false,
    .minStreakTwrEne = 0.6,
    .maxAdjacentTwrEne = 0.06,
    .minNumTwrsInStreak = 5
  };

  BeamBackgroundFilterAndQA::Config cfg_filter {
    .debug = false,
    .doQA = false,
    .doEvtAbort = false,
    .sideband = cfg_sideband
  };
  
  BeamBackgroundFilterAndQA* bgelim = new BeamBackgroundFilterAndQA(cfg_filter);
  se->registerSubsystem(bgelim);
  se->Print("NODETREE");



  RetowerCEMC *rcemc = new RetowerCEMC();
  rcemc->set_towerinfo(true);
  rcemc->Verbosity(verbosity);
  se->registerSubsystem(rcemc);
  cout << "set up retower emcal" << endl;
    
  JetReco *towerjetreco = new JetReco();
  TowerJetInput* emtji = new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,"TOWERINFO_CALIB");
  TowerJetInput* ohtji = new TowerJetInput(Jet::HCALIN_TOWERINFO,"TOWERINFO_CALIB");
  TowerJetInput* ihtji = new TowerJetInput(Jet::HCALOUT_TOWERINFO,"TOWERINFO_CALIB");
      
  towerjetreco->add_input(emtji);
  towerjetreco->add_input(ohtji);
  towerjetreco->add_input(ihtji);
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "AntiKt_Tower_HIRecoSeedsRaw_r04");
  towerjetreco->set_algo_node("ANTIKT");
  towerjetreco->set_input_node("TOWER");
  towerjetreco->Verbosity(verbosity);
  se->registerSubsystem(towerjetreco);
  string wffname = "null.root";
  Chi2checker* chi2c;
  int doall60 = 0;
  if(chi2check) chi2c = new Chi2checker(chi2filename,to_string(rn)+"_"+to_string(nproc),debug,wffname,dowf,(isdat)?true:false,doall60,1);
  if(chi2check) se->registerSubsystem(chi2c);
  se->Print("NODETREE");
  cout << "run " << nevt << endl;
  se->run(nevt);
  se->Print("NODETREE");
  cout << "ending" << endl;
  se->End();
  cout << "printing timer" << endl;
  se->PrintTimer();
  cout << "Ended server" << endl;
  return 0;

}
