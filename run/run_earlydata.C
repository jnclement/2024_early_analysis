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
#include <r24earlytreemaker/R24earlytreemaker.h>
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
R__LOAD_LIBRARY(libr24earlytreemaker.so)
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
int run_earlydata(string tag = "", int nproc = 0, int debug = 0, int nevt = 0, int rn = 0, int szs = 0, int datorsim = 1, int chi2check = 0, int sampletype = -1, string dir = ".")
{
int verbosity = 0;//debug;
  string filename = dir+"/"+to_string(datorsim?rn:nproc)+"/events_"+tag+(tag==""?"":"_");
  filename += to_string(datorsim?rn:nproc)+"_";
  filename += to_string(nproc)+"_";
  filename += to_string(nevt);
  string chi2filename = dir+"/"+to_string(rn)+"_chi2/events_"+tag+"_"+to_string(rn)+"_"+to_string(nproc)+"_"+to_string(nevt)+"_chi2file.root";

  string trigzvtxfilename = dir+"/"+to_string(rn)+"_chi2/events_"+tag+"_"+to_string(rn)+"_"+to_string(nproc)+"_"+to_string(nevt)+"_mbtree.root";
  filename += ".root";
  
  if(debug > 1) cout << "test1" << endl;
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc =  recoConsts::instance();
  if(datorsim) rc->set_StringFlag("CDB_GLOBALTAG","ProdA_2024");
  else rc->set_StringFlag("CDB_GLOBALTAG","MDC2");
  if(datorsim) rc->set_uint64Flag("TIMESTAMP",rn);
  else rc->set_uint64Flag("TIMESTAMP",21);
  
  se->Verbosity(verbosity);
  // just if we set some flags somewhere in this macro


  Trigzvtxchecker* tz;
  if(datorsim) tz = new Trigzvtxchecker(trigzvtxfilename, rn, nproc, debug, "tzvtx");
  if(datorsim) se->registerSubsystem(tz);

  Fun4AllInputManager *in_1 = new Fun4AllDstInputManager("DSTin1");
  //Fun4AllInputManager *in_2 = new Fun4AllDstInputManager("DSTin2");
  Fun4AllInputManager *in_3 = new Fun4AllDstInputManager("DSTin3");
  Fun4AllInputManager *in_4 = new Fun4AllDstInputManager("DSTin4");
  cout << "get filenames" << endl;
  ifstream list3, list2, list1;
  //if(!datorsim) list3.open("lists/dst_truth_jet.list",ifstream::in);
  //if(!datorsim) list3.open("lists/g4hits.list");
  //if(!datorsim) list2.open("lists/dst_global.list",ifstream::in);
  string line1, line2, line3, line4;
  if(datorsim) line1 = "./dsts/"+to_string(rn)+"/"+to_string(rn)+"_"+to_string(nproc)+".root";
  else line1 = "./dsts/"+to_string(nproc)+"/calo_cluster_"+to_string(nproc)+".root";
  //line2 = "./dsts/"+to_string(nproc)+"/global_"+to_string(nproc)+".root";
  if(!datorsim) line3 = "./dsts/"+to_string(nproc)+"/mbd_epd_"+to_string(nproc)+".root";
  else line3 = "./dsts/"+to_string(rn)+"/"+to_string(rn)+"_"+to_string(nproc)+"_jetcalo.root";
  //line4 = "./dsts/"+to_string(nproc)+"/g4hits_"+to_string(nproc)+".root";
  line4 = "./dsts/"+to_string(nproc)+"/truth_jet_"+to_string(nproc)+".root";
  in_1->AddFile(line1);
  //if(!datorsim) in_2->AddFile(line2);
  in_3->AddFile(line3);
  if(!datorsim) in_4->AddFile(line4);
  cout << "register managers" << endl;
  se->registerInputManager( in_1 );
  
  //if(!datorsim) se->registerInputManager( in_2 );
  se->registerInputManager( in_3 );
  if(!datorsim) se->registerInputManager(in_4);

  std::cout << "status setters" << std::endl;

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
  //auto mbddigi = new MbdDigitization();
  auto mbdreco = new MbdReco();
  GlobalVertexReco* gblvertex = new GlobalVertexReco();
  
  if (!datorsim)
    {
      //      mbddigi->Verbosity(verbosity);
      //se->registerSubsystem(mbddigi);
      mbdreco->Verbosity(verbosity);
      se->registerSubsystem(mbdreco);
      
      gblvertex->Verbosity(verbosity);
      se->registerSubsystem(gblvertex);
    }
  
  se->Print("NODETREE");



  //TriggerRunInfoReco* tana = new TriggerRunInfoReco("tana");
  //se->registerSubsystem(tana);
  RetowerCEMC *rcemc = new RetowerCEMC();
  if(!datorsim)
    {
      rcemc->set_towerinfo(true);
      rcemc->Verbosity(verbosity);
      se->registerSubsystem(rcemc);
      cout << "set up retower emcal" << endl;
    }
  /*
  JetReco *truthjetreco = new JetReco();
  
  if(!datorsim)
    {
      TruthJetInput *tji = new TruthJetInput(Jet::PARTICLE);
      tji->add_embedding_flag(0);  // changes depending on signal vs. embedded
      truthjetreco->add_input(tji);
      truthjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
      truthjetreco->set_algo_node("ANTIKT");
      truthjetreco->set_input_node("TRUTH");
      //se->registerSubsystem(truthjetreco);
    }
  */
  
  JetReco *towerjetreco = new JetReco();
  TowerJetInput* emtji = new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,"TOWERINFO_CALIB");
  TowerJetInput* ohtji = new TowerJetInput(Jet::HCALIN_TOWERINFO);
  TowerJetInput* ihtji = new TowerJetInput(Jet::HCALOUT_TOWERINFO);
  if(!datorsim)
    {
      //towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWER));
      emtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
      ohtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
      ihtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
      
      towerjetreco->add_input(emtji);
      towerjetreco->add_input(ohtji);
      towerjetreco->add_input(ihtji);
      towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "AntiKt_Tower_HIRecoSeedsRaw_r04");
      towerjetreco->set_algo_node("ANTIKT");
      towerjetreco->set_input_node("TOWER");
      towerjetreco->Verbosity(verbosity);
      se->registerSubsystem(towerjetreco);
      
      cout << "set up jetreco" << endl;
    }
  

  
  Chi2checker* chi2c;
  if(chi2check) chi2c = new Chi2checker(chi2filename,to_string(rn)+"_"+to_string(nproc),debug);
  if(chi2check) se->registerSubsystem(chi2c);
  cout << "set up chi2check" << endl;
  
  cout << "test2" << endl;
  R24earlytreemaker *tt = new R24earlytreemaker(filename, debug, datorsim, 0, sampletype);
  cout << "test3" << endl;
  if(!datorsim) se->registerSubsystem( tt );
  
  cout << "test4" << endl;
  se->Print("NODETREE");
  cout << "run " << nevt << endl;
  se->run(nevt);
  cout << "ran " << nevt << endl;
  cout << "Ran all events" << endl;
  se->Print("NODETREE");
  cout << "ending" << endl;
  se->End();
  cout << "printing timer" << endl;
  se->PrintTimer();
  cout << "Ended server" << endl;
  delete se;
  cout << "Deleted server" << endl;
  //cout << "wrote " << filename << endl;
  //gSystem->Exit(0);
  //cout << "Exited gSystem" << endl;
  return 0;

}
