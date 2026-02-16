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
//#include <jetbackground/BeamBackgroundFilterAndQA.h>
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
//#include </sphenix/user/jocl/projects/run2024_earlydata/run/Calo_Calib.C>
#include <Calo_Calib.C>
#include <GlobalVertex.h>
#include <TruthJetInput.h>//#include <G4Setup_sPHENIX.C>
//#include <MbdDigitization.h>
#include <MbdReco.h>
#include <jetbackground/TimingCut.h>
#include <jetbase/JetCalib.h>
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
int run_earlydata(string tag = "", int nproc = 0, int debug = 0, int nevt = 0, int rn = 0, int szs = 0, int isdat = 1, int chi2check = 0, int sampletype = -1, string dir = ".", int dowf = 0)
{
  int verbosity = debug;
  string filename = dir+"/"+to_string(isdat?rn:nproc)+"/events_"+tag+(tag==""?"":"_");
  cout << "test" << endl;
  filename += to_string(isdat?rn:nproc)+"_";
  filename += to_string(nproc)+"_";
  filename += to_string(nevt);
  string chi2filename = dir+"/"+to_string(rn)+"_chi2/events_"+tag+"_"+to_string(rn)+"_"+to_string(nproc)+"_"+to_string(nevt)+"_chi2file.root";
  cout << "another test" << endl;
  string wffname = dir+"/"+to_string(rn)+"_chi2/events_"+tag+"_"+to_string(rn)+"_"+to_string(nproc)+"_"+to_string(nevt)+"_waveforms.root";

  string trigzvtxfilename = dir+"/"+to_string(rn)+"_chi2/events_"+tag+"_"+to_string(rn)+"_"+to_string(nproc)+"_"+to_string(nevt)+"_mbtree.root";
  filename += ".root";
  
  cout << "test1" << endl;
  Fun4AllServer *se = Fun4AllServer::instance();
  recoConsts *rc =  recoConsts::instance();
  if(isdat) rc->set_StringFlag("CDB_GLOBALTAG","ProdA_2024");
  else rc->set_StringFlag("CDB_GLOBALTAG","MDC2");
  if(isdat) rc->set_uint64Flag("TIMESTAMP",rn);
  else rc->set_uint64Flag("TIMESTAMP",28);
  
  se->Verbosity(verbosity);
  // just if we set some flags somewhere in this macro




  Fun4AllDstInputManager *in_1 = new Fun4AllDstInputManager("DSTin1");
  Fun4AllDstInputManager *in_2 = new Fun4AllDstInputManager("DSTin2");
  Fun4AllDstInputManager *in_3 = new Fun4AllDstInputManager("DSTin3");
  Fun4AllDstInputManager *in_4 = new Fun4AllDstInputManager("DSTin4");
  in_1->Verbosity(debug);
  in_2->Verbosity(debug);
  in_3->Verbosity(debug);
  in_4->Verbosity(debug);
  cout << "get filenames" << endl;
  ifstream list3, list2, list1;
  //if(!isdat) list3.open("lists/dst_truth_jet.list",ifstream::in);
  //if(!isdat) list3.open("lists/g4hits.list");
  //if(!isdat) list2.open("lists/dst_global.list",ifstream::in);
  string line1, line2, line3, line4;
  if(isdat) line1 = "./dsts/"+to_string(rn)+"/"+to_string(rn)+"_"+to_string(nproc)+".root";
  else line1 = "./calo_cluster_temp.list";
  if(!isdat) line3 = "./mbd_epd_temp.list";
  if(!isdat) line2 = "./g4hits_temp.list";
  line4 = "./truth_jet_temp.list";
  if(isdat) in_1->AddListFile("thelist.list");
  else in_1->AddListFile(line1);
  if(!isdat) in_2->AddListFile(line2);
  if(!isdat) in_3->AddListFile(line3);
  if(!isdat) in_4->AddListFile(line4);
  cout << "register managers" << endl;
  se->registerInputManager( in_1 );
  
  
  if(!isdat)
    {
      cout << "registering special sim input managers" << endl;
      se->registerInputManager( in_2 );
      se->registerInputManager( in_3 );
      se->registerInputManager(in_4);
    }

  std::cout << "status setters" << std::endl;


  Fun4AllDstInputManager* inseb[18];

  if(dowf)
    {
      for(int i=0; i<18; ++i)
	{
	  inseb[i] = new Fun4AllDstInputManager("dstinseb"+to_string(i));
	  inseb[i]->AddFile((i<10?"./dsts/seb0":"./dsts/seb")+to_string(i)+".root");
	  se->registerInputManager(inseb[i]);
	}
    }

  
  CDBInterface::instance()->Verbosity(1);

  Process_Calo_Calib();

  //auto mbddigi = new MbdDigitization();
  auto mbdreco = new MbdReco();
  GlobalVertexReco* gblvertex = new GlobalVertexReco();
  //if (!isdat)
  //      mbddigi->Verbosity(verbosity);
  //se->registerSubsystem(mbddigi);
  mbdreco->Verbosity(verbosity);
  se->registerSubsystem(mbdreco);
  
  //Trigzvtxchecker* tz;
  //if(isdat) tz = new Trigzvtxchecker(trigzvtxfilename, rn, nproc, debug, "tzvtx");
  //if(isdat) se->registerSubsystem(tz);
  
  gblvertex->Verbosity(verbosity);
  se->registerSubsystem(gblvertex);
  
  se->Print("NODETREE");



  //TriggerRunInfoReco* tana = new TriggerRunInfoReco("tana");
  //se->registerSubsystem(tana);
  RetowerCEMC *rcemc = new RetowerCEMC();
  //if(!isdat)
    {
      rcemc->set_towerinfo(true);
      rcemc->Verbosity(verbosity);
      se->registerSubsystem(rcemc);
      cout << "set up retower emcal" << endl;
    }
  /*
  JetReco *truthjetreco = new JetReco();
  
  if(!isdat)
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
  TowerJetInput* ohtji = new TowerJetInput(Jet::HCALIN_TOWERINFO,"TOWERINFO_CALIB");
  TowerJetInput* ihtji = new TowerJetInput(Jet::HCALOUT_TOWERINFO,"TOWERINFO_CALIB");
  //if(!isdat)
  //{
      //towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWER));
  emtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  ohtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  ihtji->set_GlobalVertexType(GlobalVertex::VTXTYPE::MBD);
  
  towerjetreco->add_input(emtji);
  towerjetreco->add_input(ohtji);
  towerjetreco->add_input(ihtji);
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.2), "AntiKt_unsubtracted_r02");
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.3), "AntiKt_unsubtracted_r03");
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "AntiKt_unsubtracted_r04");
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.6), "AntiKt_unsubtracted_r06");
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.8), "AntiKt_unsubtracted_r08");
  towerjetreco->set_algo_node("ANTIKT");
  towerjetreco->set_input_node("TOWER");
  towerjetreco->Verbosity(verbosity);
  se->registerSubsystem(towerjetreco);

  JetCalib *jetCalib02 = new JetCalib("JetCalib02");
  jetCalib02->set_InputNode("AntiKt_unsubtracted_r02");
  jetCalib02->set_OutputNode("AntiKt_unsubtracted_r02_calib");
  jetCalib02->set_JetRadius(0.2);
  jetCalib02->set_ZvrtxNode("GlobalVertexMap");
  jetCalib02->set_ApplyZvrtxDependentCalib(true);
  jetCalib02->set_ApplyEtaDependentCalib(true);
  se->registerSubsystem(jetCalib02);


  JetCalib *jetCalib03 = new JetCalib("JetCalib03");
  jetCalib03->set_InputNode("AntiKt_unsubtracted_r03");
  jetCalib03->set_OutputNode("AntiKt_unsubtracted_r03_calib");
  jetCalib03->set_JetRadius(0.3);
  jetCalib03->set_ZvrtxNode("GlobalVertexMap");
  jetCalib03->set_ApplyZvrtxDependentCalib(true);
  jetCalib03->set_ApplyEtaDependentCalib(true);
  se->registerSubsystem(jetCalib03);


  JetCalib *jetCalib04 = new JetCalib("JetCalib04");
  jetCalib04->set_InputNode("AntiKt_unsubtracted_r04");
  jetCalib04->set_OutputNode("AntiKt_unsubtracted_r04_calib");
  jetCalib04->set_JetRadius(0.4);
  jetCalib04->set_ZvrtxNode("GlobalVertexMap");
  jetCalib04->set_ApplyZvrtxDependentCalib(true);
  jetCalib04->set_ApplyEtaDependentCalib(true);
  se->registerSubsystem(jetCalib04);


  JetCalib *jetCalib06 = new JetCalib("JetCalib06");
  jetCalib06->set_InputNode("AntiKt_unsubtracted_r06");
  jetCalib06->set_OutputNode("AntiKt_unsubtracted_r06_calib");
  jetCalib06->set_JetRadius(0.6);
  jetCalib06->set_ZvrtxNode("GlobalVertexMap");
  jetCalib06->set_ApplyZvrtxDependentCalib(true);
  jetCalib06->set_ApplyEtaDependentCalib(true);
  se->registerSubsystem(jetCalib06);


  JetCalib *jetCalib08 = new JetCalib("JetCalib08");
  jetCalib08->set_InputNode("AntiKt_unsubtracted_r08");
  jetCalib08->set_OutputNode("AntiKt_unsubtracted_r08_calib");
  jetCalib08->set_JetRadius(0.8);
  jetCalib08->set_ZvrtxNode("GlobalVertexMap");
  jetCalib08->set_ApplyZvrtxDependentCalib(true);
  jetCalib08->set_ApplyEtaDependentCalib(true);
  se->registerSubsystem(jetCalib08);

  cout << "set up jetreco" << endl;
      //}
  
  TimingCut* tc = new TimingCut("AntiKt_unsubtracted_r04");
  se->registerSubsystem(tc);
  
  Chi2checker* chi2c;
  int doall60 = 0;
  if(chi2check) chi2c = new Chi2checker(chi2filename,to_string(rn)+"_"+to_string(nproc),debug,wffname,dowf,(isdat)?true:false,doall60);
  if(chi2check) se->registerSubsystem(chi2c);
  cout << "set up chi2check" << endl;
  
  cout << "test2" << endl;
  //R24earlytreemaker *tt = new R24earlytreemaker(filename, debug, isdat, 0, sampletype);
  cout << "test3" << endl;
  //if(!isdat) se->registerSubsystem( tt );
  
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
  //delete se;
  cout << "Deleted server" << endl;
  //cout << "wrote " << filename << endl;
  //gSystem->Exit(0);
  //cout << "Exited gSystem" << endl;
  return 0;

}
