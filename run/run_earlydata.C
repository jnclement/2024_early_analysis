#include <calowaveformsim/EliminateBackground.h>
//#include <trigger/TriggerRunInfo.h>
//#include <trigger/TriggerRunInfoReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
//#include <caloreco/CaloTowerStatus.h>
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
//#include <caloreco/CaloTowerCalib.h>
#include <g4mbd/MbdDigitization.h>
#include <mbd/MbdReco.h>
#include <frog/FROG.h>
#include <ffamodules/CDBInterface.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <g4centrality/PHG4CentralityReco.h>
#include <centrality/CentralityReco.h>

#include <r24earlytreemaker/R24earlytreemaker.h>
#include <chi2checker/Chi2checker.h>
#include </sphenix/user/jocl/projects/macros/common/Calo_Calib.C>
#include <CaloTowerCalib.h>
#include <globalvertex/GlobalVertexReco.h>
#include <beambackgroundfilterandqa/BeamBackgroundFilterAndQA.h>
//#include <G4Setup_sPHENIX.C>
using namespace std;
R__LOAD_LIBRARY(libchi2checker.so)
R__LOAD_LIBRARY(libr24earlytreemaker.so)
R__LOAD_LIBRARY(libg4centrality.so)
R__LOAD_LIBRARY(libFROG.so)
R__LOAD_LIBRARY(libg4vertex.so)
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
R__LOAD_LIBRARY(libemulatortreemaker.so)
R__LOAD_LIBRARY(libbeambackgroundfilterandqa.so);
bool file_exists(const char* filename)
{
  std::ifstream infile(filename);
  return infile.good();
}
int run_earlydata(string tag = "", int nproc = 0, int debug = 0, int nevt = 0, int rn = 0, int szs = 0, int datorsim = 1, int chi2check = 0, string dir = ".")
{
  //cout << "test0" << endl;
  int verbosity = 0;
  string filename = dir+"/"+to_string(rn)+"/events_"+tag+(tag==""?"":"_");
  //filename += (szs?"yszs_":"nszs_")+to_string(szs)+"_"+
  filename += to_string(rn)+"_";
  filename += to_string(nproc)+"_";
  filename += to_string(nevt);
  string chi2filename = dir+"/"+to_string(rn)+"_chi2/events_"+tag+"_"+to_string(rn)+"_"+to_string(nproc)+"_"+to_string(nevt)+"_chi2file.root";
  filename += ".root";
  FROG *fr = new FROG();
  //cout << "test0.5" << endl;
  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4detectors.so");
  gSystem->Load("libjetbackground.so");
  gSystem->Load("libcalo_io.so");
  gSystem->Load("libg4dst.so");
  //cout << "test1" << endl;
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity( verbosity );
  // just if we set some flags somewhere in this macro
  recoConsts *rc =  recoConsts::instance();
  if(datorsim) rc->set_StringFlag("CDB_GLOBALTAG","2024p007");
  else rc->set_StringFlag("CDB_GLOBALTAG","MDC2");
  if(datorsim) rc->set_uint64Flag("TIMESTAMP",rn);
  else rc->set_uint64Flag("TIMESTAMP",15);
  /*
  ifstream list1;
  string line1;
  ifstream list2;
  string line2;
  
  string line3;
  list1.open(datorsim?("./lists/"+to_string(rn)+".list"):"lists/dst_calo_waveform.list", ifstream::in);
  if(!datorsim) list2.open("lists/dst_global.list",ifstream::in);
  if(!datorsim) list3.open("lists/dst_truth_jet.list",ifstream::in);

  if(!list1)
    {
      cout << "nolist!" << endl;
      exit(1);
    }
  //cout << list1 << endl;
  */
  Fun4AllInputManager *in_1 = new Fun4AllDstInputManager("DSTin1");
  Fun4AllInputManager *in_2 = new Fun4AllDstInputManager("DSTin2");
  Fun4AllInputManager *in_3 = new Fun4AllDstInputManager("DSTin3");
  /*
  for(int i=0; i<nproc+1; i++)
    {
      getline(list1, line1);
      //cout << line1 << endl;
      if(!datorsim && list2) getline(list2, line2);
      if(!datorsim) getline(list3, line3);
    }
  */
  ifstream list3, list2, list1;
  if(!datorsim) list3.open("lists/dst_truth_jet.list",ifstream::in);
  if(!datorsim && !list3) list3.open("lists/g4hits.list");
  if(!datorsim) list2.open("lists/dst_global.list",ifstream::in);
  string line1, line2, line3;
  if(datorsim) line1 = "./dsts/"+to_string(rn)+"/"+to_string(rn)+"_"+to_string(nproc)+".root";
  else line1 = "./dsts/"+to_string(rn)+"/calo_waveform_"+to_string(nproc)+".root";
  line2 = "./dsts/"+to_string(rn)+"/global_"+to_string(nproc)+".root";
  if(list3) line3 = "./dsts/"+to_string(rn)+"/truth_jet_"+to_string(nproc)+".root";
  else line3 = "./dsts/"+to_string(rn)+"/g4hits_"+to_string(nproc)+".root";
  in_1->AddFile(line1);
  if(!datorsim && list2) in_2->AddFile(line2);
  if(!datorsim) in_3->AddFile(line3);
  se->registerInputManager( in_1 );
  if(!datorsim) se->registerInputManager( in_2 );
  if(!datorsim) se->registerInputManager( in_3 );
  // this points to the global tag in the CDB
  //rc->set_StringFlag("CDB_GLOBALTAG","");//"ProdA_2023");                                     
  // The calibrations have a validity range set by the beam clock which is not read out of the prdfs as of now
  //Process_Calo_Calib();


  std::cout << "status setters" << std::endl;
  if(!datorsim)
    {
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
      
      ////////////////////
      // Calibrate towers
      
      std::cout << "Calibrating EMCal" << std::endl;
      CaloTowerCalib *calibEMC = new CaloTowerCalib("CEMCCALIB");
      calibEMC->set_detector_type(CaloTowerDefs::CEMC);
      se->registerSubsystem(calibEMC);
      
      std::cout << "Calibrating OHcal" << std::endl;
      CaloTowerCalib *calibOHCal = new CaloTowerCalib("HCALOUT");
      calibOHCal->set_detector_type(CaloTowerDefs::HCALOUT);
      se->registerSubsystem(calibOHCal);
      
      std::cout << "Calibrating IHcal" << std::endl;
      CaloTowerCalib *calibIHCal = new CaloTowerCalib("HCALIN");
      calibIHCal->set_detector_type(CaloTowerDefs::HCALIN);
      se->registerSubsystem(calibIHCal);
      
    }
  /*
  std::cout << "Calibrating ZDC" << std::endl;
  CaloTowerCalib *calibZDC = new CaloTowerCalib("ZDC");
  calibZDC->set_detector_type(CaloTowerDefs::ZDC);
  se->registerSubsystem(calibZDC);
  */
  //////////////////
  // Clusters
  /*
  std::cout << "Building clusters" << std::endl;
  RawClusterBuilderTemplate *ClusterBuilder = new RawClusterBuilderTemplate("EmcRawClusterBuilderTemplate");
  ClusterBuilder->Detector("CEMC");
  ClusterBuilder->set_threshold_energy(0.030);  // for when using basic calibration
  std::string emc_prof = getenv("CALIBRATIONROOT");
  emc_prof += "/EmcProfile/CEMCprof_Thresh30MeV.root";
  ClusterBuilder->LoadProfile(emc_prof);
  ClusterBuilder->set_UseTowerInfo(1);  // to use towerinfo objects rather than old RawTower
  se->registerSubsystem(ClusterBuilder);
  */

  CDBInterface::instance()->Verbosity(0);
  int cont = 0;
  //cout << "test1.5" << endl;  

  //EliminateBackground* bgelim = new EliminateBackground("bgelim");
  //se->registerSubsystem(bgelim);
  auto mbddigi = new MbdDigitization();
  auto mbdreco = new MbdReco();
  GlobalVertexReco* gblvertex = new GlobalVertexReco();
  if (!datorsim)
    {
      
      mbddigi->Verbosity(verbosity);
      se->registerSubsystem(mbddigi);
      mbdreco->Verbosity(verbosity);
      se->registerSubsystem(mbdreco);
      
      gblvertex->Verbosity(verbosity);
      se->registerSubsystem(gblvertex);
    }

  se->Print("NODETREE");



  //TriggerRunInfoReco* tana = new TriggerRunInfoReco("tana");
  //se->registerSubsystem(tana);
  
  RetowerCEMC *rcemc = new RetowerCEMC();
  rcemc->set_towerinfo(true);
  rcemc->Verbosity(verbosity);
  se->registerSubsystem(rcemc);
  cout << "set up retower emcal" << endl;
  JetReco *truthjetreco = new JetReco();
  if(!datorsim)
    {
      TruthJetInput *tji = new TruthJetInput(Jet::PARTICLE);
      tji->add_embedding_flag(0);  // changes depending on signal vs. embedded
      truthjetreco->add_input(tji);
      truthjetreco->add_algo(new FastJetAlgo(Jet::ANTIKT, 0.4), "AntiKt_Truth_r04");
      truthjetreco->set_algo_node("ANTIKT");
      truthjetreco->set_input_node("G4TruthInfo");
      se->registerSubsystem(truthjetreco);
    }

  JetReco *towerjetreco = new JetReco();
  //towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWER));
  towerjetreco->add_input(new TowerJetInput(Jet::CEMC_TOWERINFO_RETOWER,"TOWERINFO_CALIB"));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALIN_TOWERINFO));
  towerjetreco->add_input(new TowerJetInput(Jet::HCALOUT_TOWERINFO));
  towerjetreco->add_algo(new FastJetAlgoSub(Jet::ANTIKT, 0.4), "AntiKt_Tower_HIRecoSeedsRaw_r04");
  towerjetreco->set_algo_node("ANTIKT");
  towerjetreco->set_input_node("TOWER");
  //towerjetreco->Verbosity(verbosity);
  se->registerSubsystem(towerjetreco);
  cout << "set up jetreco" << endl;

  BeamBackgroundFilterAndQA* bfq = new BeamBackgroundFilterAndQA();
  se->registerSubsystem(bfq);

  Chi2checker* chi2c;

  if(chi2check) chi2c = new Chi2checker(chi2filename,to_string(rn)+"_"+to_string(nproc),debug);
  if(chi2check) se->registerSubsystem(chi2c);
  cout << "set up chi2check" << endl;
  //cout << "test2" << endl;
  R24earlytreemaker *tt = new R24earlytreemaker(filename, debug, datorsim, 1);
  //cout << "test3" << endl;
  //se->registerSubsystem( tt );
  cout << "test4" << endl;
  se->Print("NODETREE");
  cout << "run " << nevt << endl;
  //se->skip(1000);
  se->run(nevt);
  cout << "ran " << nevt << endl;
  cout << "Ran all events" << endl;
  se->Print("NODETREE");
  se->End();
  se->PrintTimer();
  cout << "Ended server" << endl;
  delete se;
  cout << "Deleted server" << endl;
  cout << "wrote " << filename << endl;
  gSystem->Exit(0);
  cout << "Exited gSystem" << endl;
  return 0;

}
