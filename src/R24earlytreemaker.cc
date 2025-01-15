#include "R24earlytreemaker.h"
#include <ffaobjects/EventHeaderv1.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainerSimv1.h>
#include <calobase/TowerInfoContainerSimv2.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <mbd/MbdPmtSimContainerV1.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <ffaobjects/EventHeaderv1.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <HepMC/GenEvent.h>
#include <mbd/MbdPmtHit.h>
#include <jetbackground/TowerBackgroundv1.h>
#include <cmath>
#include <mbd/MbdOut.h>
#include <TLorentzVector.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
#include <g4centrality/PHG4CentralityReco.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jet.h>
#include <calobase/RawTowerv1.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos
#include <TLorentzVector.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <iostream>
#include <phool/recoConsts.h>
#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/MinimumBiasInfov1.h>
#include <ffarawobjects/Gl1Packetv2.h>
#include <jetbase/JetMapv1.h>
#include <jetbase/JetMap.h>
using namespace std;

static const float radius_EM = 93.5;
static const float minz_EM = -130.23;
static const float maxz_EM = 130.23;

static const float radius_IH = 127.503;
static const float minz_IH = -170.299;
static const float maxz_IH = 170.299;

static const float radius_OH = 225.87;
static const float minz_OH = -301.683;
static const float maxz_OH = 301.683;

float get_emcal_mineta_zcorrected(float zvertex) {
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zertex);
  float minlimit = emcal_mineta;
  if (ihcal_mineta > minlimit) minlimit = ihcal_mineta;
  if (ohcal_mineta > minlimit) minlimit = ohcal_mineta;
  float maxlimit = emcal_maxeta;
  if (ihcal_maxeta < maxlimit) maxlimit = ihcal_maxeta;
  if (ohcal_maxeta < maxlimit) maxlimit = ohcal_maxeta;
  minlimit += jet_radius;
  maxlimit -= jet_radius;
  return jet_eta < minlimit || jet_eta > maxlimit;
}

//____________________________________________________________________________..
R24earlytreemaker::R24earlytreemaker(const std::string &name, const int debug, int datorsim, int dotow):
  SubsysReco("test")//).c_str())
{
  _rc = recoConsts::instance();
  _dotow = dotow;
  _evtnum = 0;
  _evtct = 0;
  _foutname = name;
  _debug = debug;
  mbevt = 0;
  _datorsim = datorsim;
}

//____________________________________________________________________________..
R24earlytreemaker::~R24earlytreemaker()
{

}

//____________________________________________________________________________..
int R24earlytreemaker::Init(PHCompositeNode *topNode)
{

  if(_debug > 1) cout << "Begin init: " << endl;
  _f = new TFile( _foutname.c_str(), "RECREATE");
  _jett = new TTree("jett","a tree of jets");
  _tree2 = new TTree("ttree2","another persevering date tree");
  _tree2->Branch("_evtct",&_evtct,"_evtct/I");
  _tree = new TTree("ttree","a persevering date tree");
  if(_datorsim) _tree->Branch("triggervec",&triggervec,"triggervec/g");
  _tree2->Branch("mbevt",&mbevt,"mbevt/I");
  //_tree->Branch("ismb",&ismb,"ismb/I");

  //_jett->Branch("ismb",&ismb,"ismb/I");
  //_tree->Branch("caloEfrac",caloEfrac,"caloEfrac[3]/F");
  //_tree->Branch("maxTowerChi2",maxTowerChi2,"maxTowerChi2[3]/F");
  //_tree->Branch("maxTowerET",maxTowerET,"maxTowerET[3]/F");
  if(_dotow)_jett->Branch("allcomp",&allcomp,"allcomp/I");
  if(_dotow)_jett->Branch("alcet",alcet,"alcet[allcomp]/F");
  if(_dotow)_tree->Branch("emetot",&emetot,"emetot/F");
  if(_dotow)_tree->Branch("ihetot",&ihetot,"ihetot/F");
  if(_dotow)_tree->Branch("ohetot",&ohetot,"ohetot/F");
  if(_dotow) 
    {
      _tree->Branch("sectorem",&sectorem,"sectorem/I"); //Number of hit sectors in the emcal
      _tree->Branch("sectorih",&sectorih,"sectorih/I"); // IHcal etc.
      _tree->Branch("sectoroh",&sectoroh,"sectoroh/I");
      //_tree->Branch("sectoremuc",&sectoremuc,"sectoremuc/I");
      _tree->Branch("emcalen",emcalen,"emcalen[sectorem]/F"); //energy per EMCal sector
      _tree->Branch("ihcalen",ihcalen,"ihcalen[sectorih]/F"); // per IHCal sector (etc.)
      _tree->Branch("ohcalen",ohcalen,"ohcalen[sectoroh]/F");
      //_tree->Branch("emcalchi2",emcalchi2,"emcalchi2[sectorem]/F"); //energy per EMCal sector
      //_tree->Branch("ihcalchi2",ihcalchi2,"ihcalchi2[sectorih]/F"); // per IHCal sector (etc.)
      //_tree->Branch("ohcalchi2",ohcalchi2,"ohcalchi2[sectoroh]/F");
      //_tree->Branch("emcalenuc",emcalenuc,"emcalenuc[sectoremuc]/F");
      _tree->Branch("emcaletabin",emcaletabin,"emcaletabin[sectorem]/I"); //eta of EMCal sector
      _tree->Branch("ihcaletabin",ihcaletabin,"ihcaletabin[sectorih]/I");
      _tree->Branch("ohcaletabin",ohcaletabin,"ohcaletabin[sectoroh]/I");
      _tree->Branch("emcalphibin",emcalphibin,"emcalphibin[sectorem]/I"); //phi of EMCal sector
      _tree->Branch("ihcalphibin",ihcalphibin,"ihcalphibin[sectorih]/I");
      _tree->Branch("ohcalphibin",ohcalphibin,"ohcalphibin[sectoroh]/I");
    }
  if(_dotow)_tree->Branch("sectormb",&sectormb,"sectormb/I");
  if(_dotow)_tree->Branch("mbenrgy",mbenrgy,"mbenrgy[sectormb]/F"); //MBD reported value (could be charge or time)
  //_tree->Branch("emcalt",emcalt,"emcalt[sectorem]/F"); //time value of EMCal sector
  //_tree->Branch("ihcalt",ihcalt,"ihcalt[sectorih]/F");
  //_tree->Branch("ohcalt",ohcalt,"ohcalt[sectoroh]/F");
  _tree->Branch("vtx",vtx,"vtx[3]/F");
  if(_dotow) _tree->Branch("sector_rtem",&sector_rtem,"sector_rtem/I");
  _tree->Branch("l2pcEta",&_l2pcEta,"l2pcEta/F");
  _tree->Branch("njet",&njet,"njet/I");
  _tree->Branch("frcem",_frcem,"frcem[njet]");
  _tree->Branch("frcoh",_frcoh,"frcoh[njet]");
  //_jett->Branch("njet",&njet,"njet/I");
  //_jett->Branch("aceta",aceta,"aceta[njet]/F");
  //_tree->Branch("seedD",&seedD,"seedD[njet]/F");
  _tree->Branch("jet_e",jet_e,"jet_e[njet]/F");
  //_jett->Branch("jet_e",jet_e,"jet_e[njet]/F");
  //_tree->Branch("jet_r",jet_r,"jet_r[njet]/F");
  _tree->Branch("jet_et",jet_et,"jet_et[njet]/F");
  _tree->Branch("jet_ph",jet_ph,"jet_ph[njet]/F");
  //_tree->Branch("failscut",&failscut,"failscut/I");
  //_jett->Branch("jet_et",jet_et,"jet_et[njet]/F");
  //_jett->Branch("jet_ph",jet_ph,"jet_ph[njet]/F");
  
  if(_dotow)
    {
      _tree->Branch("rtemen",rtemen,"rtemen[sector_rtem]/F");
      _tree->Branch("rtemet",rtemet,"rtemet[sector_rtem]/I");
      _tree->Branch("rtemph",rtemph,"rtemph[sector_rtem]/I");
    }
  
  if(_dotow)_tree->Branch("ehjet",ehjet,"ehjet[njet]/F");
  _tree->Branch("_evtnum",&_evtnum,"_evtnum/I");
  if(_dotow)_tree->Branch("cluster_n",&_cluster_n,"cluster_n/I");
  if(_dotow)_tree->Branch("cluster_E",_cluster_E,"cluster_E[cluster_n]/F");
  if(_dotow)_tree->Branch("cluster_Ecore",_cluster_Ecore,"cluster_Ecore[cluster_n]/F");
  if(_dotow)_tree->Branch("cluster_phi",_cluster_phi,"cluster_phi[cluster_n]/F");
  if(_dotow)_tree->Branch("cluster_eta",_cluster_eta,"cluster_eta[cluster_n]/F");
  //_tree->Branch("cluster_r",_cluster_r,"cluster_r[cluster_n]/F");
  //_tree->Branch("cluster_chi2",_cluster_chi2,"cluster_chi2[cluster_n]/F");
  //_tree->Branch("cluster_template_chi2",_cluster_template_chi2,"cluster_template_chi2[cluster_n]/F");
  if(_dotow)_tree->Branch("cluster_nTower",_cluster_nTower,"cluster_nTower[cluster_n]/I");
  //if(!_datorsim) _tree->Branch("ntj",&ntj,"ntj/I");
  //if(!_datorsim) _tree->Branch("tjet_e",tjet_e,"tjet_e[ntj]/F");
  //if(!_datorsim) _tree->Branch("tjet_eta",tjet_eta,"tjet_eta[ntj]/F");
  _tree->Branch("bbfqavec",&_bbfqavec,"bbfqavec/i");

  //_tree->Branch("nLayerEm",&_nLayerEm,"nLayerEm/I");
  //_tree->Branch("nLayerOh",&_nLayerOh,"nLayerOh/I");
  _tree->Branch("n2pc",&_n2pc,"n2pc/I");
  _tree->Branch("dPhi2pcd",_dPhi2pc,"dPhi2pc[n2pc]/F");
  _tree->Branch("dEta2pcd",_dEta2pc,"dEta2pc[n2pc]/F");
  /*
  _tree->Branch("emLayerJetPhi",_emLayerJetPhi,"emLayerJetPhi[nLayerEm]/F");
  _tree->Branch("ohLayerJetPhi",_ohLayerJetPhi,"ohLayerJetPhi[nLayerOh]/F");
  _tree->Branch("emLayerJetEta",_emLayerJetEta,"emLayerJetEta[nLayerEm]/F");
  _tree->Branch("ohLayerJetEta",_ohLayerJetEta,"ohLayerJetEta[nLayerOh]/F");
  _tree->Branch("emLayerJetET",_emLayerJetET,"emLayerJetET[nLayerEm]/F");
  _tree->Branch("ohLayerJetET",_ohLayerJetET,"ohLayerJetET[nLayerOh]/F");
  */
  _tree->Branch("dPhiLayer",_dPhiLayer,"dPhiLayer[njet]/F");

  //if(!_datorsim) _tree->Branch("tjet_phi",tjet_phi,"tjet_phi[ntj]/F");
  if(_debug > 1) cout << "Init done"  << endl;
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int R24earlytreemaker::InitRun(PHCompositeNode *topNode)
{
  if(_debug > 1) cout << "Initializing!" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
//____________________________________________________________________________..
int R24earlytreemaker::process_event(PHCompositeNode *topNode)
{
  _evtnum++;
  _evtct++;
  if(_debug > 1) cout << endl << endl << endl << "Beginning event processing" << endl;
  if(_debug > 1) cout << "Event " << _evtct << "Event " << _evtnum << endl;
  float mbdq = 0;
  //reset lengths to all zero
  allcomp = 0;
  ntj = 0;
  sectorem = 0;
  sectorih = 0;
  sectoroh = 0;
  sectormb = 0;
  sectorzd = 0;
  sectoremuc = 0;
  njet = 0;
  sector_rtem = 0;
  //Get towerinfocontainer objects from nodetree
  TowerInfoContainer *towersEM = findNode::getClass<TowerInfoContainerSimv1>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  //towersEM = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC");
  //if(!towersEM) towersEM = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC");
  

  _bbfqavec = _rc->get_IntFlag("HasBeamBackground_StreakSidebandFilter") << 5;
  TowerInfoContainer *rtem = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  if(!rtem) rtem = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");

  TowerInfoContainer *towersIH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALIN");
  if(!towersIH) towersIH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *towersOH = findNode::getClass<TowerInfoContainerSimv1>(topNode, "TOWERINFO_CALIB_HCALOUT");
  if(!towersOH) towersOH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALOUT");
  if(!towersOH) towersOH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALOUT");
  TowerInfoContainer *towersEMuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_CEMC");
  TowerInfoContainer *towersIHuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALIN");
  TowerInfoContainer *towersOHuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALOUT");
  //TowerInfoContainer *towersZD = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_ZDC");
  //TowerInfoContainer *rtem = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC");
  //if(!rtem) rtem = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC");
  JetContainer *jets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_Tower_HIRecoSeedsRaw_r04");
  MbdVertexMap* mbdvtxmap = findNode::getClass<MbdVertexMapv1>(topNode, "MbdVertexMap");
  GlobalVertexMap* gvtxmap = findNode::getClass<GlobalVertexMapv1>(topNode, "GlobalVertexMap");
  Gl1Packet *gl1;
  ismb = 0;

  MbdPmtContainer *mbdtow = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if(!mbdtow) mbdtow = findNode::getClass<MbdPmtContainerV1>(topNode, "MbdPmtContainer");
  if(!mbdtow) mbdtow = findNode::getClass<MbdPmtSimContainerV1>(topNode, "MbdPmtContainer");
  RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTERINFO_CEMC");
  if(!clusters) clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_CEMC");       
  JetContainer* truthjets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_Truth_r04");

  RawTowerGeomContainer *geom[3];
  geom[0] = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC");
  geom[1] = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALIN");
  geom[2] = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALOUT");


  if(mbdtow && !_datorsim)
    {
      int northhit = 0;
      int southhit = 0;
      sectormb = 128;//mbdtow->get_npmt();
      //if(_debug) cout << "Got " << sectormb << " mbd sectors in sim." << endl;
      for(int i=0; i<sectormb; ++i)
	{
	  MbdPmtHit *mbdhit = mbdtow->get_pmt(i);
	  //if(_debug > 2) cout << "PMT " << i << " address: " << mbdhit << " charge: " << mbdhit->get_q() << endl;
	  mbenrgy[i] = mbdhit->get_q();
	  if(mbenrgy[i] > 0.4 && i < 64) northhit = 1;
	  if(mbenrgy[i] > 0.4 && i > 63) southhit = 1;
	  
	  mbdq += mbdhit->get_q();
	}
      // if(_debug) cout << "n/s: " << northhit << "/" << southhit << endl;
      if(northhit && southhit && !_datorsim)
	{
	  ismb = 1;
	  ++mbevt;
	}
      else if(!_datorsim && (!northhit || !southhit))
	{
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
    }
  else
    {
      for(int i=0; i<128; ++i)
	{
	  mbenrgy[i] = -1;
	  
	}
      //if(_debug) cout << "No MBD info!" << endl;
    }
  if(_datorsim)
    {
      gl1 = findNode::getClass<Gl1Packetv2>(topNode, "GL1Packet");
      if(!gl1)
	{
	  cout << "No trigger info!" << endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
      if(_debug > 1) cout << "Getting gl1 trigger vector from: " << gl1 << endl;
      triggervec = gl1->getScaledVector();
      if(triggervec >> 10 & 0x1)
	{
	  ++mbevt;
	  ismb = 1;
	}
      else
	{
	  ismb = 0;
	}
    }
  //int isjettrig = ((triggervec >> 16) & 1) | ((triggervec >> 17) & 1) | ((triggervec >> 18) & 1) | ((triggervec >> 19) & 1);
  //if(!ismb &! isjettrig) return Fun4AllReturnCodes::EVENT_OK;
  ntj = 0;
  if(truthjets && !_datorsim)
    {
      for(int i=0; i<truthjets->size(); ++i)
	{
	  if(i>2) break;
	  Jet* jet = truthjets->get_jet(i);
	  tjet_e[ntj] = jet->get_e();
	  if(tjet_e[ntj] < 8) continue;
	  tjet_eta[ntj] = jet->get_eta();
	  tjet_phi[ntj] = jet->get_phi();
	  ntj++;
	}
    }

  vtx[0] = 0;
  vtx[1] = 0;
  vtx[2] = NAN;
  if(_datorsim || !gvtxmap)
    {
      for(auto iter = mbdvtxmap->begin(); iter != mbdvtxmap->end(); ++iter)
	{
	  MbdVertex* mbdvtx = iter->second;
	  vtx[2] = mbdvtx->get_z();
	  break;
	}
    }
  else if(gvtxmap)
    {
      auto iter = gvtxmap->begin();
      while(iter != gvtxmap->end())
	{
	  GlobalVertex* gvtx = iter->second;
	  vtx[2] = gvtx->get_z();
	  vtx[0] = gvtx->get_x();
	  vtx[1] = gvtx->get_y();
	  iter++;
	  break;
	}
    }
  if(std::isnan(vtx[2]) || abs(vtx[2]) > 150)
    {
      if(ismb) mbevt--;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  int outfailscut = 0;

  float zvtx = vtx[2];

  //JetContainer* emjets = findNode::getClass<JetContainerv1>(topNode, "emtowjet");
  //JetContainer* ohjets = findNode::getClass<JetContainerv1>(topNode, "ohtowjet");

  //_nLayerEm = 0;
  //_nLayerOh = 0;
  _n2pc = 0;
  /*
  if(emjets)
    {
      int tocheck = emjets->size();
      for(int i=0; i<tocheck; ++i)
        {
          Jet* jet = emjets->get_jet(i);
          if(jet)
            {
              float testJetET = jet->get_e()/cosh(jet->get_eta());
              if(testJetET < 8) continue;
              _emLayerJetEta[_nLayerEm] = jet->get_eta();
              if(check_bad_jet_eta(_emLayerJetEta[_nLayerEm],zvtx,0.4)) continue;
              _emLayerJetPhi[_nLayerEm] = jet->get_phi();
              _emLayerJetET[_nLayerEm] = testJetET;
              _nLayerEm++;
            }
        }
    }

  if(ohjets)
    {
      int tocheck = ohjets->size();
      for(int i=0; i<tocheck; ++i)
        {
          Jet* jet = ohjets->get_jet(i);
          if(jet)
            {
              float testJetET = jet->get_e()/cosh(jet->get_eta());
              if(testJetET < 8) continue;
              _ohLayerJetEta[_nLayerOh] = jet->get_eta();
              if(check_bad_jet_eta(_ohLayerJetEta[_nLayerOh],zvtx,0.4)) continue;
              _ohLayerJetPhi[_nLayerOh] = jet->get_phi();
              _ohLayerJetET[_nLayerOh] = testJetET;
              _nLayerOh++;
            }
        }
    }
  */
  int nchan = 1536;
  vector<vector<float>> emTowAbove1GeV;
  vector<vector<float>> ohTowAbove1GeV;
  float maxTowET
  if(towersEM)
    {
      for(int i=0; i<nchan; ++i)
        {
          TowerInfo* tower = towersEM->get_tower_at_channel(i);
          if(!tower->get_isGood()) continue;
          int key = towersEM->encode_key(i);
          const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towersEM->getTowerEtaBin(key), towersEM->getTowerPhiBin(key));
          RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey); //encode tower geometry                                                                                              

          float radius = 93.5;
          float ihEta = tower_geom->get_eta();
          float emZ = radius/(tan(2*atan(exp(-ihEta))));
          float newz = emZ - zvtx;
          float newTheta = atan2(radius,newz);
          float towerEta = -log(tan(0.5*newTheta));
          float towerPhi = tower_geom->get_phi();
	  float towerET = tower->get_energy()/cosh(towerEta);
	  if(towerET < 1) continue;
          if(towerET > maxTowET)
            {
              maxTowET = towerET;
              _l2pcEta = towerEta;
            }
          vector<float> toPush = {towerEta, towerPhi};
          emTowAbove1GeV.push_back(toPush);
        }
    }

  if(towersOH)
    {
      for(int i=0; i<nchan; ++i)
        {
          TowerInfo* tower = towersOH->get_tower_at_channel(i);
          if(!tower->get_isGood()) continue;
          int key = towersOH->encode_key(i);
          const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, towersOH->getTowerEtaBin(key), towersOH->getTowerPhiBin(key));
          RawTowerGeom *tower_geom = geom[2]->get_tower_geometry(geomkey); //encode tower geometry                                                                                              

          float radius = tower_geom->get_center_radius();
          float newz = tower_geom->get_center_z() - zvtx;
          float newTheta = atan2(radius,newz);
          float towerEta = -log(tan(0.5*newTheta));
          float towerPhi = tower_geom->get_phi();
	  float towerET = tower->get_energy()/cosh(towerEta);
	  if(towerET < 1) continue;
          if(towerET > maxTowET)
            {
              maxTowET = towerET;
              _l2pcEta = towerEta;
            }
          vector<float> toPush = {towerEta, towerPhi};
          ohTowAbove1GeV.push_back(toPush);
        }
    }

  for(int i=0; i<emTowAbove1GeV.size(); ++i)
    {
      for(int j=0; j<ohTowAbove1GeV.size(); ++j)
        {
          _dPhi2pc[_n2pc] = emTowAbove1GeV.at(i).at(1) - ohTowAbove1GeV.at(j).at(1);
          if(_dPhi2pc[_n2pc] > M_PI) _dPhi2pc[_n2pc] -= 2*M_PI;
          if(_dPhi2pc[_n2pc] < -M_PI) _dPhi2pc[_n2pc] += 2*M_PI;
          _dEta2pc[_n2pc] = emTowAbove1GeV.at(i).at(0) - ohTowAbove1GeV.at(j).at(0);
          ++_n2pc;
        }
    }

  allcomp = 0;
  for(int i=0; i<3; ++i)
    {
      maxTowerET[i] = 0;
      maxTowerChi2[i] = 0;
      caloEfrac[i] = 0;
    }
  if(_debug > 1) cout << "Getting jets: " << endl;
  if(jets)
    {
      int tocheck = jets->size();
      if(_debug > 2) cout << "Found " << tocheck << " jets to check..." << endl;
      //cout << towersEM << endl;
      for(int i=0; i<tocheck; ++i)
	{
	  Jet *jet = jets->get_jet(i);
	  if(jet)
	    {
	      jet_r[njet] = 0.4;
	      jet_et[njet] = jet->get_eta();
	      if(check_bad_jet_eta(jet_et[njet],vtx[2],0.4)) continue;
	      jet_ph[njet] = jet->get_phi();
	      jet_e[njet] = jet->get_e()/cosh(jet_et[njet]);
	    }
	  else
	    {
	      continue;
	    }
	  _frcem[njet] = 0;
	  _frcoh[njet] = 0;
	  if(jet_e[njet] < 8) continue;
	  
	  if(_debug > 2) cout << "found a good jet!" << endl;
	  float maxeovertot = 0;
	  float hcale = 0;
	  float ihcale = 0;
	  float ohcale = 0;
	  float ecale = 0;
	  int ncomp = 0;
	  TLorentzVector emAxis;
	  TLorentzVector ohAxis;
	  for(auto comp: jet->get_comp_vec())
	    {
	      unsigned int channel = comp.second;
	      TowerInfo* tower;
	      if(comp.first == 13 || comp.first == 28 || comp.first == 25)
		{
		  if(_debug > 3) cout << "em component" << endl;
		  tower = towersEM->get_tower_at_channel(channel);
		  int key = towersEM->encode_key(channel);
		  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towersEM->getTowerEtaBin(key), towersEM->getTowerPhiBin(key));
		  RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey);
		  if(_debug > 3) cout << " got tower geom" << endl;
		  float radius = 93.5;
		  float ihEta = tower_geom->get_eta();
		  float emZ = radius/(tan(2*atan(exp(-ihEta))));
		  float newz = emZ - vtx[2];
		  float newTheta = atan2(radius,newz);
		  float towerEta = -log(tan(0.5*newTheta));
		  TLorentzVector tempEM;
		  tempEM.SetPtEtaPhiE(tower->get_energy()/cosh(towerEta),towerEta,tower_geom->get_phi(),tower->get_energy());
		  emAxis += tempEM;
		  _frcem[njet] += tower->get_energy()/cosh(towerEta);
		  if(_debug > 3) cout << "end em component" << endl;
		}
	      if(comp.first == 7 || comp.first == 27)
		{
		  if(_debug > 3) cout << "oh component" << endl;
		  tower = towersOH->get_tower_at_channel(channel);
		  if(_debug > 3) cout << "got tower " << tower << endl;
		  int key = towersOH->encode_key(channel);
		  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, towersOH->getTowerEtaBin(key), towersOH->getTowerPhiBin(key));
		  RawTowerGeom *tower_geom = geom[2]->get_tower_geometry(geomkey);
		  if(_debug > 3) cout << "got tower geom " << tower_geom << endl;
		  float radius = tower_geom->get_center_radius();
		  float newz = tower_geom->get_center_z() - vtx[2];
		  float newTheta = atan2(radius,newz);
		  float towerEta = -log(tan(0.5*newTheta));
		  TLorentzVector tempOH;
		  tempOH.SetPtEtaPhiE(tower->get_energy()/cosh(towerEta),towerEta,tower_geom->get_phi(),tower->get_energy());
		  ohAxis += tempOH;
		  _frcoh[njet] += tower->get_energy()/cosh(towerEta);
		  if(_debug > 3) cout << "end oh component" << endl;
		}
	    }
	  _frcem[njet] /= jet_e[njet];
	  _frcoh[njet] /= jet_e[njet];
	  if(_debug > 3) cout << "end comp vector" << endl;
	  /*
	  for(auto comp: jet->get_comp_vec())
	    {
	      unsigned int channel = comp.second;
	      TowerInfo* tower;
	      if(_debug > 1) cout << towersIH << " " << towersOH << " " << towersEM << endl;
	      if(comp.first == 5 || comp.first == 26)
		{
		  if(_debug > 2) cout << "ihfail" << endl;
		  tower = towersIH->get_tower_at_channel(channel);
		  if(_debug > 2) cout << "towerad: " << tower << endl;
		  hcale += tower->get_energy();
		  ihcale += tower->get_energy();
		  int key = towersIH->encode_key(channel);
		  unsigned int etabin = towersIH->getTowerEtaBin(key);
		  aceta[njet] += (etabin-11.5)*1.1/12;
		  alcet[allcomp] = (etabin-11.5)*1.1/12;
		  if(tower->get_energy() > maxTowerET[1])
		    {
		      maxTowerET[1] = tower->get_energy();
		      maxTowerChi2[1] = tower->get_chi2();
		    }
		}
	      else if(comp.first == 7 || comp.first == 27)
		{
		  if(_debug > 2) cout << "ohfail" << endl;
		  tower = towersOH->get_tower_at_channel(channel);
		  if(_debug > 2) cout << "towerad: " << tower << endl;
		  hcale += tower->get_energy();
		  ohcale += tower->get_energy();
		  int key = towersOH->encode_key(channel);
		  unsigned int etabin = towersOH->getTowerEtaBin(key);
		  aceta[njet] += (etabin-11.5)*1.1/12;
		  alcet[allcomp] = (etabin-11.5)*1.1/12;
		  if(tower->get_energy() > maxTowerET[2])
		    {
		      maxTowerET[2] = tower->get_energy();
		      maxTowerChi2[2] = tower->get_chi2();
		    }
		}
	      else if(comp.first == 13 || comp.first == 28 || comp.first == 25)
		{
		  if(_debug > 2) cout << "emfail" << endl;
		  tower = towersEM->get_tower_at_channel(channel);
		  if(_debug > 2) cout << "towerad: " << tower << endl;
		  ecale += tower->get_energy();
		  int key = towersEM->encode_key(channel);
		  unsigned int etabin = towersEM->getTowerEtaBin(key);
		  aceta[njet] += (etabin-47.5)*1.1/48;
		  alcet[allcomp] = (etabin-47.5)*1.1/48;
		  if(tower->get_energy() > maxTowerET[0])
		    {
		      maxTowerET[0] = tower->get_energy();
		      maxTowerChi2[0] = tower->get_chi2();
		    }
		  if(_debug > 1) cout << "comp etabin: " << alcet[allcomp] << endl;
		}
	      else
		{
		  if(_debug > 1) cout << "No good detector for getting jet components" << endl;
		  continue;
		}
	      if(_debug > 1) cout << tower << endl;
	      float eval = tower->get_energy();
	      if(eval > maxeovertot) maxeovertot = eval;
	      allcomp++;
	      ncomp++;
	    }
	  */
	  //caloEfrac[0] = ecale/jet_e[njet];
	  //caloEfrac[1] = ihcale/jet_e[njet];
	  //caloEfrac[2] = ohcale/jet_e[njet];
	  //aceta[njet] /= ncomp;
	  //aceta[njet] /= jet_e[njet];
	  if(_debug > 2) cout << "Now filling some jet adjacent numbers: " << endl;
	  if(_debug && jet_e[njet] - ecale - hcale > 0.01) cout << jet_e[njet] << " " << ecale << " " << hcale << endl;
	  //ehjet[njet] = ecale/hcale;
	  //maxeovertot /= jet_e[njet];
	  //if(_debug) cout << maxeovertot << endl;
	  //if(maxeovertot > 0.7) continue;
	  //seedD[njet] = maxeovertot;
	  //jet_ph[njet] = (jet_ph[njet]>0?jet_ph[njet]-M_PI:jet_ph[njet]+M_PI);
	  _dPhiLayer[njet] = emAxis.Phi() - ohAxis.Phi();
	  if(_dPhiLayer[njet] > M_PI) _dPhiLayer[njet] -= 2*M_PI;
	  if(_dPhiLayer[njet] < -M_PI) _dPhiLayer[njet] += 2*M_PI;
	  ++njet;
	}
    }
  else
    {
      if(_debug > 1) cout << "No jet node!" << endl;
    }
  //int jetfired = ((1 & triggervec >> 20) | (1 & triggervec >> 21) | (1 & triggervec >> 22) | (1 & triggervec >> 23));
  //int phtfired = ((1 & triggervec >> 28) | (1 & triggervec >> 29) | (1 & triggervec >> 30) | (1 & triggervec >> 31));
  if(!njet && !_dotow) return Fun4AllReturnCodes::EVENT_OK;
  
  if(_debug > 1) cout << "Getting retowered EMCal towers: " << endl;
  if(_dotow) 
    {
      if(rtem)
	{ //get EMCal values
	  int nchannels = 1536; //channels in emcal
	  for(int i=0; i<nchannels; ++i) //loop over channels 
	    {
	      TowerInfo *tower = rtem->get_tower_at_channel(i); //get EMCal tower
	      if(tower->get_isHot() || tower->get_isBadChi2()) continue;
	      rtemen[sector_rtem] = tower->get_energy(); //actual tower energy (calibrated)
	      //cout << rtemen[sector_rtem] << " " << tower->get_energy() << endl;
	      int key = rtem->encode_key(i);
	      rtemet[sector_rtem] = rtem->getTowerEtaBin(key);
	      rtemph[sector_rtem] = rtem->getTowerPhiBin(key);
	      sector_rtem++;
	      //cout << sector_rtem << endl;
	    }
	  if(_debug > 1) cout << sector_rtem << endl;
	}
      else
	{
	  if(_debug > 1) cout << "No retowered towers!" << endl;
	}
    }
  
  emetot = 0;
  ohetot = 0;
  ihetot = 0;
  if(_debug > 1) cout << "Getting EMCal info" << endl;
  /*
  if(clusters)
    {
      RawClusterContainer::ConstRange begin_end = clusters->getClusters();
      _cluster_n = 0;
      for (RawClusterContainer::ConstIterator rtiter = begin_end.first; rtiter != begin_end.second; ++rtiter)
	{
	  RawCluster *cluster = rtiter->second;
	  if (cluster->get_energy() < 0.5) continue;
	  //if (cluster->get_chi2() > 10) continue;                                                                                                                                                             
	  /*
	    float best_E = -1;
	    float best_chi2 = -1;
	    
	    RawCluster::TowerConstRange begin_end_Towers = cluster->get_towers();
	    for (RawCluster::TowerConstIterator iter = begin_end_Towers.first; iter != begin_end_Towers.second; ++iter) {
	    
	    int index1 = RawTowerDefs::decode_index1( iter->first );
	    int index2 = RawTowerDefs::decode_index2( iter->first );
	    unsigned int key = TowerInfoDefs::encode_emcal( index1, index2 );
	    TowerInfo* tower = towers_EM->get_tower_at_key( key );
	    
	    if ( tower->get_energy() > best_E ) {
	    best_E = tower->get_energy();
	    best_chi2 = tower->get_chi2();
	    }
	    
	    }
  
	  //if (best_chi2 > 1e+5) continue;                                                                                                                                                                     
	  
	  //float theta = atan2( cluster->get_r() , cluster->get_z() );                                                                                                                                         
	  //float eta = -1 * log( tan( theta / 2.0 ) );                                                                                                                                                         
	  
	  //std::cout << " -> cluster with E / eta / phi = " << cluster->get_energy() << " / " << eta << " / " << cluster->get_phi() << std::endl;                                                              
	  
	  //_cluster_template_chi2[ _cluster_n ] = best_chi2;
	  
	  //_cluster_chi2[_cluster_n] = cluster->get_chi2();
	  _cluster_nTower[_cluster_n] = cluster->getNTowers();
	  
	  _cluster_Ecore[_cluster_n] = cluster->get_ecore();
	  _cluster_E[_cluster_n] = cluster->get_energy();
	  
	  _cluster_phi[_cluster_n] = cluster->get_phi();
	  _cluster_eta[_cluster_n] = asinh((cluster->get_z()-vtx[2])/cluster->get_r());
	  
	  
	  //std::cout << _cluster_z[_cluster_n] << " " << _cluster_r[_cluster_n] << " " << _cluster_E[_cluster_n] << std::endl;                                                                                 
	  
	  _cluster_n++;
	  
	}
      
    }
  else
    {
      if(_debug > 0) cout << "No cluster info!" << endl;
    }
  */
  if(towersEM && _dotow)
    { //get EMCal values
      int nchannels = 24576; //channels in emcal
      int nover = 0;
      int nneg = 0;
      int nlarge = 0;
      int nzero = 0;
      for(int i=0; i<nchannels; ++i) //loop over channels 
	{
	  TowerInfo *tower = towersEM->get_tower_at_channel(i); //get EMCal tower
	  int key = towersEM->encode_key(i);
	  int etabin, phibin;
	  if(_dotow)
	    {
	      etabin = towersEM->getTowerEtaBin(key);
	      phibin = towersEM->getTowerPhiBin(key);
	    }
	  //float time = towersEM->get_tower_at_channel(i)->get_time_float(); //get time
	  if (tower->get_isHot() || tower->get_isBadChi2())
	    { 
	      //if(_debug) cout << etabin << " " << phibin << " " << tower->get_isBadChi2() << tower->get_energy() << endl;
	      continue;
	    }
	  emcalen[sectorem] = tower->get_energy(); //actual tower energy (calibrated)
	  //emcalchi2[sectorem] = tower->get_chi2();
	  /*
	  if(emcalen[sectorem] == 0) nzero++;
	  if(emcalen[sectorem] > 0.1)
	    {
	      //if(_debug > 1) cout << emcalen[sectorem] << endl;
	      nlarge++;
	    }
	  
	  if(emcalen[sectorem] > 0.005) nover++;
	  if(emcalen[sectorem] < -0.005) nneg++;
	  */
	  emetot += emcalen[sectorem];
	  //emcalt[sectorem] = time; //store time value
	  if(_dotow) 
	    {
	      emcaletabin[sectorem] = etabin; //get eta and phi of towers
	      emcalphibin[sectorem] = phibin;
	      sectorem++;
	    }
	}
      /*
      if(_debug > 1)
	{
	  cout << "Nlarge: " << nlarge << endl;
	  cout << "Nover: " << nover << endl;
	  cout << "NNeg: " << nneg << endl;
	  cout << "Nzero: " <<nzero << endl;
	}
      if(_debug > 1) cout << sectorem << endl;
      */
    }
  else if(_dotow)
    {
      for(int i=0; i<16; ++i)
	{
	  emcalen[i] = -1;
	  emcalt[i] = -1;
	}
    }
  if(_debug > 1) cout << "total EMCal E: " << emetot << endl;

  if(towersIH && _dotow)
    { //get IHCal values
      int nchannels = 1536; //channels in ihcal
      for(int i=0; i<nchannels; ++i) //loop over channels 
	{
	  TowerInfo *tower = towersIH->get_tower_at_channel(i); //get IHCal tower
	  int key = towersIH->encode_key(i);
	  int etabin, phibin;
	  if(_dotow) 
	    {
	      etabin = towersIH->getTowerEtaBin(key);
	      phibin = towersIH->getTowerPhiBin(key);
	    }
	  //float time = towersIH->get_tower_at_channel(i)->get_time_float(); //get time
	  if (tower->get_isHot() || tower->get_isBadChi2())
	    { 
	      continue;
	    }
	  ihcalen[sectorih] = tower->get_energy(); //actual tower energy (calibrated)
	  ihetot += ihcalen[sectorih];
	  /*
	  ihcalchi2[sectorih] = tower->get_chi2();
	  ihcalt[sectorih] = time; //store time value
	  */
	  if(_dotow)
	    {
	      ihcaletabin[sectorih] = etabin; //get eta and phi of towers
	      ihcalphibin[sectorih] = phibin;
	      sectorih++;
	    }
	}
      if(_debug > 1) cout << sectorih << endl;
    }
  else if(_dotow)
    {
      for(int i=0; i<16; ++i)
	{
	  ihcalen[i] = -1;
	  ihcalt[i] = -1;
	}
    }
  if(_debug > 1) cout << "getting OHCal info" << endl;
  if(towersOH && _dotow)
    { //get OHCal values
      int nchannels = 1536; //channels in ohcal
      for(int i=0; i<nchannels; ++i) //loop over channels 
	{
	  TowerInfo *tower = towersOH->get_tower_at_channel(i); //get OHCal tower
	  int key = towersOH->encode_key(i);
	  int etabin, phibin;
	  if(_dotow) 
	    {
	      etabin = towersOH->getTowerEtaBin(key);
	      phibin = towersOH->getTowerPhiBin(key);
	    }
	  //float time = towersOH->get_tower_at_channel(i)->get_time_float(); //get time
	  if (tower->get_isHot() || tower->get_isBadChi2())
	    { 
	      continue;
	    }
	  ohcalen[sectoroh] = tower->get_energy(); //actual tower energy (calibrated)
	  ohetot += ohcalen[sectoroh];
	  //ohcalchi2[sectoroh] = tower->get_chi2();
	  //ohcalt[sectoroh] = time; //store time value
	  if(_dotow) 
	    {
	      ohcaletabin[sectoroh] = etabin; //get eta and phi of towers
	      ohcalphibin[sectoroh] = phibin;
	      sectoroh++;
	    }
	}
      if(_debug > 1) cout << sectoroh << endl;
    }
  else if(_dotow)
    {
      for(int i=0; i<16; ++i)
	{
	  ohcalen[i] = -1;
	  ohcalt[i] = -1;
	}
    }
  //if(_debug > 1) cout << "Getting MBD info" << endl;
  

  if(_debug > 1) cout << "Filling" << endl;
  //if(_debug) cout << rtemen[1535] << " " <<rtemen[sector_rtem-1] << endl;
  _tree->Fill();
  _jett->Fill();
  return Fun4AllReturnCodes::EVENT_OK;
    
}
//____________________________________________________________________________..
int R24earlytreemaker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "R24earlytreemaker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int R24earlytreemaker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "R24earlytreemaker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int R24earlytreemaker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "R24earlytreemaker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  _tree2->Fill();
  _f->cd();
  //_jett->Write();
  _tree->Write();
  _tree2->Write();
  _f->Write();
  _f->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int R24earlytreemaker::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "R24earlytreemaker::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void R24earlytreemaker::Print(const std::string &what) const
{
  std::cout << "R24earlytreemaker::Print(const std::string &what) const Printing info for " << what << std::endl;
}
