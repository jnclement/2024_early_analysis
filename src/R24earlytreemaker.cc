#include "R24earlytreemaker.h"
#include <ffaobjects/EventHeaderv1.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
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

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <iostream>

#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <ffarawobjects/Gl1Packetv2.h>
#include <jetbase/JetMapv1.h>
#include <jetbase/JetMap.h>
using namespace std;
//____________________________________________________________________________..
R24earlytreemaker::R24earlytreemaker(const std::string &name, const int debug, int datorsim, int dotow):
  SubsysReco("test")//).c_str())
{
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
  _tree->Branch("triggervec",&triggervec,"triggervec/g");
  _tree2->Branch("mbevt",&mbevt,"mbevt/I");
  _tree->Branch("ismb",&ismb,"ismb/I");
  _jett->Branch("ismb",&ismb,"ismb/I");
  _jett->Branch("allcomp",&allcomp,"allcomp/I");
  _jett->Branch("alcet",alcet,"alcet[allcomp]/F");
  _tree->Branch("emetot",&emetot,"emetot/F");
  _tree->Branch("ihetot",&ihetot,"ihetot/F");
  _tree->Branch("ohetot",&ohetot,"ohetot/F");
  if(_dotow) 
    {
      //_tree->Branch("sectorem",&sectorem,"sectorem/I"); //Number of hit sectors in the emcal
      _tree->Branch("sectorih",&sectorih,"sectorih/I"); // IHcal etc.
      _tree->Branch("sectoroh",&sectoroh,"sectoroh/I");
      //_tree->Branch("sectoremuc",&sectoremuc,"sectoremuc/I");
      //_tree->Branch("emcalen",emcalen,"emcalen[sectorem]/F"); //energy per EMCal sector
      _tree->Branch("ihcalen",ihcalen,"ihcalen[sectorih]/F"); // per IHCal sector (etc.)
      _tree->Branch("ohcalen",ohcalen,"ohcalen[sectoroh]/F");
      //_tree->Branch("emcalchi2",emcalchi2,"emcalchi2[sectorem]/F"); //energy per EMCal sector
      //_tree->Branch("ihcalchi2",ihcalchi2,"ihcalchi2[sectorih]/F"); // per IHCal sector (etc.)
      //_tree->Branch("ohcalchi2",ohcalchi2,"ohcalchi2[sectoroh]/F");
      //_tree->Branch("emcalenuc",emcalenuc,"emcalenuc[sectoremuc]/F");
      //_tree->Branch("emcaletabin",emcaletabin,"emcaletabin[sectorem]/I"); //eta of EMCal sector
      _tree->Branch("ihcaletabin",ihcaletabin,"ihcaletabin[sectorih]/I");
      _tree->Branch("ohcaletabin",ohcaletabin,"ohcaletabin[sectoroh]/I");
      //_tree->Branch("emcalphibin",emcalphibin,"emcalphibin[sectorem]/I"); //phi of EMCal sector
      _tree->Branch("ihcalphibin",ihcalphibin,"ihcalphibin[sectorih]/I");
      _tree->Branch("ohcalphibin",ohcalphibin,"ohcalphibin[sectoroh]/I");
    }
  _tree->Branch("sectormb",&sectormb,"sectormb/I");
  _tree->Branch("mbenrgy",mbenrgy,"mbenrgy[sectormb]/F"); //MBD reported value (could be charge or time)
  //_tree->Branch("emcalt",emcalt,"emcalt[sectorem]/F"); //time value of EMCal sector
  //_tree->Branch("ihcalt",ihcalt,"ihcalt[sectorih]/F");
  //_tree->Branch("ohcalt",ohcalt,"ohcalt[sectoroh]/F");
  _tree->Branch("vtx",vtx,"vtx[3]/F");
  if(_dotow) _tree->Branch("sector_rtem",&sector_rtem,"sector_rtem/I");
  _tree->Branch("njet",&njet,"njet/I");
  _jett->Branch("njet",&njet,"njet/I");
  _jett->Branch("aceta",aceta,"aceta[njet]/F");
  _tree->Branch("seedD",&seedD,"seedD[njet]/F");
  _tree->Branch("jet_e",jet_e,"jet_e[njet]/F");
  _jett->Branch("jet_e",jet_e,"jet_e[njet]/F");
  //_tree->Branch("jet_r",jet_r,"jet_r[njet]/F");
  _tree->Branch("jet_et",jet_et,"jet_et[njet]/F");
  _tree->Branch("jet_ph",jet_ph,"jet_ph[njet]/F");

  _jett->Branch("jet_et",jet_et,"jet_et[njet]/F");
  _jett->Branch("jet_ph",jet_ph,"jet_ph[njet]/F");

  if(_dotow)
    {
      _tree->Branch("rtemen",rtemen,"rtemen[sector_rtem]/F");
      _tree->Branch("rtemet",rtemet,"rtemet[sector_rtem]/I");
      _tree->Branch("rtemph",rtemph,"rtemph[sector_rtem]/I");
    }
  _tree->Branch("ehjet",ehjet,"ehjet[njet]/F");
  _tree->Branch("_evtnum",&_evtnum,"_evtnum/I");
  _tree->Branch("cluster_n",&_cluster_n,"cluster_n/I");
  _tree->Branch("cluster_E",_cluster_E,"cluster_E[cluster_n]/F");
  _tree->Branch("cluster_Ecore",_cluster_Ecore,"cluster_Ecore[cluster_n]/F");
  _tree->Branch("cluster_phi",_cluster_phi,"cluster_phi[cluster_n]/F");
  _tree->Branch("cluster_eta",_cluster_eta,"cluster_eta[cluster_n]/F");
  //_tree->Branch("cluster_r",_cluster_r,"cluster_r[cluster_n]/F");
  //_tree->Branch("cluster_chi2",_cluster_chi2,"cluster_chi2[cluster_n]/F");
  //_tree->Branch("cluster_template_chi2",_cluster_template_chi2,"cluster_template_chi2[cluster_n]/F");
  _tree->Branch("cluster_nTower",_cluster_nTower,"cluster_nTower[cluster_n]/I");
  if(!_datorsim) _tree->Branch("ntj",&ntj,"ntj/I");
  if(!_datorsim) _tree->Branch("tjet_e",tjet_e,"tjet_e[ntj]/F");
  if(!_datorsim) _tree->Branch("tjet_eta",tjet_eta,"tjet_eta[ntj]/F");
  if(!_datorsim) _tree->Branch("tjet_phi",tjet_phi,"tjet_phi[ntj]/F");
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
  TowerInfoContainer *towersEM = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC");
  if(!towersEM) towersEM = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC");
  TowerInfoContainer *towersIH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALIN");
  if(!towersIH) towersIH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *towersOH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALOUT");
  if(!towersOH) towersOH = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_HCALOUT");
  TowerInfoContainer *towersEMuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_CEMC");
  TowerInfoContainer *towersIHuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALIN");
  TowerInfoContainer *towersOHuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALOUT");
  //TowerInfoContainer *towersZD = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_ZDC");
  TowerInfoContainer *rtem = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC");
  if(!rtem) rtem = findNode::getClass<TowerInfoContainerv1>(topNode, "TOWERINFO_CALIB_CEMC");
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

  if(mbdtow)
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
  int isjettrig = ((triggervec >> 16) & 1) | ((triggervec >> 17) & 1) | ((triggervec >> 18) & 1) | ((triggervec >> 19) & 1);
  //if(!ismb &! isjettrig) return Fun4AllReturnCodes::EVENT_OK;
  ntj = 0;
  if(truthjets && !_datorsim)
    {
      for(int i=0; i<truthjets->size(); ++i)
	{
	  Jet* jet = truthjets->get_jet(i);
	  tjet_e[ntj] = jet->get_e();
	  if(tjet_e[ntj] < 4) continue;
	  tjet_eta[ntj] = jet->get_eta();
	  tjet_phi[ntj] = jet->get_phi();
	  ntj++;
	}
    }

  vtx[0] = 0;
  vtx[1] = 0;
  vtx[2] = NAN;
  if(_datorsim)
    {
      for(auto iter = mbdvtxmap->begin(); iter != mbdvtxmap->end(); ++iter)
	{
	  MbdVertex* mbdvtx = iter->second;
	  vtx[2] = mbdvtx->get_z();
	  break;
	}
    }
  else
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
  if(std::isnan(vtx[2]) || abs(vtx[2]) > 1000)
    {
      if(ismb) mbevt--;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  /*
  else
    {
      ismb = 1;
    }
  */
  allcomp = 0;
  if(_debug > 1) cout << "Getting jets: " << endl;
  if(jets)
    {
      int tocheck = jets->size();
      if(_debug > 2) cout << "Found " << tocheck << " jets to check..." << endl;
      for(int i=0; i<tocheck; ++i)
	{
	  Jet *jet = jets->get_jet(i);
	  if(jet)
	    {
	      jet_e[njet] = jet->get_e();
	      jet_r[njet] = 0.4;
	      jet_et[njet] = jet->get_eta();
	      jet_ph[njet] = jet->get_phi();
	    }
	  else
	    {
	      continue;
	    }
	  if(jet_e[njet] < 4) continue;
	  if(_debug > 2) cout << "found a good jet!" << endl;
	  float maxeovertot = 0;
	  float hcale = 0;
	  float ecale = 0;
	  int ncomp = 0;
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
		  int key = towersIH->encode_key(channel);
		  unsigned int etabin = towersIH->getTowerEtaBin(key);
		  aceta[njet] += (etabin-11.5)*1.1/12;
		  alcet[allcomp] = (etabin-11.5)*1.1/12;
		}
	      else if(comp.first == 7 || comp.first == 27)
		{
		  if(_debug > 2) cout << "ohfail" << endl;
		  tower = towersOH->get_tower_at_channel(channel);
		  if(_debug > 2) cout << "towerad: " << tower << endl;
		  hcale += tower->get_energy();
		  int key = towersOH->encode_key(channel);
		  unsigned int etabin = towersOH->getTowerEtaBin(key);
		  aceta[njet] += (etabin-11.5)*1.1/12;
		  alcet[allcomp] = (etabin-11.5)*1.1/12;
		}
	      else if(comp.first == 13 || comp.first == 28 || comp.first == 25)
		{
		  if(_debug > 2) cout << "emfail" << endl;
		  tower = rtem->get_tower_at_channel(channel);
		  if(_debug > 2) cout << "towerad: " << tower << endl;
		  ecale += tower->get_energy();
		  int key = rtem->encode_key(channel);
		  unsigned int etabin = rtem->getTowerEtaBin(key);
		  aceta[njet] += (etabin-47.5)*1.1/48;
		  alcet[allcomp] = (etabin-47.5)*1.1/48;
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
	  aceta[njet] /= ncomp;
	  //aceta[njet] /= jet_e[njet];
	  if(_debug > 2) cout << "Now filling some jet adjacent numbers: " << endl;
	  if(_debug && jet_e[njet] - ecale - hcale > 0.01) cout << jet_e[njet] << " " << ecale << " " << hcale << endl;
	  ehjet[njet] = ecale/hcale;
	  maxeovertot /= jet_e[njet];
	  //if(_debug) cout << maxeovertot << endl;
	  //if(maxeovertot > 0.7) continue;
	  seedD[njet] = maxeovertot;
	  jet_ph[njet] = (jet_ph[njet]>0?jet_ph[njet]-M_PI:jet_ph[njet]+M_PI);
	  ++njet;
	}
    }
  else
    {
      if(_debug > 1) cout << "No jet node!" << endl;
    }
  //int jetfired = ((1 & triggervec >> 20) | (1 & triggervec >> 21) | (1 & triggervec >> 22) | (1 & triggervec >> 23));
  //int phtfired = ((1 & triggervec >> 28) | (1 & triggervec >> 29) | (1 & triggervec >> 30) | (1 & triggervec >> 31));
  //if(!njet)return Fun4AllReturnCodes::EVENT_OK;

  if(_debug > 1) cout << "Getting retowered EMCal towers: " << endl;
  if(_dotow) 
    {
      if(rtem)
	{ //get EMCal values
	  int nchannels = 24576; //channels in emcal
	  for(int i=0; i<nchannels; ++i) //loop over channels 
	    {
	      TowerInfo *tower = rtem->get_tower_at_channel(i); //get EMCal tower
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
	  */
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
  if(towersEM)
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
  else
    {
      for(int i=0; i<16; ++i)
	{
	  emcalen[i] = -1;
	  emcalt[i] = -1;
	}
    }
  if(_debug > 1) cout << "total EMCal E: " << emetot << endl;

  if(towersIH)
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
  else
    {
      for(int i=0; i<16; ++i)
	{
	  ihcalen[i] = -1;
	  ihcalt[i] = -1;
	}
    }
  if(_debug > 1) cout << "getting OHCal info" << endl;
  if(towersOH)
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
  else
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
  _jett->Write();
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
