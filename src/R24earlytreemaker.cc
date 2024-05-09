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
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertex.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <mbd/MbdPmtContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <globalvertex/MbdVertexMap.h>
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

using namespace std;
//____________________________________________________________________________..
R24earlytreemaker::R24earlytreemaker(const std::string &name, const int debug):
  SubsysReco("test")//).c_str())
{
  _evtct = 0;
  _foutname = name;
  _debug = debug;
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
  
  _tree = new TTree("ttree","a persevering date tree");
  //_tree->Branch("emetot",&emetot,"emetot/F");
  //_tree->Branch("ihetot",&ihetot,"ihetot/F");
  //_tree->Branch("ohetot",&ohetot,"ohetot/F");
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
  //_tree->Branch("sectormb",&sectormb,"sectormb/I");
  //_tree->Branch("mbenrgy",mbenrgy,"mbenrgy[sectormb]/F"); //MBD reported value (could be charge or time)
  //_tree->Branch("emcalt",emcalt,"emcalt[sectorem]/F"); //time value of EMCal sector
  //_tree->Branch("ihcalt",ihcalt,"ihcalt[sectorih]/F");
  //_tree->Branch("ohcalt",ohcalt,"ohcalt[sectoroh]/F");
  //_tree->Branch("vtx",vtx,"vtx[3]/F");
  _tree->Branch("sector_rtem",&sector_rtem,"sector_rtem/I");
  _tree->Branch("njet",&njet,"njet/I");
  _tree->Branch("seedD",&seedD,"seedD[njet]/F");
  _tree->Branch("jet_e",jet_e,"jet_e[njet]/F");
  //_tree->Branch("jet_r",jet_r,"jet_r[njet]/F");
  _tree->Branch("jet_et",jet_et,"jet_et[njet]/F");
  _tree->Branch("jet_ph",jet_ph,"jet_ph[njet]/F");
  _tree->Branch("rtemen",rtemen,"rtemen[sector_rtem]/F");
  _tree->Branch("rtemet",rtemet,"rtemet[sector_rtem]/I");
  _tree->Branch("rtemph",rtemph,"rtemph[sector_rtem]/I");
  _tree->Branch("ehjet",ehjet,"ehjet[njet]/F");
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
  _evtct++;
  if(_debug > 1) cout << endl << endl << endl << "Beginning event processing" << endl;
  if(_debug > 1) cout << "Event " << _evtct << endl;
  float mbdq = 0;
  //reset lengths to all zero
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
  TowerInfoContainer *towersIH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALIN");
  TowerInfoContainer *towersOH = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_HCALOUT");
  TowerInfoContainer *towersEMuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_CEMC");
  TowerInfoContainer *towersIHuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALIN");
  TowerInfoContainer *towersOHuc = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_HCALOUT");
  TowerInfoContainer *towersZD = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERS_ZDC");
  TowerInfoContainer *rtem = findNode::getClass<TowerInfoContainerv2>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  JetContainer *jets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_Tower_HIRecoSeedsRaw_r02");
  if(_debug > 1) cout << "Getting jets: " << endl;
  if(jets)
    {
      int tocheck = jets->size();
      for(int i=0; i<tocheck; ++i)
	{
	  Jet *jet = jets->get_jet(i);
	  if(jet)
	    {
	      jet_e[njet] = jet->get_e();
	      jet_r[njet] = 0.2;
	      jet_et[njet] = jet->get_eta();
	      jet_ph[njet] = jet->get_phi();
	    }
	  if(jet_e[njet] < 4) continue;
	  float maxeovertot = 0;
	  float hcale = 0;
	  float ecale = 0;
	  for(auto comp: jet->get_comp_vec())
	    {
	      unsigned int channel = comp.second;
	      TowerInfo* tower;
	      if(comp.first == 5 || comp.first == 26)
		{
		  tower = towersIH->get_tower_at_channel(channel);
		  hcale += tower->get_energy();
		}
	      else if(comp.first == 7 || comp.first == 27)
		{
		  tower = towersOH->get_tower_at_channel(channel);
		  hcale += tower->get_energy();
		}
	      else if(comp.first == 13 || comp.first == 28)
		{
		  tower = rtem->get_tower_at_channel(channel);
		  ecale += tower->get_energy();
		}
	      else
		{
		  if(_debug > 1) cout << "No good detector for getting jet components" << endl;
		  continue;
		}
	      float eval = tower->get_energy();
	      if(eval > maxeovertot) maxeovertot = eval;
	    }
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

  if(!njet) return Fun4AllReturnCodes::EVENT_OK;

  if(_debug > 1) cout << "Getting retowered EMCal towers: " << endl;

  if(rtem)
    { //get EMCal values
      int nchannels = 1536; //channels in emcal
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

  if(towersEMuc)
    { //get EMCal values
      int nchannels = 24576; //channels in emcal
      int nneg = 0;
      int nover = 0;
      int nint = 0;
      for(int i=0; i<nchannels; ++i) //loop over channels 
	{
	  TowerInfo *tower = towersEMuc->get_tower_at_channel(i); //get EMCal tower
	  emcalenuc[sectoremuc] = tower->get_energy(); //actual tower energy (calibrated)
	  //if(_debug > 1 && (emcalenuc[sectoremuc] > 10)) cout << emcalenuc[sectoremuc] << endl;
	  if(emcalenuc[sectoremuc] < 0) nneg++;
	  if(emcalenuc[sectoremuc] > 10) nover++;
	  if(emcalenuc[sectoremuc] - floor(emcalenuc[sectoremuc]) == 0)
	    {
	      nint++;
	    } 
	  sectoremuc++;
	}
      if(_debug > 1)
	{
	  cout << "NNegADC: " << nneg << endl;
	  cout << "NoverADC: " << nover << endl;
	  cout << "Nint: " << nint << endl;
	}
      if(_debug > 1) cout << sectoremuc << endl;
    }
  

  if(towersZD)
    {
      int nchannels = 16;
      for(int i = 0; i < nchannels; ++i)
	{
	  TowerInfo *tower = towersZD->get_tower_at_channel(i);
	  float time = tower->get_time_float(); //get time
	  if (tower->get_isHot() || tower->get_isBadChi2()) 
	    { 
	      continue;
	    }
	  zdcalen[sectorzd] = tower->get_energy(); //actual tower energy (calibrated)
	  zdcalt[sectorzd] = time; //store time value
	  sectorzd++;
	}
      if(_debug > 1) cout << sectorzd << endl;
    }
  else
    {
      for(int i=0; i<16; ++i)
	{
	  zdcalen[i] = -1;
	  zdcalt[i] = -1;
	}
    }
  emetot = 0;
  ohetot = 0;
  ihetot = 0;
  if(_debug > 1) cout << "Getting EMCal info" << endl;
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
	  int etabin = towersEM->getTowerEtaBin(key);
	  int phibin = towersEM->getTowerPhiBin(key);
	  float time = towersEM->get_tower_at_channel(i)->get_time_float(); //get time
	  if (tower->get_isHot() || tower->get_isBadChi2()) { 
	    //if(_debug) cout << etabin << " " << phibin << " " << tower->get_isBadChi2() << tower->get_energy() << endl;
	    continue;
	  }
	  emcalen[sectorem] = tower->get_energy(); //actual tower energy (calibrated)
	  emcalchi2[sectorem] = tower->get_chi2();
	  if(emcalen[sectorem] == 0) nzero++;
	  if(emcalen[sectorem] > 0.1)
	    {
	      //if(_debug > 1) cout << emcalen[sectorem] << endl;
	      nlarge++;
	    }
	  if(emcalen[sectorem] > 0.005) nover++;
	  if(emcalen[sectorem] < -0.005) nneg++;
	  emetot += emcalen[sectorem];
	  emcalt[sectorem] = time; //store time value
	  emcaletabin[sectorem] = etabin; //get eta and phi of towers
	  emcalphibin[sectorem] = phibin;
	  sectorem++;
	}
      if(_debug > 1)
	{
	  cout << "Nlarge: " << nlarge << endl;
	  cout << "Nover: " << nover << endl;
	  cout << "NNeg: " << nneg << endl;
	  cout << "Nzero: " <<nzero << endl;
	}
      if(_debug > 1) cout << sectorem << endl;
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
	  int etabin = towersIH->getTowerEtaBin(key);
	  int phibin = towersIH->getTowerPhiBin(key);
	  float time = towersIH->get_tower_at_channel(i)->get_time_float(); //get time
	  if (tower->get_isHot() || tower->get_isBadChi2()) { 
	    continue;
	  }
	  ihcalen[sectorih] = tower->get_energy(); //actual tower energy (calibrated)
	  ihetot += ihcalen[sectorih];
	  ihcalchi2[sectorih] = tower->get_chi2();
	  ihcalt[sectorih] = time; //store time value
	  ihcaletabin[sectorih] = etabin; //get eta and phi of towers
	  ihcalphibin[sectorih] = phibin;
	  sectorih++;
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
	  int etabin = towersOH->getTowerEtaBin(key);
	  int phibin = towersOH->getTowerPhiBin(key);
	  float time = towersOH->get_tower_at_channel(i)->get_time_float(); //get time
	  if (tower->get_isHot() || tower->get_isBadChi2()) { 
	    continue;
	  }
	  ohcalen[sectoroh] = tower->get_energy(); //actual tower energy (calibrated)
	  ohetot += ohcalen[sectoroh];
	  ohcalchi2[sectoroh] = tower->get_chi2();
	  ohcalt[sectoroh] = time; //store time value
	  ohcaletabin[sectoroh] = etabin; //get eta and phi of towers
	  ohcalphibin[sectoroh] = phibin;
	  sectoroh++;
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
  if(_debug > 1) cout << "Getting MBD info" << endl;
  
  MbdPmtContainer *mbdtow = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  if(mbdtow)
  {
    sectormb = 128;//mbdtow->get_npmt();
    //if(_debug) cout << "Got " << sectormb << " mbd sectors in sim." << endl;
    for(int i=0; i<sectormb; ++i)
      {
	MbdPmtHit *mbdhit = mbdtow->get_pmt(i);
	//if(_debug > 2) cout << "PMT " << i << " address: " << mbdhit << " charge: " << mbdhit->get_q() << endl;
	mbenrgy[i] = mbdhit->get_q();
	mbdq += mbdhit->get_q();
      }
  }
  else
    {
      for(int i=0; i<128; ++i)
	{
	  mbenrgy[i] = -1;
	  if(_debug > 1) cout << "No MBD info!" << endl;
	}
    }
  if(_debug > 1) cout << "Filling" << endl;
  //if(_debug) cout << rtemen[1535] << " " <<rtemen[sector_rtem-1] << endl;
  _tree->Fill();
  
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

  _f->cd();
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
