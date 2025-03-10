#ifndef R24TREEMAKER_H
#define R24TREEMAKER_H

#include <fun4all/SubsysReco.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include <phool/recoConsts.h>
#include <phparameter/PHParameters.h>
class PHCompositeNode;
class CentralityInfo;
class R24earlytreemaker : public SubsysReco
{
 public:

  R24earlytreemaker(const std::string &name = "R24earlytreemaker", const int debug = 0, int datorsim = 1, int dotow = 0, int sampleType = 1);

  virtual ~R24earlytreemaker();

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;


 private:
  PHParameters _cutParams;
  int ismb = 0;
  long unsigned int _bbfqavec;
  float _frcem[1000];
  float _frcoh[1000];
  int _cluster_n;
  int _cluster_nTower[10000];
  float _cluster_Ecore[10000];
  float _cluster_E[10000];
  float _cluster_phi[10000];
  float _cluster_eta[10000];
  int _dotow;
  int _datorsim;
  int _evtnum;
  int _evtct;
  int _debug;
  TFile *_f;
  TTree *_tree;
  TTree *_tree2;
  TTree *_jett;
  int failscut;
  std::string _foutname;
  float emetot, ihetot, ohetot;
  int sectorem;
  int sectorih;
  int sectoroh;
  int sectormb;
  int sectorzd;
  int sectoremuc;
  int njet;
  int allcomp;
  long long unsigned int triggervec;
  int mbevt;
  int sector_rtem;
  float seedD[1000];
  float jet_e[1000];
  float jet_pt[1000];
  float jet_r[1000];
  float jet_et[1000];
  float jet_ph[1000];
  int ntj;
  float tjet_pt[1000];
  float tjet_et[1000];
  float tjet_eta[1000];
  float tjet_phi[1000];
  float ehjet[1000];
  float aceta[1000];
  float alcet[24576];
  float rtemen[24576];
  int rtemet[24576];
  int rtemph[24576];
  float emcalen[24576];
  float ihcalen[1536];
  float ohcalen[1536];
  float emcalchi2[24576];
  float ihcalchi2[1536];
  float ohcalchi2[1536];
  float zdcalen[16];
  float emcalenuc[24576];
  float emcalt[24576];
  float ihcalt[1536];
  float ohcalt[1536];
  float zdcalt[16];
  float maxTowerET;
  float subTowerET;
  float maxTowerChi2[3]; 
  float caloEfrac[3];
  int _sampleType;
  float _dPhi2pc[1000];
  float _dEta2pc[1000];
  /*
  float _emLayerJetPhi[10];
  float _ohLayerJetPhi[10];
  float _emLayerJetEta[10];
  float _ohLayerJetEta[10];
  float _emLayerJetET[10];
  float _ohLayerJetET[10];
  int _nLayerEm;
  int _nLayerOh;
  */
  int _n2pc;
  float _dPhiLayer[10];
  float _l2pcEta;
  //int emcaladc[24576];
  //int ihcaladc[1536];
  //int ohcaladc[1536];
  //int emcalzsadc[24576];
  //int ihcalzsadc[1536];
  //int ohcalzsadc[1536];
  //float emcalzs[24576];
  //float ihcalzs[1536];
  //float ohcalzs[1536];
  //float emcalpos[24576][3];
  //float ihcalpos[1536][3];
  //float ohcalpos[1536][3];
  int emcaletabin[24576];
  int ihcaletabin[1536];
  int ohcaletabin[1536];
  int emcalphibin[24576];
  int ihcalphibin[1536];
  int ohcalphibin[1536];
  float mbenrgy[256];
  //int centbin;
  //bool isMinBias;
  //int npart;
  //int ncoll;
  //float bimp;
  //float track_vtx[3];
  //float svtx_vtx[3];
  float vtx[3];
  //float truth_vtx[3];
  //float emetacor[24576];
  //float ihetacor[1536];
  //float ohetacor[1536];
  //float emchi2[24576];
  //float ihchi2[1536];
  //float ohchi2[1536];
  //bool emishot[24576];
  //bool ihishot[1536];
  //bool ohishot[1536];
  //int truthpar_n;
  //int truthpar_nh;
  //float truthpar_pz[100000];
  //float truthpar_pt[100000];
  //float truthpar_e[100000];
  //float truthpar_eta[100000];
  //float truthpar_phi[100000];
  //float truthparh_e[100000];
  //float truthparh_pt[100000];
  //float truthparh_pz[100000];
  //float truthparh_eta[100000];
  //float truthparh_phi[100000];
  //int truthparh_id[100000];
  //int hotmap[3][96][256] = {0};
  //float centrality;
  //std::vector<int> baryons{2212,2112,2224,2214,2114,1114,3122,3222,3212,3112,
  //    3224,3214,3114,3322,3312,3324,3314,3334,4122,4222,4212,4112,4224,4214,
  //    4114,4232,4312,4324,4314,4332,4334,4412,4422,4414,4424,4432,4434,4444,
  //    5122,5112,5212,5222,5114,5214,5224,5132,5232,5312,5322,5314,5324,5332,
  //    5334,5142,5242,5412,5422,5414,5424,5342,5432,5434,5442,5444,5512,5522,
  //    5514,5524,5532,5534,5542,5544,5554};
};

#endif // R24TREEMAKER
