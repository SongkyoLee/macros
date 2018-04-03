#ifndef __TRUTHJETTRIGGER_H__
#define __TRUTHJETTRIGGER_H__

// --- need to check all these includes...
#include <fun4all/SubsysReco.h>
#include <vector>

#include "TTree.h"
#include "TFile.h"

class PHCompositeNode;
class PHG4Particle;
namespace HepMC {
class GenParticle;
} /* namespace HepMC */

class KyoTreeMaker: public SubsysReco
{

 public:

  KyoTreeMaker(const std::string &name, bool fill_g4particle_tree=true, bool fill_track_tree=true, bool fill_tower_tree=true, bool fill_cluster_tree=true, bool fill_truthjet4_tree=true, bool fill_towerjet4_tree=true, bool fill_otherjetR_tree=true, bool fill_bhhit_tree=true);

  int Init(PHCompositeNode*);
  int process_event(PHCompositeNode*);
  int End(PHCompositeNode*);

 private:

  float PositionCorrectedNxNTowerE( PHCompositeNode *topNode, unsigned int cl_id, float eNN);
  void getClusterByIndex( PHCompositeNode *topNode, std::string detName, unsigned int cl_id, float teta, float tphi);
  void getClusterPCByIndex( PHCompositeNode *topNode, unsigned int cl_id);

  void GetJetPrimaryContributors( PHCompositeNode *topNode, PHG4Particle* g4ptl);
  void TraceHepMCParticle( HepMC::GenParticle *ptcle);

  TFile *_f;
  std::string _foutname;

  bool _fill_g4particle_tree;
  bool _fill_track_tree;
  bool _fill_tower_tree;
  //bool _fill_cluster_tree;
  bool _fill_cemc_cluster_tree;
  bool _fill_ihcal_cluster_tree;
  bool _fill_ohcal_cluster_tree;
  bool _fill_truthjet2_tree;
  bool _fill_truthjet3_tree;
  bool _fill_truthjet4_tree;
  bool _fill_truthjet5_tree;
  bool _fill_towerjet2_tree;
  bool _fill_towerjet3_tree;
  bool _fill_towerjet4_tree;
  bool _fill_towerjet5_tree;
  bool _fill_bhhit_tree;
  
  TTree *g4particle_tree;
  TTree *track_tree;
  TTree *tower_tree;
  //TTree *cluster_tree;
  TTree *cemc_cluster_tree;
  TTree *ihcal_cluster_tree;
  TTree *ohcal_cluster_tree;
  TTree *truthjet2_tree;
  TTree *truthjet3_tree;
  TTree *truthjet4_tree;
  TTree *truthjet5_tree;
  TTree *towerjet2_tree;
  TTree *towerjet3_tree;
  TTree *towerjet4_tree;
  TTree *towerjet5_tree;
  TTree *bhhit_tree;
  
  ///////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////// 
  int _b_event;
  float _b_vtx_x_1;
  float _b_vtx_y_1;
  float _b_vtx_z_1;
  float _b_vtx_x_2;
  float _b_vtx_y_2;
  float _b_vtx_z_2;
  //float _b_vtx_t;
  
  ///////////////////////////////////////////////////////////////////////// 
  int _b_g4particle_n;
  std::vector<float> _b_g4particle_e;
  std::vector<float> _b_g4particle_p;
  std::vector<float> _b_g4particle_pt;
  std::vector<float> _b_g4particle_eta;
  std::vector<float> _b_g4particle_phi;
  std::vector<int> _b_g4particle_pid;

  ///////////////////////////////////////////////////////////////////////// 
  int _b_track_n;
  std::vector<float> _b_track_p;
  std::vector<float> _b_track_pt;
  std::vector<float> _b_track_eta;
  std::vector<float> _b_track_phi;
  std::vector<int> _b_track_charge;
  std::vector<float> _b_track_chisq;
  std::vector<int> _b_track_ndf;
  std::vector<float> _b_track_dca2d;
  std::vector<float> _b_cemc_E3x3;
  std::vector<float> _b_cemc_E3x3_PC;
  std::vector<float> _b_ihcal_E3x3;
  std::vector<float> _b_ohcal_E3x3;
  std::vector<float> _b_cemc_E5x5;
  std::vector<float> _b_cemc_E5x5_PC;
  std::vector<float> _b_ihcal_E5x5;
  std::vector<float> _b_ohcal_E5x5;
  std::vector<int> _b_cemc_clNtwr;
  std::vector<float> _b_cemc_cldeta;
  std::vector<float> _b_cemc_cldphi;
  std::vector<float> _b_cemc_clE;
  std::vector<float> _b_cemc_clE_PC;
  std::vector<int> _b_ihcal_clNtwr;
  std::vector<float> _b_ihcal_cldeta;
  std::vector<float> _b_ihcal_cldphi;
  std::vector<float> _b_ihcal_clE;
  std::vector<int> _b_ohcal_clNtwr;
  std::vector<float> _b_ohcal_cldeta;
  std::vector<float> _b_ohcal_cldphi;
  std::vector<float> _b_ohcal_clE;
  
  ///////////////////////////////////////////////////////////////////////// 
  int _b_cemc_n;
  std::vector<float> _b_cemc_e;
  std::vector<float> _b_cemc_p;
  std::vector<float> _b_cemc_pt;
  std::vector<float> _b_cemc_eta;
  std::vector<float> _b_cemc_phi;
  std::vector<int> _b_cemc_id;
  int _b_ihcal_n;
  std::vector<float> _b_ihcal_e;
  std::vector<float> _b_ihcal_p;
  std::vector<float> _b_ihcal_pt;
  std::vector<float> _b_ihcal_eta;
  std::vector<float> _b_ihcal_phi;
  int _b_ohcal_n;
  std::vector<float> _b_ohcal_e;
  std::vector<float> _b_ohcal_p;
  std::vector<float> _b_ohcal_pt;
  std::vector<float> _b_ohcal_eta;
  std::vector<float> _b_ohcal_phi;
  
  ///////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////// 
  float _b_cl_cemc_e;
  float _b_cl_cemc_p;
  float _b_cl_cemc_pt;
  float _b_cl_cemc_eta;
  float _b_cl_cemc_phi;
  int _b_cl_cemc_id;
  int _b_cl_cemc_ntwr;
  float _b_cl_cemc_twrEsum;
  float _b_cl_cemc_prob;
  float _b_cl_cemc_chi2;
  std::vector<float> _b_cl_twr_cemc_e;
  std::vector<int> _b_cl_twr_cemc_id;
/*
  float _b_cl_cemcMatched_e;
  float _b_cl_cemcMatched_p;
  float _b_cl_cemcMatched_pt;
  float _b_cl_cemcMatched_eta;
  float _b_cl_cemcMatched_phi;
  int _b_cl_cemcMatched_id;
  int _b_cl_cemcMatched_ntwr;
  float _b_cl_cemcMatched_prob;
  float _b_cl_cemcMatched_chi2;
  float _b_cl_cemcMatched_track_p;
  float _b_cl_cemcMatched_track_deta;
  float _b_cl_cemcMatched_track_dphi;
  float _b_cl_cemcMatched_track_dr;
*/
  ///////////////////////////////////////////////////////////////////////// 
  float _b_cl_ihcal_e;
  float _b_cl_ihcal_p;
  float _b_cl_ihcal_pt;
  float _b_cl_ihcal_eta;
  float _b_cl_ihcal_phi;
  int _b_cl_ihcal_id;
  int _b_cl_ihcal_ntwr;
  std::vector<float> _b_cl_twr_ihcal_e;
  std::vector<int> _b_cl_twr_ihcal_id;
  ///////////////////////////////////////////////////////////////////////// 
  float _b_cl_ohcal_e;
  float _b_cl_ohcal_p;
  float _b_cl_ohcal_pt;
  float _b_cl_ohcal_eta;
  float _b_cl_ohcal_phi;
  int _b_cl_ohcal_id;
  int _b_cl_ohcal_ntwr;
  std::vector<float> _b_cl_twr_ohcal_e;
  std::vector<int> _b_cl_twr_ohcal_id;

  ///////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////// 
  float _b_truthjet2_e;
  float _b_truthjet2_p;
  float _b_truthjet2_pt;
  float _b_truthjet2_eta;
  float _b_truthjet2_phi;
  int _b_truthjet2_cons_n; //constituents
  std::vector<float> _b_truthjet2_cons_e;
  std::vector<float> _b_truthjet2_cons_p;
  std::vector<float> _b_truthjet2_cons_pt;
  std::vector<float> _b_truthjet2_cons_eta;
  std::vector<float> _b_truthjet2_cons_phi;
  std::vector<float> _b_truthjet2_cons_dR;
  std::vector<int> _b_truthjet2_cons_pid;
  int _b_truthjet2_pg_n; //progenitors
  std::vector<int> _b_truthjet2_pg_id;
  std::vector<float> _b_truthjet2_pg_fract;
  std::vector<int> _b_truthjet2_pg_status;
  
  float _b_truthjet3_e;
  float _b_truthjet3_p;
  float _b_truthjet3_pt;
  float _b_truthjet3_eta;
  float _b_truthjet3_phi;
  int _b_truthjet3_cons_n; //constituents
  std::vector<float> _b_truthjet3_cons_e;
  std::vector<float> _b_truthjet3_cons_p;
  std::vector<float> _b_truthjet3_cons_pt;
  std::vector<float> _b_truthjet3_cons_eta;
  std::vector<float> _b_truthjet3_cons_phi;
  std::vector<float> _b_truthjet3_cons_dR;
  std::vector<int> _b_truthjet3_cons_pid;
  int _b_truthjet3_pg_n; //progenitors
  std::vector<int> _b_truthjet3_pg_id;
  std::vector<float> _b_truthjet3_pg_fract;
  std::vector<int> _b_truthjet3_pg_status;
  
  float _b_truthjet4_e;
  float _b_truthjet4_p;
  float _b_truthjet4_pt;
  float _b_truthjet4_eta;
  float _b_truthjet4_phi;
  int _b_truthjet4_cons_n; //constituents
  std::vector<float> _b_truthjet4_cons_e;
  std::vector<float> _b_truthjet4_cons_p;
  std::vector<float> _b_truthjet4_cons_pt;
  std::vector<float> _b_truthjet4_cons_eta;
  std::vector<float> _b_truthjet4_cons_phi;
  std::vector<float> _b_truthjet4_cons_dR;
  std::vector<int> _b_truthjet4_cons_pid;
  int _b_truthjet4_pg_n; //progenitors
  std::vector<int> _b_truthjet4_pg_id;
  std::vector<float> _b_truthjet4_pg_fract;
  std::vector<int> _b_truthjet4_pg_status;
  
  float _b_truthjet5_e;
  float _b_truthjet5_p;
  float _b_truthjet5_pt;
  float _b_truthjet5_eta;
  float _b_truthjet5_phi;
  int _b_truthjet5_cons_n; //constituents
  std::vector<float> _b_truthjet5_cons_e;
  std::vector<float> _b_truthjet5_cons_p;
  std::vector<float> _b_truthjet5_cons_pt;
  std::vector<float> _b_truthjet5_cons_eta;
  std::vector<float> _b_truthjet5_cons_phi;
  std::vector<float> _b_truthjet5_cons_dR;
  std::vector<int> _b_truthjet5_cons_pid;
  int _b_truthjet5_pg_n; //progenitors
  std::vector<int> _b_truthjet5_pg_id;
  std::vector<float> _b_truthjet5_pg_fract;
  std::vector<int> _b_truthjet5_pg_status;

  ///////////////////////////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////////////////////////// 
  float _b_towerjet2_e;
  float _b_towerjet2_p;
  float _b_towerjet2_pt;
  float _b_towerjet2_eta;
  float _b_towerjet2_phi;
  int _b_towerjet2_cemc_n; //cemc towers
  std::vector<float> _b_towerjet2_cemc_e;
  std::vector<float> _b_towerjet2_cemc_p;
  std::vector<float> _b_towerjet2_cemc_pt;
  std::vector<float> _b_towerjet2_cemc_eta;
  std::vector<float> _b_towerjet2_cemc_phi;
  std::vector<float> _b_towerjet2_cemc_dR;
  int _b_towerjet2_ihcal_n; //ihcal towers
  std::vector<float> _b_towerjet2_ihcal_e;
  std::vector<float> _b_towerjet2_ihcal_p;
  std::vector<float> _b_towerjet2_ihcal_pt;
  std::vector<float> _b_towerjet2_ihcal_eta;
  std::vector<float> _b_towerjet2_ihcal_phi;
  std::vector<float> _b_towerjet2_ihcal_dR;
  int _b_towerjet2_ohcal_n; //ohcal towers
  std::vector<float> _b_towerjet2_ohcal_e;
  std::vector<float> _b_towerjet2_ohcal_p;
  std::vector<float> _b_towerjet2_ohcal_pt;
  std::vector<float> _b_towerjet2_ohcal_eta;
  std::vector<float> _b_towerjet2_ohcal_phi;
  std::vector<float> _b_towerjet2_ohcal_dR;
  
  float _b_towerjet3_e;
  float _b_towerjet3_p;
  float _b_towerjet3_pt;
  float _b_towerjet3_eta;
  float _b_towerjet3_phi;
  int _b_towerjet3_cemc_n; //cemc towers
  std::vector<float> _b_towerjet3_cemc_e;
  std::vector<float> _b_towerjet3_cemc_p;
  std::vector<float> _b_towerjet3_cemc_pt;
  std::vector<float> _b_towerjet3_cemc_eta;
  std::vector<float> _b_towerjet3_cemc_phi;
  std::vector<float> _b_towerjet3_cemc_dR;
  int _b_towerjet3_ihcal_n; //ihcal towers
  std::vector<float> _b_towerjet3_ihcal_e;
  std::vector<float> _b_towerjet3_ihcal_p;
  std::vector<float> _b_towerjet3_ihcal_pt;
  std::vector<float> _b_towerjet3_ihcal_eta;
  std::vector<float> _b_towerjet3_ihcal_phi;
  std::vector<float> _b_towerjet3_ihcal_dR;
  int _b_towerjet3_ohcal_n; //ohcal towers
  std::vector<float> _b_towerjet3_ohcal_e;
  std::vector<float> _b_towerjet3_ohcal_p;
  std::vector<float> _b_towerjet3_ohcal_pt;
  std::vector<float> _b_towerjet3_ohcal_eta;
  std::vector<float> _b_towerjet3_ohcal_phi;
  std::vector<float> _b_towerjet3_ohcal_dR;
  
  float _b_towerjet4_e;
  float _b_towerjet4_p;
  float _b_towerjet4_pt;
  float _b_towerjet4_eta;
  float _b_towerjet4_phi;
  int _b_towerjet4_cemc_n; //cemc towers
  std::vector<float> _b_towerjet4_cemc_e;
  std::vector<float> _b_towerjet4_cemc_p;
  std::vector<float> _b_towerjet4_cemc_pt;
  std::vector<float> _b_towerjet4_cemc_eta;
  std::vector<float> _b_towerjet4_cemc_phi;
  std::vector<float> _b_towerjet4_cemc_dR;
  std::vector<int> _b_towerjet4_cemc_id;
  int _b_towerjet4_clcemc_n; //cemc clusters
  std::vector<float> _b_towerjet4_clcemc_e;
  std::vector<float> _b_towerjet4_clcemc_p;
  std::vector<float> _b_towerjet4_clcemc_pt;
  std::vector<float> _b_towerjet4_clcemc_eta;
  std::vector<float> _b_towerjet4_clcemc_phi;
  std::vector<float> _b_towerjet4_clcemc_dR;
  std::vector<int> _b_towerjet4_clcemc_id;
  std::vector<int> _b_towerjet4_clcemc_ntwr;
  std::vector<float> _b_towerjet4_clcemc_prob;
  std::vector<float> _b_towerjet4_clcemc_chi2;
  std::vector<float> _b_towerjet4_clcemc_twrEsum;
  std::vector<float> _b_towerjet4_clcemc_twrNfrac;
  int _b_towerjet4_clcemcMatched_n; //cemc clusters matched to track!
  std::vector<float> _b_towerjet4_clcemcMatched_e;
  std::vector<float> _b_towerjet4_clcemcMatched_p;
  std::vector<float> _b_towerjet4_clcemcMatched_pt;
  std::vector<float> _b_towerjet4_clcemcMatched_eta;
  std::vector<float> _b_towerjet4_clcemcMatched_phi;
  std::vector<float> _b_towerjet4_clcemcMatched_dR;
  std::vector<int> _b_towerjet4_clcemcMatched_id;
  std::vector<int> _b_towerjet4_clcemcMatched_ntwr;
  std::vector<float> _b_towerjet4_clcemcMatched_prob;
  std::vector<float> _b_towerjet4_clcemcMatched_chi2;
  std::vector<float> _b_towerjet4_clcemcMatched_track_p;
  std::vector<float> _b_towerjet4_clcemcMatched_track_deta;
  std::vector<float> _b_towerjet4_clcemcMatched_track_dphi;
  std::vector<float> _b_towerjet4_clcemcMatched_track_dr;
  int _b_towerjet4_ihcal_n; //ihcal towers
  std::vector<float> _b_towerjet4_ihcal_e;
  std::vector<float> _b_towerjet4_ihcal_p;
  std::vector<float> _b_towerjet4_ihcal_pt;
  std::vector<float> _b_towerjet4_ihcal_eta;
  std::vector<float> _b_towerjet4_ihcal_phi;
  std::vector<float> _b_towerjet4_ihcal_dR;
  std::vector<int> _b_towerjet4_ihcal_id;
  int _b_towerjet4_clihcal_n; //ihcal clusters
  std::vector<float> _b_towerjet4_clihcal_e;
  std::vector<float> _b_towerjet4_clihcal_p;
  std::vector<float> _b_towerjet4_clihcal_pt;
  std::vector<float> _b_towerjet4_clihcal_eta;
  std::vector<float> _b_towerjet4_clihcal_phi;
  std::vector<float> _b_towerjet4_clihcal_dR;
  std::vector<int> _b_towerjet4_clihcal_id;
  std::vector<int> _b_towerjet4_clihcal_ntwr;
  std::vector<float> _b_towerjet4_clihcal_twrEsum;
  std::vector<float> _b_towerjet4_clihcal_twrNfrac;
  int _b_towerjet4_ohcal_n; //ohcal towers
  std::vector<float> _b_towerjet4_ohcal_e;
  std::vector<float> _b_towerjet4_ohcal_p;
  std::vector<float> _b_towerjet4_ohcal_pt;
  std::vector<float> _b_towerjet4_ohcal_eta;
  std::vector<float> _b_towerjet4_ohcal_phi;
  std::vector<float> _b_towerjet4_ohcal_dR;
  std::vector<int> _b_towerjet4_ohcal_id;
  int _b_towerjet4_clohcal_n; //ohcal clusters
  std::vector<float> _b_towerjet4_clohcal_e;
  std::vector<float> _b_towerjet4_clohcal_p;
  std::vector<float> _b_towerjet4_clohcal_pt;
  std::vector<float> _b_towerjet4_clohcal_eta;
  std::vector<float> _b_towerjet4_clohcal_phi;
  std::vector<float> _b_towerjet4_clohcal_dR;
  std::vector<int> _b_towerjet4_clohcal_id;
  std::vector<int> _b_towerjet4_clohcal_ntwr;
  std::vector<float> _b_towerjet4_clohcal_twrEsum;
  std::vector<float> _b_towerjet4_clohcal_twrNfrac;
  
  float _b_towerjet5_e;
  float _b_towerjet5_p;
  float _b_towerjet5_pt;
  float _b_towerjet5_eta;
  float _b_towerjet5_phi;
  int _b_towerjet5_cemc_n; //cemc towers
  std::vector<float> _b_towerjet5_cemc_e;
  std::vector<float> _b_towerjet5_cemc_p;
  std::vector<float> _b_towerjet5_cemc_pt;
  std::vector<float> _b_towerjet5_cemc_eta;
  std::vector<float> _b_towerjet5_cemc_phi;
  std::vector<float> _b_towerjet5_cemc_dR;
  int _b_towerjet5_ihcal_n; //ihcal towers
  std::vector<float> _b_towerjet5_ihcal_e;
  std::vector<float> _b_towerjet5_ihcal_p;
  std::vector<float> _b_towerjet5_ihcal_pt;
  std::vector<float> _b_towerjet5_ihcal_eta;
  std::vector<float> _b_towerjet5_ihcal_phi;
  std::vector<float> _b_towerjet5_ihcal_dR;
  int _b_towerjet5_ohcal_n; //ohcal towers
  std::vector<float> _b_towerjet5_ohcal_e;
  std::vector<float> _b_towerjet5_ohcal_p;
  std::vector<float> _b_towerjet5_ohcal_pt;
  std::vector<float> _b_towerjet5_ohcal_eta;
  std::vector<float> _b_towerjet5_ohcal_phi;
  std::vector<float> _b_towerjet5_ohcal_dR;
  
  ///////////////////////////////////////////////////////////////////////// 
  int _b_bhhit_n;
  float _b_bhhit_edepTot;
  std::vector<float> _b_bhhit_edep;
  std::vector<float> _b_bhhit_p;
  std::vector<float> _b_bhhit_pt;
  std::vector<float> _b_bhhit_eta;
  std::vector<float> _b_bhhit_phi;

};

#endif // __TRUTHJETTRIGGER_H__
