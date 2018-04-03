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

//  KyoTreeMaker(const std::string &name, bool fill_g4particle_tree=true, bool fill_track_tree=true, bool fill_tower_tree=true, bool fill_truthjet4_tree=true, bool fill_towerjet4_tree=true, bool fill_bhhit_tree=true);
  KyoTreeMaker(const std::string &name, bool fill_g4particle_tree=true, bool fill_track_tree=true, bool fill_tower_tree=true, bool fill_cluster_tree=true, bool fill_truthjet4_tree=true, bool fill_towerjet4_tree=true, bool fill_otherjetR_tree=true, bool fill_bhhit_tree=true);

  int Init(PHCompositeNode*);
  int process_event(PHCompositeNode*);
  int End(PHCompositeNode*);

 private:

  float PositionCorrectedNxNTowerE( PHCompositeNode *topNode, unsigned int cl_id, float eNN);
  void getClusterByIndex( PHCompositeNode *topNode, std::string detName, unsigned int cl_id);
  void getClusterPCByIndex( PHCompositeNode *topNode, unsigned int cl_id);

  void GetJetPrimaryContributors( PHCompositeNode *topNode, PHG4Particle* g4ptl);
  void TraceHepMCParticle( HepMC::GenParticle *ptcle);

  TFile *_f;
  std::string _foutname;

  bool _fill_g4particle_tree;
  bool _fill_track_tree;
  bool _fill_tower_tree;
  bool _fill_cluster_tree;
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
  TTree *cluster_tree;
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
  float _b_vtx_x;
  float _b_vtx_y;
  float _b_vtx_z;
  //float _b_vtx_t;
  
  ///////////////////////////////////////////////////////////////////////// 
  int _b_g4particle_n;
  float _b_g4particle_e[1000];
  float _b_g4particle_p[1000];
  float _b_g4particle_pt[1000];
  float _b_g4particle_eta[1000];
  float _b_g4particle_phi[1000];
  int _b_g4particle_pid[1000];
  
  ///////////////////////////////////////////////////////////////////////// 
  int _b_track_n;
  float _b_track_e[1000];
  float _b_track_p[1000];
  float _b_track_pt[1000];
  float _b_track_eta[1000];
  float _b_track_phi[1000];
  int _b_track_charge[1000];
  float _b_track_chisq[1000];
  int _b_track_ndf[1000];
  float _b_track_dca2d[1000];
  float _b_cemc_E3x3[1000];
  float _b_cemc_E3x3_PC[1000];
  float _b_ihcal_E3x3[1000];
  float _b_ohcal_E3x3[1000];
  float _b_cemc_E5x5[1000];
  float _b_cemc_E5x5_PC[1000];
  float _b_ihcal_E5x5[1000];
  float _b_ohcal_E5x5[1000];
  int _b_cemc_clN[1000];
  float _b_cemc_clDeta[1000];
  float _b_cemc_clDphi[1000];
  float _b_cemc_clE[1000];
  float _b_cemc_clE_PC[1000];
  int _b_ihcal_clN[1000];
  float _b_ihcal_clDeta[1000];
  float _b_ihcal_clDphi[1000];
  float _b_ihcal_clE[1000];
  int _b_ohcal_clN[1000];
  float _b_ohcal_clDeta[1000];
  float _b_ohcal_clDphi[1000];
  float _b_ohcal_clE[1000];
  
  ///////////////////////////////////////////////////////////////////////// 
  int _b_cemc_n;
  float _b_cemc_e[1000];
  float _b_cemc_p[1000];
  float _b_cemc_pt[1000];
  float _b_cemc_eta[1000];
  float _b_cemc_phi[1000];
  int _b_ihcal_n;
  float _b_ihcal_e[1000];
  float _b_ihcal_p[1000];
  float _b_ihcal_pt[1000];
  float _b_ihcal_eta[1000];
  float _b_ihcal_phi[1000];
  int _b_ohcal_n;
  float _b_ohcal_e[1000];
  float _b_ohcal_p[1000];
  float _b_ohcal_pt[1000];
  float _b_ohcal_eta[1000];
  float _b_ohcal_phi[1000];
  
  ///////////////////////////////////////////////////////////////////////// 
  int _b_cl_cemc_n;
  float _b_cl_cemc_e[1000];
  float _b_cl_cemc_p[1000];
  float _b_cl_cemc_pt[1000];
  float _b_cl_cemc_eta[1000];
  float _b_cl_cemc_phi[1000];
  int _b_cl_cemc_ntwr[1000];
  float _b_cl_cemc_prob[1000];
  float _b_cl_cemc_chi2[1000];
  int _b_cl_ihcal_n;
  float _b_cl_ihcal_e[1000];
  float _b_cl_ihcal_p[1000];
  float _b_cl_ihcal_pt[1000];
  float _b_cl_ihcal_eta[1000];
  float _b_cl_ihcal_phi[1000];
  int _b_cl_ihcal_ntwr[1000];
  //float _b_cl_ihcal_prob[1000];
  //float _b_cl_ihcal_chi2[1000];
  int _b_cl_ohcal_n;
  float _b_cl_ohcal_e[1000];
  float _b_cl_ohcal_p[1000];
  float _b_cl_ohcal_pt[1000];
  float _b_cl_ohcal_eta[1000];
  float _b_cl_ohcal_phi[1000];
  int _b_cl_ohcal_ntwr[1000];
  //float _b_cl_ohcal_prob[1000];
  //float _b_cl_ohcal_chi2[1000];

  ///////////////////////////////////////////////////////////////////////// 
  float _b_truthjet2_e;
  float _b_truthjet2_p;
  float _b_truthjet2_pt;
  float _b_truthjet2_eta;
  float _b_truthjet2_phi;
  int _b_truthjet2_cons_n; //constituents
  float _b_truthjet2_cons_e[200];
  float _b_truthjet2_cons_p[200];
  float _b_truthjet2_cons_pt[200];
  float _b_truthjet2_cons_eta[200];
  float _b_truthjet2_cons_phi[200];
  float _b_truthjet2_cons_dR[200];
  int _b_truthjet2_cons_pid[200];
  int _b_truthjet2_pg_n; //progenitors
  int _b_truthjet2_pg_id[1000];
  float _b_truthjet2_pg_fract[1000];
  int _b_truthjet2_pg_status[1000];
  
  float _b_truthjet3_e;
  float _b_truthjet3_p;
  float _b_truthjet3_pt;
  float _b_truthjet3_eta;
  float _b_truthjet3_phi;
  int _b_truthjet3_cons_n; //constituents
  float _b_truthjet3_cons_e[200];
  float _b_truthjet3_cons_p[200];
  float _b_truthjet3_cons_pt[200];
  float _b_truthjet3_cons_eta[200];
  float _b_truthjet3_cons_phi[200];
  float _b_truthjet3_cons_dR[200];
  int _b_truthjet3_cons_pid[200];
  int _b_truthjet3_pg_n; //progenitors
  int _b_truthjet3_pg_id[1000];
  float _b_truthjet3_pg_fract[1000];
  int _b_truthjet3_pg_status[1000];
  
  float _b_truthjet4_e;
  float _b_truthjet4_p;
  float _b_truthjet4_pt;
  float _b_truthjet4_eta;
  float _b_truthjet4_phi;
  int _b_truthjet4_cons_n; //constituents
  float _b_truthjet4_cons_e[200];
  float _b_truthjet4_cons_p[200];
  float _b_truthjet4_cons_pt[200];
  float _b_truthjet4_cons_eta[200];
  float _b_truthjet4_cons_phi[200];
  float _b_truthjet4_cons_dR[200];
  int _b_truthjet4_cons_pid[200];
  int _b_truthjet4_pg_n; //progenitors
  int _b_truthjet4_pg_id[1000];
  float _b_truthjet4_pg_fract[1000];
  int _b_truthjet4_pg_status[1000];
  
  float _b_truthjet5_e;
  float _b_truthjet5_p;
  float _b_truthjet5_pt;
  float _b_truthjet5_eta;
  float _b_truthjet5_phi;
  int _b_truthjet5_cons_n; //constituents
  float _b_truthjet5_cons_e[200];
  float _b_truthjet5_cons_p[200];
  float _b_truthjet5_cons_pt[200];
  float _b_truthjet5_cons_eta[200];
  float _b_truthjet5_cons_phi[200];
  float _b_truthjet5_cons_dR[200];
  int _b_truthjet5_cons_pid[200];
  int _b_truthjet5_pg_n; //progenitors
  int _b_truthjet5_pg_id[1000];
  float _b_truthjet5_pg_fract[1000];
  int _b_truthjet5_pg_status[1000];

  ///////////////////////////////////////////////////////////////////////// 
  float _b_towerjet2_e;
  float _b_towerjet2_p;
  float _b_towerjet2_pt;
  float _b_towerjet2_eta;
  float _b_towerjet2_phi;
  int _b_towerjet2_cemc_n; //cemc towers
//  float _b_towerjet2_cemc_scale;
  float _b_towerjet2_cemc_e[200];
  float _b_towerjet2_cemc_p[200];
  float _b_towerjet2_cemc_pt[200];
  float _b_towerjet2_cemc_eta[200];
  float _b_towerjet2_cemc_phi[200];
  float _b_towerjet2_cemc_dR[200];
  int _b_towerjet2_ihcal_n; //ihcal towers
//  float _b_towerjet2_ihcal_scale;
  float _b_towerjet2_ihcal_e[200];
  float _b_towerjet2_ihcal_p[200];
  float _b_towerjet2_ihcal_pt[200];
  float _b_towerjet2_ihcal_eta[200];
  float _b_towerjet2_ihcal_phi[200];
  float _b_towerjet2_ihcal_dR[200];
  int _b_towerjet2_ohcal_n; //ohcal towers
//  float _b_towerjet2_ohcal_scale;
  float _b_towerjet2_ohcal_e[200];
  float _b_towerjet2_ohcal_p[200];
  float _b_towerjet2_ohcal_pt[200];
  float _b_towerjet2_ohcal_eta[200];
  float _b_towerjet2_ohcal_phi[200];
  float _b_towerjet2_ohcal_dR[200];
  
  float _b_towerjet3_e;
  float _b_towerjet3_p;
  float _b_towerjet3_pt;
  float _b_towerjet3_eta;
  float _b_towerjet3_phi;
  int _b_towerjet3_cemc_n; //cemc towers
//  float _b_towerjet3_cemc_scale;
  float _b_towerjet3_cemc_e[200];
  float _b_towerjet3_cemc_p[200];
  float _b_towerjet3_cemc_pt[200];
  float _b_towerjet3_cemc_eta[200];
  float _b_towerjet3_cemc_phi[200];
  float _b_towerjet3_cemc_dR[200];
  int _b_towerjet3_ihcal_n; //ihcal towers
//  float _b_towerjet3_ihcal_scale;
  float _b_towerjet3_ihcal_e[200];
  float _b_towerjet3_ihcal_p[200];
  float _b_towerjet3_ihcal_pt[200];
  float _b_towerjet3_ihcal_eta[200];
  float _b_towerjet3_ihcal_phi[200];
  float _b_towerjet3_ihcal_dR[200];
  int _b_towerjet3_ohcal_n; //ohcal towers
//  float _b_towerjet3_ohcal_scale;
  float _b_towerjet3_ohcal_e[200];
  float _b_towerjet3_ohcal_p[200];
  float _b_towerjet3_ohcal_pt[200];
  float _b_towerjet3_ohcal_eta[200];
  float _b_towerjet3_ohcal_phi[200];
  float _b_towerjet3_ohcal_dR[200];
  
  float _b_towerjet4_e;
  float _b_towerjet4_p;
  float _b_towerjet4_pt;
  float _b_towerjet4_eta;
  float _b_towerjet4_phi;
  int _b_towerjet4_cemc_n; //cemc towers
//  float _b_towerjet4_cemc_scale;
  float _b_towerjet4_cemc_e[200];
  float _b_towerjet4_cemc_p[200];
  float _b_towerjet4_cemc_pt[200];
  float _b_towerjet4_cemc_eta[200];
  float _b_towerjet4_cemc_phi[200];
  float _b_towerjet4_cemc_dR[200];
  int _b_towerjet4_ihcal_n; //ihcal towers
//  float _b_towerjet4_ihcal_scale;
  float _b_towerjet4_ihcal_e[200];
  float _b_towerjet4_ihcal_p[200];
  float _b_towerjet4_ihcal_pt[200];
  float _b_towerjet4_ihcal_eta[200];
  float _b_towerjet4_ihcal_phi[200];
  float _b_towerjet4_ihcal_dR[200];
  int _b_towerjet4_ohcal_n; //ohcal towers
//  float _b_towerjet4_ohcal_scale;
  float _b_towerjet4_ohcal_e[200];
  float _b_towerjet4_ohcal_p[200];
  float _b_towerjet4_ohcal_pt[200];
  float _b_towerjet4_ohcal_eta[200];
  float _b_towerjet4_ohcal_phi[200];
  float _b_towerjet4_ohcal_dR[200];
  
  float _b_towerjet5_e;
  float _b_towerjet5_p;
  float _b_towerjet5_pt;
  float _b_towerjet5_eta;
  float _b_towerjet5_phi;
  int _b_towerjet5_cemc_n; //cemc towers
//  float _b_towerjet5_cemc_scale;
  float _b_towerjet5_cemc_e[200];
  float _b_towerjet5_cemc_p[200];
  float _b_towerjet5_cemc_pt[200];
  float _b_towerjet5_cemc_eta[200];
  float _b_towerjet5_cemc_phi[200];
  float _b_towerjet5_cemc_dR[200];
  int _b_towerjet5_ihcal_n; //ihcal towers
//  float _b_towerjet5_ihcal_scale;
  float _b_towerjet5_ihcal_e[200];
  float _b_towerjet5_ihcal_p[200];
  float _b_towerjet5_ihcal_pt[200];
  float _b_towerjet5_ihcal_eta[200];
  float _b_towerjet5_ihcal_phi[200];
  float _b_towerjet5_ihcal_dR[200];
  int _b_towerjet5_ohcal_n; //ohcal towers
//  float _b_towerjet5_ohcal_scale;
  float _b_towerjet5_ohcal_e[200];
  float _b_towerjet5_ohcal_p[200];
  float _b_towerjet5_ohcal_pt[200];
  float _b_towerjet5_ohcal_eta[200];
  float _b_towerjet5_ohcal_phi[200];
  float _b_towerjet5_ohcal_dR[200];
  
  ///////////////////////////////////////////////////////////////////////// 
  int _b_bhhit_n;
  float _b_bhhit_edepTot;
  float _b_bhhit_edep[1000];
  float _b_bhhit_p[1000];
  float _b_bhhit_pt[1000];
  float _b_bhhit_eta[1000];
  float _b_bhhit_phi[1000];

};

#endif // __TRUTHJETTRIGGER_H__
