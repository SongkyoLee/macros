#include "KyoTreeMaker.h"

#include <phool/getClass.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include "TLorentzVector.h"
#include <iostream>

#include <calobase/RawClusterContainer.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterUtility.h>

#include <calobase/RawTowerGeom.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
//#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomContainer.h>

#include <g4jets/JetMap.h>
#include <g4jets/Jet.h>

//#include <g4jets/JetUtility.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

#include <g4main/PHG4VtxPoint.h>
#include <g4vertex/GlobalVertexMap.h>
#include <g4vertex/GlobalVertex.h>
//#include <jetbackground/TowerBackground.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <HepMC/GenEvent.h>
#include <HepMC/GenParticle.h>

#include <g4hough/SvtxTrackMap.h>
#include <g4hough/SvtxTrack.h>

#include <cmath>
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// Progenitor data structures
typedef struct {
  int pid; 
  int status; 
  int barcode; 
  float energy;     // constituent energy assoc. with progenitor
  float pg_energy;  // progenitor energy 
} JetProgenitor; 

std::vector<JetProgenitor> progenitors; //for a constituent
std::vector<JetProgenitor> combined_progenitors; //for a jet

bool sortfunction( JetProgenitor i, JetProgenitor j){
  return (i.energy > j.energy);
}

KyoTreeMaker::KyoTreeMaker(const std::string &name, bool fill_g4particle_tree, bool fill_track_tree, bool fill_tower_tree, bool fill_cluster_tree, bool fill_truthjet4_tree, bool fill_towerjet4_tree, bool fill_bhhit_tree) : SubsysReco()
{
  //initialize
  _b_event = 0;
  _foutname = name;
  
  _fill_g4particle_tree=fill_g4particle_tree;
  _fill_track_tree=fill_track_tree;
  _fill_tower_tree=fill_tower_tree;
  //
  _fill_cemc_cluster_tree=fill_cluster_tree;
  _fill_ihcal_cluster_tree=fill_cluster_tree;
  _fill_ohcal_cluster_tree=fill_cluster_tree;
  //
  _fill_truthjet4_tree=fill_truthjet4_tree; 
  _fill_towerjet4_tree=fill_towerjet4_tree;
  //
  _fill_bhhit_tree=fill_bhhit_tree;
}

////////////////////////////////////////////////////////////////////////////
int KyoTreeMaker::Init(PHCompositeNode *topNode)
{

  _f = new TFile( _foutname.c_str(), "RECREATE");

  if(_fill_g4particle_tree) g4particle_tree = new TTree("g4particle_tree","all (true) g4particle");
  if (_fill_track_tree) track_tree = new TTree("track_tree","all reco track");
  if (_fill_tower_tree) tower_tree = new TTree("tower_tree","all raw towers");
  if (_fill_cemc_cluster_tree) cemc_cluster_tree = new TTree("cemc_cluster_tree","all raw cemc_clusters");
  if (_fill_ihcal_cluster_tree) ihcal_cluster_tree = new TTree("ihcal_cluster_tree","all raw ihcal_clusters");
  if (_fill_ohcal_cluster_tree) ohcal_cluster_tree = new TTree("ohcal_cluster_tree","all raw ohcal_clusters");
  if (_fill_truthjet4_tree) truthjet4_tree = new TTree("truthjet4_tree","truth jet anti-Kt, R=0.4");
  if (_fill_towerjet4_tree) towerjet4_tree = new TTree("towerjet4_tree","tower jet anti-Kt, R=0.4");
  if (_fill_bhhit_tree) bhhit_tree = new TTree("bhhit_tree","hits in black hole");

  if (_fill_g4particle_tree){
    g4particle_tree->Branch("event", &_b_event, "event/I");
    g4particle_tree->Branch("gvtx_x", &_b_gvtx_x, "gvtx_x/F");
    g4particle_tree->Branch("gvtx_y", &_b_gvtx_y, "gvtx_y/F");
    g4particle_tree->Branch("gvtx_z", &_b_gvtx_z, "gvtx_z/F");
    g4particle_tree->Branch("vtx_x", &_b_vtx_x, "vtx_x/F");
    g4particle_tree->Branch("vtx_y", &_b_vtx_y, "vtx_y/F");
    g4particle_tree->Branch("vtx_z", &_b_vtx_z, "vtx_z/F");
    g4particle_tree->Branch("g4particle_n", &_b_g4particle_n, "g4particle_n/I");
    g4particle_tree->Branch("g4particle_e", &_b_g4particle_e);
    g4particle_tree->Branch("g4particle_p", &_b_g4particle_p);
    g4particle_tree->Branch("g4particle_pt", &_b_g4particle_pt);
    g4particle_tree->Branch("g4particle_eta", &_b_g4particle_eta);
    g4particle_tree->Branch("g4particle_phi", &_b_g4particle_phi);
    g4particle_tree->Branch("g4particle_pid", &_b_g4particle_pid);

  }
  if (_fill_track_tree){
    track_tree->Branch("event", &_b_event, "event/I");
    track_tree->Branch("track_n", &_b_track_n,"track_n/I");
    track_tree->Branch("track_p", &_b_track_p);
    track_tree->Branch("track_pt", &_b_track_pt);
    track_tree->Branch("track_eta", &_b_track_eta);
    track_tree->Branch("track_phi", &_b_track_phi);
    track_tree->Branch("track_charge", &_b_track_charge);
    track_tree->Branch("track_chisq", &_b_track_chisq);
    track_tree->Branch("track_ndf", &_b_track_ndf);
    track_tree->Branch("track_dca2d", &_b_track_dca2d);
    track_tree->Branch("cemc_E3x3", &_b_cemc_E3x3);
    track_tree->Branch("cemc_E3x3_PC", &_b_cemc_E3x3_PC);
    track_tree->Branch("ihcal_E3x3", &_b_ihcal_E3x3);
    track_tree->Branch("ohcal_E3x3", &_b_ohcal_E3x3);
    track_tree->Branch("cemc_E5x5", &_b_cemc_E5x5);
    track_tree->Branch("cemc_E5x5_PC", &_b_cemc_E5x5_PC);
    track_tree->Branch("ihcal_E5x5", &_b_ihcal_E5x5);
    track_tree->Branch("ohcal_E5x5", &_b_ohcal_E5x5);
    track_tree->Branch("cemc_clNtwr", &_b_cemc_clNtwr);
    track_tree->Branch("cemc_clE", &_b_cemc_clE);
    track_tree->Branch("cemc_cldeta", &_b_cemc_cldeta);
    track_tree->Branch("cemc_cldphi", &_b_cemc_cldphi);
    track_tree->Branch("cemc_clE_PC", &_b_cemc_clE_PC);
    track_tree->Branch("ihcal_clNtwr", &_b_ihcal_clNtwr);
    track_tree->Branch("ihcal_clE", &_b_ihcal_clE);
    track_tree->Branch("ihcal_cldeta", &_b_ihcal_cldeta);
    track_tree->Branch("ihcal_cldphi", &_b_ihcal_cldphi);
    track_tree->Branch("ohcal_clNtwr", &_b_ohcal_clNtwr);
    track_tree->Branch("ohcal_clE", &_b_ohcal_clE);
    track_tree->Branch("ohcal_cldeta", &_b_ohcal_cldeta);
    track_tree->Branch("ohcal_cldphi", &_b_ohcal_cldphi);
  }  
  if (_fill_tower_tree){
    tower_tree->Branch("event", &_b_event, "event/I");
    tower_tree->Branch("cemc_n", &_b_cemc_n, "cemc_n/I");
    tower_tree->Branch("cemc_e", &_b_cemc_e);
    tower_tree->Branch("cemc_p", &_b_cemc_p);
    tower_tree->Branch("cemc_pt", &_b_cemc_pt);
    tower_tree->Branch("cemc_eta", &_b_cemc_eta);
    tower_tree->Branch("cemc_phi", &_b_cemc_phi);
    tower_tree->Branch("cemc_id", &_b_cemc_id);
    tower_tree->Branch("ihcal_n", &_b_ihcal_n, "ihcal_n/I");
    tower_tree->Branch("ihcal_e", &_b_ihcal_e);
    tower_tree->Branch("ihcal_p", &_b_ihcal_p);
    tower_tree->Branch("ihcal_pt", &_b_ihcal_pt);
    tower_tree->Branch("ihcal_eta", &_b_ihcal_eta);
    tower_tree->Branch("ihcal_phi", &_b_ihcal_phi);
    tower_tree->Branch("ohcal_n", &_b_ohcal_n, "ohcal_n/I");
    tower_tree->Branch("ohcal_e", &_b_ohcal_e);
    tower_tree->Branch("ohcal_p", &_b_ohcal_p);
    tower_tree->Branch("ohcal_pt", &_b_ohcal_pt);
    tower_tree->Branch("ohcal_eta", &_b_ohcal_eta);
    tower_tree->Branch("ohcal_phi", &_b_ohcal_phi);
  }
  if (_fill_cemc_cluster_tree){
    cemc_cluster_tree->Branch("event", &_b_event, "event/I");
    cemc_cluster_tree->Branch("cl_cemc_e", &_b_cl_cemc_e, "cl_cemc_e/F");
    cemc_cluster_tree->Branch("cl_cemc_p", &_b_cl_cemc_p, "cl_cemc_p/F");
    cemc_cluster_tree->Branch("cl_cemc_pt", &_b_cl_cemc_pt, "cl_cemc_pt/F");
    cemc_cluster_tree->Branch("cl_cemc_eta", &_b_cl_cemc_eta, "cl_cemc_eta/F");
    cemc_cluster_tree->Branch("cl_cemc_phi", &_b_cl_cemc_phi, "cl_cemc_phi/F");
    cemc_cluster_tree->Branch("cl_cemc_id", &_b_cl_cemc_id, "cl_cemc_id/I");
    cemc_cluster_tree->Branch("cl_cemc_ntwr", &_b_cl_cemc_ntwr, "cl_cemc_ntwr/I");
    cemc_cluster_tree->Branch("cl_cemc_twrEsum", &_b_cl_cemc_twrEsum, "cl_cemc_twrEsum/F");
    cemc_cluster_tree->Branch("cl_cemc_prob", &_b_cl_cemc_prob, "cl_cemc_prob/F"); //only for CEMC
    cemc_cluster_tree->Branch("cl_cemc_chi2", &_b_cl_cemc_chi2, "cl_cemc_chi2/F");
    cemc_cluster_tree->Branch("cl_twr_cemc_e", &_b_cl_twr_cemc_e);
    cemc_cluster_tree->Branch("cl_twr_cemc_id", &_b_cl_twr_cemc_id);
  }
  if (_fill_ihcal_cluster_tree){
    ihcal_cluster_tree->Branch("event", &_b_event, "event/I");
    ihcal_cluster_tree->Branch("cl_ihcal_e", &_b_cl_ihcal_e, "cl_ihcal_e/F");
    ihcal_cluster_tree->Branch("cl_ihcal_p", &_b_cl_ihcal_p, "cl_ihcal_p/F");
    ihcal_cluster_tree->Branch("cl_ihcal_pt", &_b_cl_ihcal_pt, "cl_ihcal_pt/F");
    ihcal_cluster_tree->Branch("cl_ihcal_eta", &_b_cl_ihcal_eta, "cl_ihcal_eta/F");
    ihcal_cluster_tree->Branch("cl_ihcal_phi", &_b_cl_ihcal_phi, "cl_ihcal_phi/F");
    ihcal_cluster_tree->Branch("cl_ihcal_id", &_b_cl_ihcal_id, "cl_ihcal_id/I");
    ihcal_cluster_tree->Branch("cl_ihcal_ntwr", &_b_cl_ihcal_ntwr, "cl_ihcal_ntwr/I");
    ihcal_cluster_tree->Branch("cl_twr_ihcal_e", &_b_cl_twr_ihcal_e);
    ihcal_cluster_tree->Branch("cl_twr_ihcal_id", &_b_cl_twr_ihcal_id);
  }
  if (_fill_ohcal_cluster_tree){
    ohcal_cluster_tree->Branch("event", &_b_event, "event/I");
    ohcal_cluster_tree->Branch("cl_ohcal_e", &_b_cl_ohcal_e, "cl_ohcal_e/F");
    ohcal_cluster_tree->Branch("cl_ohcal_p", &_b_cl_ohcal_p, "cl_ohcal_p/F");
    ohcal_cluster_tree->Branch("cl_ohcal_pt", &_b_cl_ohcal_pt, "cl_ohcal_pt/F");
    ohcal_cluster_tree->Branch("cl_ohcal_eta", &_b_cl_ohcal_eta, "cl_ohcal_eta/F");
    ohcal_cluster_tree->Branch("cl_ohcal_phi", &_b_cl_ohcal_phi, "cl_ohcal_phi/F");
    ohcal_cluster_tree->Branch("cl_ohcal_id", &_b_cl_ohcal_id, "cl_ohcal_id/I");
    ohcal_cluster_tree->Branch("cl_ohcal_ntwr", &_b_cl_ohcal_ntwr, "cl_ohcal_ntwr/I");
    ohcal_cluster_tree->Branch("cl_twr_ohcal_e", &_b_cl_twr_ohcal_e);
    ohcal_cluster_tree->Branch("cl_twr_ohcal_id", &_b_cl_twr_ohcal_id);
  }
  if (_fill_truthjet4_tree){
    truthjet4_tree->Branch("event", &_b_event, "event/I");
    truthjet4_tree->Branch("gvtx_x", &_b_gvtx_x, "gvtx_x/F");
    truthjet4_tree->Branch("gvtx_y", &_b_gvtx_y, "gvtx_y/F");
    truthjet4_tree->Branch("gvtx_z", &_b_gvtx_z, "gvtx_z/F");
    truthjet4_tree->Branch("vtx_x", &_b_vtx_x, "vtx_x/F");
    truthjet4_tree->Branch("vtx_y", &_b_vtx_y, "vtx_y/F");
    truthjet4_tree->Branch("vtx_z", &_b_vtx_z, "vtx_z/F");
    truthjet4_tree->Branch("truthjet4_e", &_b_truthjet4_e, "truthjet4_e/F");
    truthjet4_tree->Branch("truthjet4_p", &_b_truthjet4_p, "truthjet4_p/F");
    truthjet4_tree->Branch("truthjet4_pt", &_b_truthjet4_pt, "truthjet4_pt/F");
    truthjet4_tree->Branch("truthjet4_eta", &_b_truthjet4_eta, "truthjet4_eta/F");
    truthjet4_tree->Branch("truthjet4_phi", &_b_truthjet4_phi, "truthjet4_phi/F");
    truthjet4_tree->Branch("truthjet4_cons_n", &_b_truthjet4_cons_n, "truthjet4_cons_n/I");
    truthjet4_tree->Branch("truthjet4_cons_e", &_b_truthjet4_cons_e);
    truthjet4_tree->Branch("truthjet4_cons_p", &_b_truthjet4_cons_p);
    truthjet4_tree->Branch("truthjet4_cons_pt", &_b_truthjet4_cons_pt);
    truthjet4_tree->Branch("truthjet4_cons_eta", &_b_truthjet4_cons_eta);
    truthjet4_tree->Branch("truthjet4_cons_phi", &_b_truthjet4_cons_phi);
    truthjet4_tree->Branch("truthjet4_cons_dR", &_b_truthjet4_cons_dR);
    truthjet4_tree->Branch("truthjet4_cons_pid", &_b_truthjet4_cons_pid);
    truthjet4_tree->Branch("truthjet4_pg_n", &_b_truthjet4_pg_n, "truthjet4_pg_n/I");
    truthjet4_tree->Branch("truthjet4_pg_id", &_b_truthjet4_pg_id);
    truthjet4_tree->Branch("truthjet4_pg_fract", &_b_truthjet4_pg_fract);
    truthjet4_tree->Branch("truthjet4_pg_status", &_b_truthjet4_pg_status);
  }  
  if (_fill_towerjet4_tree){
    towerjet4_tree->Branch("event", &_b_event, "event/I");
    towerjet4_tree->Branch("gvtx_x", &_b_gvtx_x, "gvtx_x/F");
    towerjet4_tree->Branch("gvtx_y", &_b_gvtx_y, "gvtx_y/F");
    towerjet4_tree->Branch("gvtx_z", &_b_gvtx_z, "gvtx_z/F");
    towerjet4_tree->Branch("vtx_x", &_b_vtx_x, "vtx_x/F");
    towerjet4_tree->Branch("vtx_y", &_b_vtx_y, "vtx_y/F");
    towerjet4_tree->Branch("vtx_z", &_b_vtx_z, "vtx_z/F");
    towerjet4_tree->Branch("towerjet4_e", &_b_towerjet4_e, "towerjet4_e/F");
    towerjet4_tree->Branch("towerjet4_p", &_b_towerjet4_p, "towerjet4_p/F");
    towerjet4_tree->Branch("towerjet4_pt", &_b_towerjet4_pt, "towerjet4_pt/F");
    towerjet4_tree->Branch("towerjet4_eta", &_b_towerjet4_eta, "towerjet4_eta/F");
    towerjet4_tree->Branch("towerjet4_phi", &_b_towerjet4_phi, "towerjet4_phi/F");
    towerjet4_tree->Branch("towerjet4_cemc_n", &_b_towerjet4_cemc_n, "towerjet4_cemc_n/I");
    towerjet4_tree->Branch("towerjet4_cemc_e", &_b_towerjet4_cemc_e);
    towerjet4_tree->Branch("towerjet4_cemc_p", &_b_towerjet4_cemc_p);
    towerjet4_tree->Branch("towerjet4_cemc_pt", &_b_towerjet4_cemc_pt);
    towerjet4_tree->Branch("towerjet4_cemc_eta", &_b_towerjet4_cemc_eta);
    towerjet4_tree->Branch("towerjet4_cemc_phi", &_b_towerjet4_cemc_phi);
    towerjet4_tree->Branch("towerjet4_cemc_dR", &_b_towerjet4_cemc_dR);
    towerjet4_tree->Branch("towerjet4_cemc_id", &_b_towerjet4_cemc_id);
    towerjet4_tree->Branch("towerjet4_clcemc_n", &_b_towerjet4_clcemc_n, "towerjet4_clcemc_n/I");
    towerjet4_tree->Branch("towerjet4_clcemc_e", &_b_towerjet4_clcemc_e);
    towerjet4_tree->Branch("towerjet4_clcemc_p", &_b_towerjet4_clcemc_p);
    towerjet4_tree->Branch("towerjet4_clcemc_pt", &_b_towerjet4_clcemc_pt);
    towerjet4_tree->Branch("towerjet4_clcemc_eta", &_b_towerjet4_clcemc_eta);
    towerjet4_tree->Branch("towerjet4_clcemc_phi", &_b_towerjet4_clcemc_phi);
    towerjet4_tree->Branch("towerjet4_clcemc_dR", &_b_towerjet4_clcemc_dR);
    towerjet4_tree->Branch("towerjet4_clcemc_id", &_b_towerjet4_clcemc_id);
    towerjet4_tree->Branch("towerjet4_clcemc_ntwr", &_b_towerjet4_clcemc_ntwr);
    towerjet4_tree->Branch("towerjet4_clcemc_prob", &_b_towerjet4_clcemc_prob);
    towerjet4_tree->Branch("towerjet4_clcemc_chi2", &_b_towerjet4_clcemc_chi2);
    towerjet4_tree->Branch("towerjet4_clcemc_twrEsum", &_b_towerjet4_clcemc_twrEsum);
    towerjet4_tree->Branch("towerjet4_clcemc_twrNfrac", &_b_towerjet4_clcemc_twrNfrac);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_n", &_b_towerjet4_clcemcMatched_n, "towerjet4_clcemcMatched_n/I");
    towerjet4_tree->Branch("towerjet4_clcemcMatched_e", &_b_towerjet4_clcemcMatched_e);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_p", &_b_towerjet4_clcemcMatched_p);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_pt", &_b_towerjet4_clcemcMatched_pt);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_eta", &_b_towerjet4_clcemcMatched_eta);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_phi", &_b_towerjet4_clcemcMatched_phi);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_dR", &_b_towerjet4_clcemcMatched_dR);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_id", &_b_towerjet4_clcemcMatched_id);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_ntwr", &_b_towerjet4_clcemcMatched_ntwr);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_prob", &_b_towerjet4_clcemcMatched_prob);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_chi2", &_b_towerjet4_clcemcMatched_chi2);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_track_p", &_b_towerjet4_clcemcMatched_track_p);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_track_deta", &_b_towerjet4_clcemcMatched_track_deta);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_track_dphi", &_b_towerjet4_clcemcMatched_track_dphi);
    towerjet4_tree->Branch("towerjet4_clcemcMatched_track_dr", &_b_towerjet4_clcemcMatched_track_dr);
    towerjet4_tree->Branch("towerjet4_ihcal_n", &_b_towerjet4_ihcal_n, "towerjet4_ihcal_n/I");
    towerjet4_tree->Branch("towerjet4_ihcal_e", &_b_towerjet4_ihcal_e);
    towerjet4_tree->Branch("towerjet4_ihcal_p", &_b_towerjet4_ihcal_p);
    towerjet4_tree->Branch("towerjet4_ihcal_pt", &_b_towerjet4_ihcal_pt);
    towerjet4_tree->Branch("towerjet4_ihcal_eta", &_b_towerjet4_ihcal_eta);
    towerjet4_tree->Branch("towerjet4_ihcal_phi", &_b_towerjet4_ihcal_phi);
    towerjet4_tree->Branch("towerjet4_ihcal_dR", &_b_towerjet4_ihcal_dR);
    towerjet4_tree->Branch("towerjet4_ihcal_id", &_b_towerjet4_ihcal_id);
    towerjet4_tree->Branch("towerjet4_clihcal_n", &_b_towerjet4_clihcal_n, "towerjet4_clihcal_n/I");
    towerjet4_tree->Branch("towerjet4_clihcal_e", &_b_towerjet4_clihcal_e);
    towerjet4_tree->Branch("towerjet4_clihcal_p", &_b_towerjet4_clihcal_p);
    towerjet4_tree->Branch("towerjet4_clihcal_pt", &_b_towerjet4_clihcal_pt);
    towerjet4_tree->Branch("towerjet4_clihcal_eta", &_b_towerjet4_clihcal_eta);
    towerjet4_tree->Branch("towerjet4_clihcal_phi", &_b_towerjet4_clihcal_phi);
    towerjet4_tree->Branch("towerjet4_clihcal_dR", &_b_towerjet4_clihcal_dR);
    towerjet4_tree->Branch("towerjet4_clihcal_id", &_b_towerjet4_clihcal_id);
    towerjet4_tree->Branch("towerjet4_clihcal_ntwr", &_b_towerjet4_clihcal_ntwr);
    towerjet4_tree->Branch("towerjet4_clihcal_twrEsum", &_b_towerjet4_clihcal_twrEsum);
    towerjet4_tree->Branch("towerjet4_clihcal_twrNfrac", &_b_towerjet4_clihcal_twrNfrac);
    towerjet4_tree->Branch("towerjet4_ohcal_n", &_b_towerjet4_ohcal_n, "towerjet4_ohcal_n/I");
    towerjet4_tree->Branch("towerjet4_ohcal_e", &_b_towerjet4_ohcal_e);
    towerjet4_tree->Branch("towerjet4_ohcal_p", &_b_towerjet4_ohcal_p);
    towerjet4_tree->Branch("towerjet4_ohcal_pt", &_b_towerjet4_ohcal_pt);
    towerjet4_tree->Branch("towerjet4_ohcal_eta", &_b_towerjet4_ohcal_eta);
    towerjet4_tree->Branch("towerjet4_ohcal_phi", &_b_towerjet4_ohcal_phi);
    towerjet4_tree->Branch("towerjet4_ohcal_dR", &_b_towerjet4_ohcal_dR);
    towerjet4_tree->Branch("towerjet4_ohcal_id", &_b_towerjet4_ohcal_id);
    towerjet4_tree->Branch("towerjet4_clohcal_n", &_b_towerjet4_clohcal_n, "towerjet4_clohcal_n/I");
    towerjet4_tree->Branch("towerjet4_clohcal_e", &_b_towerjet4_clohcal_e);
    towerjet4_tree->Branch("towerjet4_clohcal_p", &_b_towerjet4_clohcal_p);
    towerjet4_tree->Branch("towerjet4_clohcal_pt", &_b_towerjet4_clohcal_pt);
    towerjet4_tree->Branch("towerjet4_clohcal_eta", &_b_towerjet4_clohcal_eta);
    towerjet4_tree->Branch("towerjet4_clohcal_phi", &_b_towerjet4_clohcal_phi);
    towerjet4_tree->Branch("towerjet4_clohcal_dR", &_b_towerjet4_clohcal_dR);
    towerjet4_tree->Branch("towerjet4_clohcal_id", &_b_towerjet4_clohcal_id);
    towerjet4_tree->Branch("towerjet4_clohcal_ntwr", &_b_towerjet4_clohcal_ntwr);
    towerjet4_tree->Branch("towerjet4_clohcal_twrEsum", &_b_towerjet4_clohcal_twrEsum);
    towerjet4_tree->Branch("towerjet4_clohcal_twrNfrac", &_b_towerjet4_clohcal_twrNfrac);
  }
  if (_fill_bhhit_tree){
    bhhit_tree->Branch("event", &_b_event, "event/I");
    bhhit_tree->Branch("bhhit_n", &_b_bhhit_n, "bhhit_n/I");
    bhhit_tree->Branch("bhhit_edepTot", &_b_bhhit_edepTot,"bhhit_edepTot/F");
    bhhit_tree->Branch("bhhit_edep", &_b_bhhit_edep);
    bhhit_tree->Branch("bhhit_p", &_b_bhhit_p);
    bhhit_tree->Branch("bhhit_pt", &_b_bhhit_pt);
    bhhit_tree->Branch("bhhit_eta", &_b_bhhit_eta);
    bhhit_tree->Branch("bhhit_phi", &_b_bhhit_phi);
  }
  return 0;
}
////////////////////////////////////////////////////////////////////////////
int KyoTreeMaker::process_event(PHCompositeNode *topNode)
{
  _b_event++;
  std::cout << "*** KyoTreeMaker:: Processing Event: " << _b_event << std::endl;
    
  //// *** node: g4truthinfo
  PHG4TruthInfoContainer* truthinfos = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!truthinfos) {
    std::cout << PHWHERE << "G4TruthInfo node not found on node tree"<< std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  PHG4TruthInfoContainer::Range range = truthinfos->GetPrimaryParticleRange();
  //// *** true primary vertex
  PHG4VtxPoint *gvtx = truthinfos->GetPrimaryVtx(truthinfos->GetPrimaryVertexIndex());
  _b_gvtx_x =  gvtx->get_x(); 
  _b_gvtx_y =  gvtx->get_y(); 
  _b_gvtx_z =  gvtx->get_z(); 
  //// *** reco primary vertex
  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode,"GlobalVertexMap");
  if (!vertexmap) {
    std::cout << PHWHERE << "GlobalVertexMap node not found on node tree"<< std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  GlobalVertex* vtx = vertexmap->begin()->second;
  if (vtx) {
    _b_vtx_x =  vtx->get_x(); 
    _b_vtx_y =  vtx->get_y(); 
    _b_vtx_z =  vtx->get_z(); 
  } else {
    _b_vtx_x = 0;
    _b_vtx_y = 0;
    _b_vtx_z = 0;
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if (_fill_g4particle_tree){
    
    _b_g4particle_n = 0;
    _b_g4particle_e.clear();
    _b_g4particle_p.clear();
    _b_g4particle_pt.clear();
    _b_g4particle_eta.clear();
    _b_g4particle_phi.clear();
    _b_g4particle_pid.clear();

    for ( PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter ) {
      
      PHG4Particle* g4particle = iter->second;
    
      TLorentzVector trueptl_v; 
      trueptl_v.SetPxPyPzE( g4particle->get_px(), g4particle->get_py(), g4particle->get_pz(), g4particle->get_e() );
    
      if (fabs(trueptl_v.Eta()) > 1.1) continue;
   
      _b_g4particle_e.push_back( (float)trueptl_v.E() );
      _b_g4particle_p.push_back( (float)trueptl_v.P() );
      _b_g4particle_pt.push_back( (float)trueptl_v.Pt() );
      _b_g4particle_eta.push_back( (float) trueptl_v.Eta() );
      _b_g4particle_phi.push_back( (float) trueptl_v.Phi() );
      _b_g4particle_pid.push_back( (int)g4particle->get_pid() );
      
      _b_g4particle_n++;
    }
    g4particle_tree->Fill();
  }
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if (_fill_track_tree){
    
    _b_track_n = 0;
    _b_track_p.clear();
    _b_track_pt.clear();
    _b_track_eta.clear();
    _b_track_phi.clear();
    _b_track_charge.clear();
    _b_track_chisq.clear();
    _b_track_ndf.clear();
    _b_track_dca2d.clear();
    _b_cemc_E3x3.clear();
    _b_cemc_E3x3_PC.clear();
    _b_ihcal_E3x3.clear();
    _b_ohcal_E3x3.clear();
    _b_cemc_E5x5.clear();
    _b_cemc_E5x5_PC.clear();
    _b_ihcal_E5x5.clear();
    _b_ohcal_E5x5.clear();
    _b_cemc_clNtwr.clear();
    _b_cemc_clE.clear();
    _b_cemc_cldeta.clear();
    _b_cemc_cldphi.clear();
    _b_cemc_clE_PC.clear();
    _b_ihcal_clNtwr.clear();
    _b_ihcal_clE.clear();
    _b_ihcal_cldeta.clear();
    _b_ihcal_cldphi.clear();
    _b_ohcal_clNtwr.clear();
    _b_ohcal_clE.clear();
    _b_ohcal_cldeta.clear();
    _b_ohcal_cldphi.clear();
    
    SvtxTrackMap* trackmaps = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
    if (!trackmaps) {
      std::cout << PHWHERE << "SvtxTrackMap node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    for (SvtxTrackMap::ConstIter trk_iter =trackmaps->begin(); trk_iter != trackmaps->end(); trk_iter++){
      const SvtxTrack* track = trk_iter->second;
      
      _b_track_p.push_back( (float)track->get_p() );
      _b_track_pt.push_back( (float)track->get_pt() );
      _b_track_eta.push_back( (float)track->get_eta() );
      _b_track_phi.push_back( (float)track->get_phi() );
      _b_track_charge.push_back( (int)track->get_charge() );
      _b_track_chisq.push_back( (float)track->get_chisq() );
      _b_track_ndf.push_back( (int)track->get_ndf() );
      _b_track_dca2d.push_back( (float)track->get_dca2d() );
      //// 3x3 towers energy
      _b_cemc_E3x3.push_back( (float)track->get_cal_energy_3x3(SvtxTrack::CEMC) );    
      _b_cemc_E3x3_PC.push_back( PositionCorrectedNxNTowerE(topNode,track->get_cal_cluster_id(SvtxTrack::CEMC),track->get_cal_energy_3x3(SvtxTrack::CEMC)) );    
      _b_ihcal_E3x3.push_back( (float)track->get_cal_energy_3x3(SvtxTrack::HCALIN) );
      _b_ohcal_E3x3.push_back( (float)track->get_cal_energy_3x3(SvtxTrack::HCALOUT) );
      //// 5x5 towers energy
      _b_cemc_E5x5.push_back( (float)track->get_cal_energy_5x5(SvtxTrack::CEMC) );    
      _b_cemc_E5x5_PC.push_back( PositionCorrectedNxNTowerE(topNode,track->get_cal_cluster_id(SvtxTrack::CEMC),track->get_cal_energy_5x5(SvtxTrack::CEMC)) );    
      _b_ihcal_E5x5.push_back( (float)track->get_cal_energy_5x5(SvtxTrack::HCALIN) );
      _b_ohcal_E5x5.push_back( (float)track->get_cal_energy_5x5(SvtxTrack::HCALOUT) );
      //// cluster energy
      //_b_cemc_clE.push_back( track->get_cal_cluster_e(SvtxTrack::CEMC) ); // let's fill it at getClusterByIndex
      _b_cemc_cldeta.push_back( (float)track->get_cal_deta(SvtxTrack::CEMC) ); //a track extrapolated to calorimeter planes
      _b_cemc_cldphi.push_back( (float)track->get_cal_dphi(SvtxTrack::CEMC) );
      getClusterByIndex(topNode, "CEMC", track->get_cal_cluster_id(SvtxTrack::CEMC),track->get_eta(),track->get_phi());
      getClusterPCByIndex(topNode, track->get_cal_cluster_id(SvtxTrack::CEMC));
      //_b_ihcal_clE.push_back( track->get_cal_cluster_e(SvtxTrack::HCALIN) );
      _b_ihcal_cldeta.push_back( track->get_cal_deta(SvtxTrack::HCALIN) );
      _b_ihcal_cldphi.push_back( track->get_cal_dphi(SvtxTrack::HCALIN) );
      getClusterByIndex(topNode, "HCALIN", track->get_cal_cluster_id(SvtxTrack::HCALIN),track->get_eta(),track->get_phi());
      //_b_ohcal_clE.push_back( track->get_cal_cluster_e(SvtxTrack::HCALOUT) );
      _b_ohcal_cldeta.push_back( track->get_cal_deta(SvtxTrack::HCALOUT) );
      _b_ohcal_cldphi.push_back( track->get_cal_dphi(SvtxTrack::HCALOUT) );
      getClusterByIndex(topNode, "HCALOUT", track->get_cal_cluster_id(SvtxTrack::HCALOUT),track->get_eta(),track->get_phi());

      _b_track_n++;
    } 
    track_tree->Fill();
  }
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if (_fill_tower_tree){
    
    //// 1) CEMC
    _b_cemc_n = 0;
    _b_cemc_e.clear();
    _b_cemc_p.clear();
    _b_cemc_pt.clear();
    _b_cemc_eta.clear();
    _b_cemc_phi.clear();
    
    RawTowerContainer *towers_cemc = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_CEMC");
    if (!towers_cemc) {
      std::cout << PHWHERE << "TOWER_CALIB_CEMC node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerGeomContainer *towersgeom_cemc = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if (!towersgeom_cemc) {
      std::cout << PHWHERE << "TOWERGEOM_CEMC node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerContainer::ConstRange range_cemc = towers_cemc->getTowers(); 
      
    for ( RawTowerContainer::ConstIterator tower_iter = range_cemc.first; tower_iter != range_cemc.second; ++tower_iter ) {
      RawTower *twr_cemc = tower_iter->second; 
      RawTowerGeom *tgeom_cemc = towersgeom_cemc->get_tower_geometry(twr_cemc -> get_key()); 
      float r = tgeom_cemc->get_center_radius();
      float phi = atan2(tgeom_cemc->get_center_y(), tgeom_cemc->get_center_x()); 
      float z0 = tgeom_cemc->get_center_z();
      float z = z0 - _b_vtx_z;
      float eta = asinh(z/r); // eta after shift from vertex
      float energy = twr_cemc->get_energy();
      float pt = energy / cosh(eta);
      float px = pt * cos(phi);
      float py = pt * sin(phi);
      float pz = pt * sinh(eta);
      float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
      int id = (int)twr_cemc->get_id();
      
      _b_cemc_e.push_back( (float)energy );
      _b_cemc_p.push_back( (float)p );
      _b_cemc_pt.push_back( (float)pt );
      _b_cemc_eta.push_back( (float)eta );
      _b_cemc_phi.push_back( (float)phi );
      _b_cemc_id.push_back( (int)id );
      
      _b_cemc_n++;
   } 
    //// 2) HCALIN
    _b_ihcal_n = 0;
    _b_ihcal_e.clear();
    _b_ihcal_p.clear();
    _b_ihcal_pt.clear();
    _b_ihcal_eta.clear();
    _b_ihcal_phi.clear();
    
    RawTowerContainer *towers_ihcal = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_HCALIN");
    if (!towers_ihcal) {
      std::cout << PHWHERE << "TOWER_CALIB_HCALIN node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerGeomContainer *towersgeom_ihcal = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if (!towersgeom_ihcal) {
      std::cout << PHWHERE << "TOWERGEOM_HCALIN node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerContainer::ConstRange range_ihcal = towers_ihcal->getTowers(); 
      
    for ( RawTowerContainer::ConstIterator tower_iter = range_ihcal.first; tower_iter != range_ihcal.second; ++tower_iter ) {
      RawTower *twr_ihcal = tower_iter->second; 
      RawTowerGeom *tgeom_ihcal = towersgeom_ihcal->get_tower_geometry(twr_ihcal -> get_key()); 
      float r = tgeom_ihcal->get_center_radius();
      float phi = atan2(tgeom_ihcal->get_center_y(), tgeom_ihcal->get_center_x()); 
      float z0 = tgeom_ihcal->get_center_z();
      float z = z0 - _b_vtx_z;
      float eta = asinh(z/r); // eta after shift from vertex
      float energy = twr_ihcal->get_energy();
      float pt = energy / cosh(eta);
      float px = pt * cos(phi);
      float py = pt * sin(phi);
      float pz = pt * sinh(eta);
      float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
      
      _b_ihcal_e.push_back( (float)energy );
      _b_ihcal_p.push_back( (float)p );
      _b_ihcal_pt.push_back( (float)pt );
      _b_ihcal_eta.push_back( (float)eta );
      _b_ihcal_phi.push_back( (float)phi );
      
      _b_ihcal_n++;
   } 
    //// 3) HCALOUT
    _b_ohcal_n = 0;
    _b_ohcal_e.clear();
    _b_ohcal_p.clear();
    _b_ohcal_pt.clear();
    _b_ohcal_eta.clear();
    _b_ohcal_phi.clear();
    
    RawTowerContainer *towers_ohcal = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_HCALOUT");
    if (!towers_ohcal) {
      std::cout << PHWHERE << "TOWER_CALIB_HCALOUT node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerGeomContainer *towersgeom_ohcal = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if (!towersgeom_ohcal) {
      std::cout << PHWHERE << "TOWERGEOM_HCALOUT node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerContainer::ConstRange range_ohcal = towers_ohcal->getTowers(); 
      
    for ( RawTowerContainer::ConstIterator tower_iter = range_ohcal.first; tower_iter != range_ohcal.second; ++tower_iter ) {
      RawTower *twr_ohcal = tower_iter->second; 
      RawTowerGeom *tgeom_ohcal = towersgeom_ohcal->get_tower_geometry(twr_ohcal -> get_key()); 
      float r = tgeom_ohcal->get_center_radius();
      float phi = atan2(tgeom_ohcal->get_center_y(), tgeom_ohcal->get_center_x()); 
      float z0 = tgeom_ohcal->get_center_z();
      float z = z0 - _b_vtx_z;
      float eta = asinh(z/r); // eta after shift from vertex
      float energy = twr_ohcal->get_energy();
      float pt = energy / cosh(eta);
      float px = pt * cos(phi);
      float py = pt * sin(phi);
      float pz = pt * sinh(eta);
      float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
             
      _b_ohcal_e.push_back( (float)energy );
      _b_ohcal_p.push_back( (float)p );
      _b_ohcal_pt.push_back( (float)pt );
      _b_ohcal_eta.push_back( (float)eta );
      _b_ohcal_phi.push_back( (float)phi );
      
      _b_ohcal_n++;
   } 
    tower_tree->Fill();
  }
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if (_fill_cemc_cluster_tree){
    
    RawClusterContainer *clusters_cemc = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_CEMC");
    if (!clusters_cemc) {
      std::cout << PHWHERE << "CLUSTER_CEMC node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerContainer *towers_cemc = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_CEMC");
    if (!towers_cemc) {
      std::cout << PHWHERE << "TOWER_CALIB_CEMC node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerGeomContainer *clgm_cemc = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
    if (!clgm_cemc) {
      std::cout << PHWHERE << "TOWERGEOM_CEMC node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawClusterContainer::ConstRange clrange_cemc = clusters_cemc->getClusters(); 
    //RawTowerContainer::ConstRange range_cemc = towers_cemc->getTowers(); 
    
    for ( RawClusterContainer::ConstIterator cluster_iter = clrange_cemc.first; cluster_iter != clrange_cemc.second; ++cluster_iter ) {
      RawCluster *cl_cemc = cluster_iter->second; 
      float r = clgm_cemc->get_radius();
      float phi = cl_cemc->get_phi();
      float z0 = cl_cemc->get_z();
      float z = z0 - _b_vtx_z;
      float eta = asinh(z/r); // eta after shift from vertex

      float energy = cl_cemc->get_energy();
      float pt = energy / cosh(eta);
      float px = pt * cos(phi);
      float py = pt * sin(phi);
      float pz = pt * sinh(eta);
      float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
      int id = (int)cluster_iter->first;
      int ntwr = cl_cemc->getNTowers();
      float prob = cl_cemc->get_prob();
      float chi2 = cl_cemc->get_chi2();
      
      _b_cl_cemc_e = energy;
      _b_cl_cemc_p = p;
      _b_cl_cemc_pt = pt;
      _b_cl_cemc_eta = eta;
      _b_cl_cemc_phi = phi;
      _b_cl_cemc_id = id;
      _b_cl_cemc_ntwr = ntwr;
      _b_cl_cemc_prob = prob;
      _b_cl_cemc_chi2 = chi2;
     
      _b_cl_twr_cemc_e.clear(); 
      _b_cl_twr_cemc_id.clear(); 
      
      //loop over towers in the cluster
      float esum = 0;

      RawCluster::TowerConstRange cl_twrRange_cemc = cl_cemc->get_towers();
      //RawCluster::TowerConstIterator cl_twr_iter;
      for (RawCluster::TowerConstIterator cl_twr_iter = cl_twrRange_cemc.first; cl_twr_iter != cl_twrRange_cemc.second; ++ cl_twr_iter){
        RawTower *cl_twr_cemc = towers_cemc->getTower(cl_twr_iter->first);
        esum += cl_twr_cemc->get_energy(); 
        _b_cl_twr_cemc_e.push_back( (float)cl_twr_cemc->get_energy() );
        _b_cl_twr_cemc_id.push_back( (int)cl_twr_cemc->get_id() );
      
      }
      _b_cl_cemc_twrEsum = esum;
      
      /*
      //loop over tracks associated to the cluster
      SvtxTrackMap* trackmaps = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
      if (!trackmaps) {
        std::cout << PHWHERE << "SvtxTrackMap node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }    
      for (SvtxTrackMap::ConstIter trk_iter =trackmaps->begin(); trk_iter != trackmaps->end(); trk_iter++){
        const SvtxTrack* track = trk_iter->second;
        if ((int)track->get_cal_cluster_id(SvtxTrack::CEMC) == (int)cluster_iter->first) { 
          _b_cl_cemcMatched_e = energy;
          _b_cl_cemcMatched_p = p;
          _b_cl_cemcMatched_pt = pt;
          _b_cl_cemcMatched_eta = eta;
          _b_cl_cemcMatched_phi = phi;
          _b_cl_cemcMatched_id = id;
          _b_cl_cemcMatched_ntwr = ntwr;
          _b_cl_cemcMatched_prob = prob;
          _b_cl_cemcMatched_chi2 = chi2;
          _b_cl_cemcMatched_track_p = track->get_p();
          float deta = eta-track->get_eta();
          float dphi = phi-track->get_phi();
          _b_cl_cemcMatched_track_deta = deta;
          _b_cl_cemcMatched_track_dphi = dphi;
          float dr = sqrt( pow(deta,2)+pow(dphi,2) );
          _b_cl_cemcMatched_track_dr = dr;
        }
      }
      */
      cemc_cluster_tree->Fill();
    } 
  }
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if (_fill_ihcal_cluster_tree){
    RawClusterContainer *clusters_ihcal = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_HCALIN");
    if (!clusters_ihcal) {
      std::cout << PHWHERE << "CLUSTER_HCALIN node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerContainer *towers_ihcal = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_HCALIN");
    if (!towers_ihcal) {
      std::cout << PHWHERE << "TOWER_CALIB_HCALIN node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerGeomContainer *clgm_ihcal = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
    if (!clgm_ihcal) {
      std::cout << PHWHERE << "TOWERGEOM_HCALIN node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawClusterContainer::ConstRange clrange_ihcal = clusters_ihcal->getClusters(); 
      
    for ( RawClusterContainer::ConstIterator cluster_iter = clrange_ihcal.first; cluster_iter != clrange_ihcal.second; ++cluster_iter ) {
      RawCluster *cl_ihcal = cluster_iter->second; 
      float r = clgm_ihcal->get_radius();
      float phi = cl_ihcal->get_phi();
      float z0 = cl_ihcal->get_z();
      float z = z0 - _b_vtx_z;
      float eta = asinh(z/r); // eta after shift from vertex

      float energy = cl_ihcal->get_energy();
      float pt = energy / cosh(eta);
      float px = pt * cos(phi);
      float py = pt * sin(phi);
      float pz = pt * sinh(eta);
      float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
      int id = (int)cluster_iter->first;
      int ntwr = cl_ihcal->getNTowers();
      
      _b_cl_ihcal_e = energy;
      _b_cl_ihcal_p = p;
      _b_cl_ihcal_pt = pt;
      _b_cl_ihcal_eta = eta;
      _b_cl_ihcal_phi = phi;
      _b_cl_ihcal_id = id;
      _b_cl_ihcal_ntwr = ntwr;
       
      _b_cl_twr_ihcal_e.clear(); 
      _b_cl_twr_ihcal_id.clear(); 
      
      //loop over towers in the cluster
      RawCluster::TowerConstRange cl_twrRange_ihcal = cl_ihcal->get_towers();
      //RawCluster::TowerConstIterator cl_twr_iter;
      for (RawCluster::TowerConstIterator cl_twr_iter = cl_twrRange_ihcal.first; cl_twr_iter != cl_twrRange_ihcal.second; ++ cl_twr_iter){
        RawTower *cl_twr_ihcal = towers_ihcal->getTower(cl_twr_iter->first);
      
        _b_cl_twr_ihcal_e.push_back( (float)cl_twr_ihcal->get_energy() );
        _b_cl_twr_ihcal_id.push_back( (int)cl_twr_ihcal->get_id() );
      
      }
      ihcal_cluster_tree->Fill();

    }
  }
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if (_fill_ohcal_cluster_tree){
    RawClusterContainer *clusters_ohcal = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_HCALOUT");
    if (!clusters_ohcal) {
      std::cout << PHWHERE << "CLUSTER_HCALOUT node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerContainer *towers_ohcal = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_HCALOUT");
    if (!towers_ohcal) {
      std::cout << PHWHERE << "TOWER_CALIB_HCALOUT node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawTowerGeomContainer *clgm_ohcal = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
    if (!clgm_ohcal) {
      std::cout << PHWHERE << "TOWERGEOM_HCALOUT node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    RawClusterContainer::ConstRange clrange_ohcal = clusters_ohcal->getClusters(); 
      
    for ( RawClusterContainer::ConstIterator cluster_iter = clrange_ohcal.first; cluster_iter != clrange_ohcal.second; ++cluster_iter ) {
      RawCluster *cl_ohcal = cluster_iter->second; 
      float r = clgm_ohcal->get_radius();
      float phi = cl_ohcal->get_phi();
      float z0 = cl_ohcal->get_z();
      float z = z0 - _b_vtx_z;
      float eta = asinh(z/r); // eta after shift from vertex

      float energy = cl_ohcal->get_energy();
      float pt = energy / cosh(eta);
      float px = pt * cos(phi);
      float py = pt * sin(phi);
      float pz = pt * sinh(eta);
      float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
      int id = (int)cluster_iter->first;
      int ntwr = cl_ohcal->getNTowers();
       
      _b_cl_ohcal_e = energy;
      _b_cl_ohcal_p = p;
      _b_cl_ohcal_pt = pt;
      _b_cl_ohcal_eta = eta;
      _b_cl_ohcal_phi = phi;
      _b_cl_ohcal_id = id;
      _b_cl_ohcal_ntwr = ntwr;
       
      _b_cl_twr_ohcal_e.clear(); 
      _b_cl_twr_ohcal_id.clear(); 
      
      //loop over towers in the cluster
      RawCluster::TowerConstRange cl_twrRange_ohcal = cl_ohcal->get_towers();
      //RawCluster::TowerConstIterator cl_twr_iter;
      for (RawCluster::TowerConstIterator cl_twr_iter = cl_twrRange_ohcal.first; cl_twr_iter != cl_twrRange_ohcal.second; ++ cl_twr_iter){
        RawTower *cl_twr_ohcal = towers_ohcal->getTower(cl_twr_iter->first);
      
        _b_cl_twr_ohcal_e.push_back( (float)cl_twr_ohcal->get_energy() );
        _b_cl_twr_ohcal_id.push_back( (int)cl_twr_ohcal->get_id() );
      
      }
      ohcal_cluster_tree->Fill();
    } 
  }
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if (_fill_truthjet4_tree){
    
    JetMap* truth_jets = findNode::getClass<JetMap>(topNode,"AntiKt_Truth_r04");
    if (!truth_jets) {
      std::cout << PHWHERE << "AntiKt_Truth_r04 node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    
    for (JetMap::Iter jet_iter = truth_jets->begin(); jet_iter != truth_jets->end(); ++jet_iter) {
      Jet* this_jet = jet_iter->second;
      //// fill jet info
      float jet_e = this_jet->get_e();
      float jet_p = this_jet->get_p();
      float jet_pt = this_jet->get_pt();
      float jet_eta = this_jet->get_eta();
      float jet_phi = this_jet->get_phi();
      if (jet_pt < 5 || fabs(jet_eta) > 5) continue;
      _b_truthjet4_e = jet_e;
      _b_truthjet4_p = jet_p;
      _b_truthjet4_pt = jet_pt;
      _b_truthjet4_eta = jet_eta;
      _b_truthjet4_phi = jet_phi;
            
      ////  loop truth info for constituents
      //// *** node: g4truthinfo
      PHG4TruthInfoContainer* truthinfos_cons = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
      PHG4TruthInfoContainer::Range range_cons = truthinfos_cons->GetPrimaryParticleRange();
      if (!truthinfos_cons) {
        std::cout << PHWHERE << "G4TruthInfo node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      //// *** node: hepmcgenevt
      PHHepMCGenEventMap *genevtmap_cons = findNode::getClass<PHHepMCGenEventMap>(topNode,"PHHepMCGenEventMap"); 
      if (!genevtmap_cons) {
        std::cout << PHWHERE << "PHHepMCGenEventMap node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      
      combined_progenitors.clear(); //combined for one jet
      _b_truthjet4_cons_n = 0;
      _b_truthjet4_cons_e.clear();
      _b_truthjet4_cons_p.clear();
      _b_truthjet4_cons_pt.clear();
      _b_truthjet4_cons_eta.clear();
      _b_truthjet4_cons_phi.clear();
      _b_truthjet4_cons_dR.clear();
      _b_truthjet4_cons_pid.clear();
      
      for ( PHG4TruthInfoContainer::ConstIterator truth_iter = range_cons.first; truth_iter != range_cons.second; ++truth_iter ) {
        
        PHG4Particle* g4ptl_cons = truth_iter->second;
        
        //// loop jet components and get_track_id
        for (Jet::ConstIter comp_iter = this_jet->begin_comp(); comp_iter != this_jet->end_comp(); ++comp_iter) {
          
          if (g4ptl_cons->get_track_id()>=0 
          && (unsigned int)g4ptl_cons->get_track_id()==comp_iter->second){
    
            TLorentzVector cons_v; 
            cons_v.SetPxPyPzE( g4ptl_cons->get_px(), g4ptl_cons->get_py(), g4ptl_cons->get_pz(), g4ptl_cons->get_e() );
    
            float cons_e = cons_v.E();
            float cons_p = cons_v.P();
            float cons_pt = cons_v.Pt();
            float cons_eta = cons_v.Eta(); 
            float cons_phi = cons_v.Phi();
            
            float cons_dphi = cons_phi-jet_phi;
            if (cons_dphi > +M_PI) { cons_dphi -= 2 * M_PI; }
            if (cons_dphi < -M_PI) { cons_dphi += 2 * M_PI; }
            float cons_dR = sqrt(pow(cons_eta-jet_eta,2)+pow(cons_dphi,2)); 
            int cons_pid = g4ptl_cons->get_pid();

            _b_truthjet4_cons_e.push_back( (float)cons_e );
            _b_truthjet4_cons_p.push_back( (float)cons_p );
            _b_truthjet4_cons_pt.push_back( (float)cons_pt );
            _b_truthjet4_cons_eta.push_back( (float)cons_eta );
            _b_truthjet4_cons_phi.push_back( (float)cons_phi );
            _b_truthjet4_cons_dR.push_back( (float)cons_dR );
            _b_truthjet4_cons_pid.push_back( (int)cons_pid );
 
            //// Let's see if we can find the constituent's progenitor(s)
            progenitors.clear(); // for constituents
            GetJetPrimaryContributors( topNode, g4ptl_cons );

            //// Add this constituent's progenitors to the list of jet progenitors
            float prog_tot_e = 0.0;
            for(unsigned int n=0; n<progenitors.size(); n++){
              prog_tot_e += progenitors[n].pg_energy;
            }
            for(unsigned int n=0; n<progenitors.size(); n++){
              //// Is this progenitor already in the list?
              bool inList = false;
              unsigned int match = 0;
              for(unsigned int m=0; m<combined_progenitors.size(); m++){
                if(progenitors[n].barcode == combined_progenitors[m].barcode) {
                  inList = true;
                  match = m;
                  break;
                } 
              }
              if(inList){
                combined_progenitors[match].energy += cons_e*(progenitors[n].pg_energy/prog_tot_e);
              }
              else{
                JetProgenitor newProg = progenitors[n];
                newProg.energy = cons_e*(progenitors[n].pg_energy/prog_tot_e);
                combined_progenitors.push_back(newProg);
              }  
            }  
            _b_truthjet4_cons_n++;
          }//g4ptl == jet_component matched
        } // end loop of jet comp
      } // end loop of (all) g4ptl
            
      // Sort from highest energy to lowest energy 
      std::sort(combined_progenitors.begin(), combined_progenitors.end(), sortfunction);
      _b_truthjet4_pg_n = combined_progenitors.size();
      _b_truthjet4_pg_id.clear();
      _b_truthjet4_pg_fract.clear();
      _b_truthjet4_pg_status.clear();
            
      ////// Now pg info (for a jet) 
      for(unsigned int ip=0; ip<combined_progenitors.size(); ip++){ 
        _b_truthjet4_pg_id.push_back( (int)combined_progenitors[ip].pid ); 
        _b_truthjet4_pg_fract.push_back( (float)(combined_progenitors[ip].energy/jet_e) ); 
        _b_truthjet4_pg_status.push_back( (int)combined_progenitors[ip].status ); 
      }
      truthjet4_tree->Fill();
    }
  }    
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if (_fill_towerjet4_tree){
    
    JetMap* tower_jets = findNode::getClass<JetMap>(topNode,"AntiKt_Tower_r04");
    if (!tower_jets) {
      std::cout << PHWHERE << "AntiKt_Tower_r04 node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    for (JetMap::Iter jet_iter = tower_jets->begin(); jet_iter != tower_jets->end(); ++jet_iter) {
      
      Jet* this_jet = jet_iter->second;
     
      //// fill jet info
      float jet_e = this_jet->get_e();
      float jet_p = this_jet->get_p();
      float jet_pt = this_jet->get_pt();
      float jet_eta = this_jet->get_eta();
      float jet_phi = this_jet->get_phi();
      if (jet_pt < 5 || fabs(jet_eta) > 5) continue;
      _b_towerjet4_e = jet_e;
      _b_towerjet4_p = jet_p;
      _b_towerjet4_pt = jet_pt;
      _b_towerjet4_eta = jet_eta;
      _b_towerjet4_phi = jet_phi;
      //std::cout << ""<<std::endl;
      //std::cout << "************* EVENT ************** " << _b_event << std::endl;
      //std::cout << "jet_e = " << jet_e <<std::endl;
      //std::cout << "jet_eta = " << jet_eta <<std::endl;
      //std::cout << "jet_phi = " << jet_phi <<std::endl;
            
      ////////// 1) CEMC
      
      ////  a) loop tower info for constituents
      RawTowerContainer *towers_cemc = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_CEMC");
      if (!towers_cemc) {
        std::cout << PHWHERE << "TOWER_CALIB_CEMC node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      RawTowerGeomContainer *towersgeom_cemc = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
      if (!towersgeom_cemc) {
        std::cout << PHWHERE << "TOWERGEOM_CEMC node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      RawTowerContainer::ConstRange range_cemc = towers_cemc->getTowers(); 
      
      _b_towerjet4_cemc_n = 0;
      _b_towerjet4_cemc_e.clear(); 
      _b_towerjet4_cemc_p.clear(); 
      _b_towerjet4_cemc_pt.clear(); 
      _b_towerjet4_cemc_eta.clear(); 
      _b_towerjet4_cemc_phi.clear(); 
      _b_towerjet4_cemc_dR.clear(); 
      _b_towerjet4_cemc_id.clear(); 
      
      for ( RawTowerContainer::ConstIterator tower_iter = range_cemc.first; tower_iter != range_cemc.second; ++tower_iter ) {
        
        RawTower *twr_cemc = tower_iter->second; 
        RawTowerGeom *tgeom_cemc = towersgeom_cemc->get_tower_geometry(twr_cemc -> get_key()); 
        //// loop jet components and get_id
        for (Jet::ConstIter comp_iter = this_jet->begin_comp(); comp_iter != this_jet->end_comp(); ++comp_iter) {
          
          if (twr_cemc->get_id()>=0 
          && (unsigned int)twr_cemc->get_id()==comp_iter->second){
    
            float r = tgeom_cemc->get_center_radius();
            float phi = atan2(tgeom_cemc->get_center_y(), tgeom_cemc->get_center_x()); 
            float z0 = tgeom_cemc->get_center_z();
            float z = z0 - _b_vtx_z;
            float eta = asinh(z/r); // eta after shift from vertex
            float energy = twr_cemc->get_energy();
            float pt = energy / cosh(eta);
            float px = pt * cos(phi);
            float py = pt * sin(phi);
            float pz = pt * sinh(eta);
            float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
             
            float cemc_e = energy;
            float cemc_p = p;
            float cemc_pt = pt;
            float cemc_eta = eta;
            float cemc_phi = phi;
            
            float cemc_dphi = cemc_phi-jet_phi;
            if (cemc_dphi > +M_PI) { cemc_dphi -= 2 * M_PI; }
            if (cemc_dphi < -M_PI) { cemc_dphi += 2 * M_PI; }
            float cemc_dR = sqrt(pow(cemc_eta-jet_eta,2)+pow(cemc_dphi,2)); 
            int cemc_id = (int)twr_cemc->get_id(); 

            _b_towerjet4_cemc_e.push_back( (float)cemc_e );
            _b_towerjet4_cemc_p.push_back( (float)cemc_p );
            _b_towerjet4_cemc_pt.push_back( (float)cemc_pt );
            _b_towerjet4_cemc_eta.push_back( (float)cemc_eta );
            _b_towerjet4_cemc_phi.push_back( (float)cemc_phi );
            _b_towerjet4_cemc_dR.push_back( (float)cemc_dR );
            _b_towerjet4_cemc_id.push_back( (int)cemc_id );
            
            _b_towerjet4_cemc_n++;
          }
        } //EOL over jet cons
      } // EOL over tower
      
      //// b)  loop cluster info
      RawClusterContainer *clusters_cemc = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_CEMC"); 
      if (!clusters_cemc) {
        std::cout << PHWHERE << "CLUSTER_CEMC node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      RawClusterContainer::ConstRange clrange_cemc = clusters_cemc->getClusters();

      _b_towerjet4_clcemc_n = 0;
      _b_towerjet4_clcemc_e.clear(); 
      _b_towerjet4_clcemc_p.clear(); 
      _b_towerjet4_clcemc_pt.clear(); 
      _b_towerjet4_clcemc_eta.clear(); 
      _b_towerjet4_clcemc_phi.clear(); 
      _b_towerjet4_clcemc_dR.clear(); 
      _b_towerjet4_clcemc_id.clear(); 
      _b_towerjet4_clcemc_ntwr.clear(); 
      _b_towerjet4_clcemc_prob.clear(); 
      _b_towerjet4_clcemc_chi2.clear(); 
      _b_towerjet4_clcemc_twrEsum.clear(); 
      _b_towerjet4_clcemc_twrNfrac.clear(); 
      
      _b_towerjet4_clcemcMatched_n = 0;
      _b_towerjet4_clcemcMatched_e.clear();
      _b_towerjet4_clcemcMatched_p.clear();
      _b_towerjet4_clcemcMatched_pt.clear();
      _b_towerjet4_clcemcMatched_eta.clear();
      _b_towerjet4_clcemcMatched_phi.clear();
      _b_towerjet4_clcemcMatched_dR.clear(); 
      _b_towerjet4_clcemcMatched_id.clear();
      _b_towerjet4_clcemcMatched_ntwr.clear(); 
      _b_towerjet4_clcemcMatched_prob.clear(); 
      _b_towerjet4_clcemcMatched_chi2.clear(); 
      _b_towerjet4_clcemcMatched_track_p.clear();
      
      for ( RawClusterContainer::ConstIterator cluster_iter = clrange_cemc.first; cluster_iter != clrange_cemc.second; ++cluster_iter ) {
        
        RawCluster *cl_cemc = cluster_iter->second; 
        float r = towersgeom_cemc->get_radius();
        float phi = cl_cemc->get_phi();
        float z0 = cl_cemc->get_z();
        float z = z0 - _b_vtx_z;
        float eta = asinh(z/r); // eta after shift from vertex
        //float eta2 = RawClusterUtility::GetPseudorapidity(*cl_cemc, CLHEP::Hep3Vector(_b_vtx_x,_b_vtx_y,_b_vtx_z)); // eta after shift from vertex
        //std::cout << "cl_cemc eta = " << eta << std::endl; 
        //std::cout << "cl_cemc eta2 = " << eta2 << std::endl; 
        
        //// match with this_jet    
        float clcemc_dphi = phi-jet_phi;
        if (clcemc_dphi > +M_PI) { clcemc_dphi -= 2 * M_PI; }
        if (clcemc_dphi < -M_PI) { clcemc_dphi += 2 * M_PI; }
        float clcemc_dR = sqrt(pow(eta-jet_eta,2)+pow(clcemc_dphi,2)); 
        
        //std::cout << "   ********** ALL CLUSTER ************* " << std::endl;
        //std::cout << "   _b_towerjet4_clcemc_n = " << _b_towerjet4_clcemc_n <<std::endl;
        //std::cout << "   towerjet4_clcemc_id = " << (int)cluster_iter->first <<std::endl;
        //std::cout << "   towerjet4_clcemc_e = " << cl_cemc->get_energy() <<std::endl;
        //std::cout << "   towerjet4_clcemc_dR = " << clcemc_dR <<std::endl;

        if (clcemc_dR > 0.4) continue;
        
        float energy = cl_cemc->get_energy();
        float pt = energy / cosh(eta);
        float px = pt * cos(phi);
        float py = pt * sin(phi);
        float pz = pt * sinh(eta);
        float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
        int ntwr = cl_cemc->getNTowers();
        int id = (int)cluster_iter->first;
        float prob = cl_cemc->get_prob();
        float chi2 = cl_cemc->get_chi2();
      
        _b_towerjet4_clcemc_e.push_back( (float)energy );
        _b_towerjet4_clcemc_p.push_back( (float)p );
        _b_towerjet4_clcemc_pt.push_back( (float)pt );
        _b_towerjet4_clcemc_eta.push_back( (float)eta );
        _b_towerjet4_clcemc_phi.push_back( (float)phi );
        _b_towerjet4_clcemc_dR.push_back( (float)clcemc_dR );
        _b_towerjet4_clcemc_id.push_back( (int)id );
        _b_towerjet4_clcemc_ntwr.push_back( (float)ntwr );
        _b_towerjet4_clcemc_prob.push_back( (float)prob );
        _b_towerjet4_clcemc_chi2.push_back( (float)chi2 );
        //std::cout << "*** KYO Cl*** " << std::endl;
        //std::cout << "   towerjet4_clcemc_e = " << energy <<std::endl;
        //std::cout << "   towerjet4_clcemc_eta = " << eta <<std::endl;
        //std::cout << "   towerjet4_clcemc_phi = " << phi <<std::endl;
         
        ////loop over towers in the cluster
        float esum = 0;
        int numDen = 0;
        int numNum = 0;

        RawCluster::TowerConstRange cltwrRange_cemc = cl_cemc->get_towers();
        //RawCluster::TowerConstIterator cl_twr_iter;
        for (RawCluster::TowerConstIterator cltwr_iter = cltwrRange_cemc.first; cltwr_iter != cltwrRange_cemc.second; ++ cltwr_iter){
          
          RawTower *cltwr_cemc = towers_cemc->getTower(cltwr_iter->first);
          numDen++;  
          for (int itwr =0; itwr < (int)_b_towerjet4_cemc_id.size(); itwr ++) { // for cemc tower id
            if ( (int)_b_towerjet4_cemc_id.at(itwr) == (int)cltwr_iter->first ) {
              esum += cltwr_cemc->get_energy();
              numNum ++;  
            }
          }  
        }//EOL over towers in a cluster
        _b_towerjet4_clcemc_twrEsum.push_back( (float)esum );
        _b_towerjet4_clcemc_twrNfrac.push_back( (float)numNum/numDen );
        
        ////loop over tracks associated to the cluster
        SvtxTrackMap* trackmaps = findNode::getClass<SvtxTrackMap>(topNode,"SvtxTrackMap");
        if (!trackmaps) {
          std::cout << PHWHERE << "SvtxTrackMap node not found on node tree"<< std::endl;
          return Fun4AllReturnCodes::ABORTEVENT;
        }
        for (SvtxTrackMap::ConstIter trk_iter =trackmaps->begin(); trk_iter != trackmaps->end(); trk_iter++){
          const SvtxTrack* track = trk_iter->second;
          //std::cout << "   ********** ALL TRACK ************* " << std::endl;
          //std::cout << "   _b_towerjet4_clcemc_n = " << _b_towerjet4_clcemc_n <<std::endl;
          //if (_b_towerjet4_clcemc_n==0) std::cout << " *** track->get_cal_cluster_id(SvtxTrack::CEMC) : "<<track->get_cal_cluster_id(SvtxTrack::CEMC)<< std::endl;
          if (track->get_cal_cluster_id(SvtxTrack::CEMC) == cluster_iter->first) {
            //std::cout << " **************** MATCHED!!! track->get_cal_cluster_id(SvtxTrack::CEMC) : "<<track->get_cal_cluster_id(SvtxTrack::CEMC)<< std::endl;
            //std::cout << "track->get_p() = " << track->get_p() <<std::endl;
            //std::cout << "track->get_eta() = " << track->get_eta() <<std::endl;
            //std::cout << "track->get_phi() = " << track->get_phi() <<std::endl;
            //std::cout << "my deta = " << (track->get_eta()-eta) <<std::endl;
            //std::cout << "my dphi = " << (track->get_phi()-phi) <<std::endl;
            //std::cout << "my dr = " << sqrt(pow((track->get_eta()-eta),2)+pow((track->get_phi()-phi),2)) <<std::endl;
            //std::cout << "_clcemcMatched_id = " << id << std::endl;
            //std::cout << "_clcemcMatched_e = " << energy << std::endl;
            //std::cout << "_clcemcMatched_eta = " << eta << std::endl;
            //std::cout << "_clcemcMatched_phi = " << phi << std::endl;
            //std::cout << "_clcemcMatched_dR = " << clcemc_dR << std::endl;
            //std::cout << "track->get_cal_deta(SvtxTrack::CEMC) = " << track->get_cal_deta(SvtxTrack::CEMC) <<std::endl;
            //std::cout << "track->get_cal_dphi(SvtxTrack::CEMC) = " << track->get_cal_dphi(SvtxTrack::CEMC) <<std::endl;
            //std::cout << "track dr = " << sqrt(pow(track->get_cal_deta(SvtxTrack::CEMC),2)+pow(track->get_cal_dphi(SvtxTrack::CEMC),2)) <<std::endl;
            
            //// *** At this stage, store all clusters (double counted in case of two tracks being merged into the same cluster)
            //// *** It will be sorted at rootAna_sorting.C later
            //bool isInList = false;
            //for (int il =0; il< (int)_b_towerjet4_clcemcMatched_id.size(); il++){ //loop over clcemcMatched filled until now
            //  if ( ( (int)track->get_cal_cluster_id(SvtxTrack::CEMC) == _b_towerjet4_clcemcMatched_id.at(il) ) 
            //  && (track->get_p() < _b_towerjet4_clcemcMatched_track_p.at(il)) ) {
            //    isInList = true;
            //  }
            //}
            //if (!isInList){
            _b_towerjet4_clcemcMatched_n++;
            _b_towerjet4_clcemcMatched_e.push_back( (float)energy );
            _b_towerjet4_clcemcMatched_p.push_back( (float)p );
            _b_towerjet4_clcemcMatched_pt.push_back( (float)pt );
            _b_towerjet4_clcemcMatched_eta.push_back( (float)eta );
            _b_towerjet4_clcemcMatched_phi.push_back( (float)phi );
            _b_towerjet4_clcemcMatched_dR.push_back( (float)clcemc_dR );
            _b_towerjet4_clcemcMatched_id.push_back( (int)track->get_cal_cluster_id(SvtxTrack::CEMC) );
            _b_towerjet4_clcemcMatched_ntwr.push_back( (float)ntwr );
            _b_towerjet4_clcemcMatched_prob.push_back( (float)prob );
            _b_towerjet4_clcemcMatched_chi2.push_back( (float)chi2 );
            _b_towerjet4_clcemcMatched_track_p.push_back( (float)track->get_p() );
            _b_towerjet4_clcemcMatched_track_deta.push_back( (float)track->get_cal_deta(SvtxTrack::CEMC) );
            _b_towerjet4_clcemcMatched_track_dphi.push_back( (float)track->get_cal_dphi(SvtxTrack::CEMC) );
            float dr = sqrt(pow(track->get_cal_deta(SvtxTrack::CEMC),2)+pow(track->get_cal_dphi(SvtxTrack::CEMC),2)); 
            _b_towerjet4_clcemcMatched_track_dr.push_back( (float)dr );
            //float deta = eta-track->get_eta();
            //float dphi = phi-track->get_phi();
            //_b_towerjet4_clcemcMatched_track_deta.push_back( (float)deta );
            //_b_towerjet4_clcemcMatched_track_dphi.push_back( (float)dphi );
            //float dr = sqrt( pow(deta,2)+pow(dphi,2) ); 
            //_b_towerjet4_clcemcMatched_track_dr.push_back( (float)dr );
            //}
          }
        }
        _b_towerjet4_clcemc_n++;
      }//EOL over clusters within a jet 
      
      ////////// 2) HCALIN
      
      ////  a) loop tower info for constituents
      RawTowerContainer *towers_ihcal = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_HCALIN");
      if (!towers_ihcal) {
        std::cout << PHWHERE << "TOWER_CALIB_HCALIN node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      RawTowerGeomContainer *towersgeom_ihcal = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
      if (!towersgeom_ihcal) {
        std::cout << PHWHERE << "TOWERGEOM_HCALIN node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      RawTowerContainer::ConstRange range_ihcal = towers_ihcal->getTowers(); 
      
      _b_towerjet4_ihcal_n = 0;
      _b_towerjet4_ihcal_e.clear(); 
      _b_towerjet4_ihcal_p.clear(); 
      _b_towerjet4_ihcal_pt.clear(); 
      _b_towerjet4_ihcal_eta.clear(); 
      _b_towerjet4_ihcal_phi.clear(); 
      _b_towerjet4_ihcal_dR.clear(); 
      _b_towerjet4_ihcal_id.clear(); 
      
      for ( RawTowerContainer::ConstIterator tower_iter = range_ihcal.first; tower_iter != range_ihcal.second; ++tower_iter ) {
        
        RawTower *twr_ihcal = tower_iter->second; 
        RawTowerGeom *tgeom_ihcal = towersgeom_ihcal->get_tower_geometry(twr_ihcal -> get_key()); 
        //// loop jet components and get_id
        for (Jet::ConstIter comp_iter = this_jet->begin_comp(); comp_iter != this_jet->end_comp(); ++comp_iter) {
          
          if (twr_ihcal->get_id()>=0 
          && (unsigned int)twr_ihcal->get_id()==comp_iter->second){
    
            float r = tgeom_ihcal->get_center_radius();
            float phi = atan2(tgeom_ihcal->get_center_y(), tgeom_ihcal->get_center_x()); 
            float z0 = tgeom_ihcal->get_center_z();
            float z = z0 - _b_vtx_z;
            float eta = asinh(z/r); // eta after shift from vertex
            float energy = twr_ihcal->get_energy();
            float pt = energy / cosh(eta);
            float px = pt * cos(phi);
            float py = pt * sin(phi);
            float pz = pt * sinh(eta);
            float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
             
            float ihcal_e = energy;
            float ihcal_p = p;
            float ihcal_pt = pt;
            float ihcal_eta = eta;
            float ihcal_phi = phi;
            
            float ihcal_dphi = ihcal_phi-jet_phi;
            if (ihcal_dphi > +M_PI) { ihcal_dphi -= 2 * M_PI; }
            if (ihcal_dphi < -M_PI) { ihcal_dphi += 2 * M_PI; }
            float ihcal_dR = sqrt(pow(ihcal_eta-jet_eta,2)+pow(ihcal_dphi,2)); 
            int ihcal_id = (int)twr_ihcal->get_id(); 

            _b_towerjet4_ihcal_e.push_back( (float)ihcal_e );
            _b_towerjet4_ihcal_p.push_back( (float)ihcal_p );
            _b_towerjet4_ihcal_pt.push_back( (float)ihcal_pt );
            _b_towerjet4_ihcal_eta.push_back( (float)ihcal_eta );
            _b_towerjet4_ihcal_phi.push_back( (float)ihcal_phi );
            _b_towerjet4_ihcal_dR.push_back( (float)ihcal_dR );
            _b_towerjet4_ihcal_id.push_back( (float)ihcal_id );
            
            _b_towerjet4_ihcal_n++;
          }
        } // EOL over jet cons
      } // EOL over tower
      
      //// b)  loop cluster info
      RawClusterContainer *clusters_ihcal = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_HCALIN");
      if (!clusters_ihcal) {
        std::cout << PHWHERE << "CLUSTER_HCALIN node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      RawClusterContainer::ConstRange clrange_ihcal = clusters_ihcal->getClusters();

      _b_towerjet4_clihcal_n = 0;
      _b_towerjet4_clihcal_e.clear(); 
      _b_towerjet4_clihcal_p.clear(); 
      _b_towerjet4_clihcal_pt.clear(); 
      _b_towerjet4_clihcal_eta.clear(); 
      _b_towerjet4_clihcal_phi.clear(); 
      _b_towerjet4_clihcal_dR.clear(); 
      _b_towerjet4_clihcal_id.clear(); 
      _b_towerjet4_clihcal_ntwr.clear(); 
      _b_towerjet4_clihcal_twrEsum.clear(); 
      _b_towerjet4_clihcal_twrNfrac.clear(); 
      
      for ( RawClusterContainer::ConstIterator cluster_iter = clrange_ihcal.first; cluster_iter != clrange_ihcal.second; ++cluster_iter ) {
        
        RawCluster *cl_ihcal = cluster_iter->second; 
        float r = towersgeom_ihcal->get_radius();
        float phi = cl_ihcal->get_phi();
        float z0 = cl_ihcal->get_z();
        float z = z0 - _b_vtx_z;
        float eta = asinh(z/r); // eta after shift from vertex
        
        //// match with this_jet    
        float clihcal_dphi = phi-jet_phi;
        if (clihcal_dphi > +M_PI) { clihcal_dphi -= 2 * M_PI; }
        if (clihcal_dphi < -M_PI) { clihcal_dphi += 2 * M_PI; }
        float clihcal_dR = sqrt(pow(eta-jet_eta,2)+pow(clihcal_dphi,2)); 

        if (clihcal_dR > 0.4) continue;
        
        float energy = cl_ihcal->get_energy();
        float pt = energy / cosh(eta);
        float px = pt * cos(phi);
        float py = pt * sin(phi);
        float pz = pt * sinh(eta);
        float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
        int ntwr = cl_ihcal->getNTowers();
        int id = (int)cluster_iter->first;
         
        _b_towerjet4_clihcal_e.push_back( (float)energy );
        _b_towerjet4_clihcal_p.push_back( (float)p );
        _b_towerjet4_clihcal_pt.push_back( (float)pt );
        _b_towerjet4_clihcal_eta.push_back( (float)eta );
        _b_towerjet4_clihcal_phi.push_back( (float)phi );
        _b_towerjet4_clihcal_dR.push_back( (float)clihcal_dR );
        _b_towerjet4_clihcal_id.push_back( (float)id );
        _b_towerjet4_clihcal_ntwr.push_back( (float)ntwr );
     
        ////loop over towers in the cluster
        float esum = 0;
        int numDen = 0;
        int numNum = 0;

        RawCluster::TowerConstRange cltwrRange_ihcal = cl_ihcal->get_towers();
        //RawCluster::TowerConstIterator cl_twr_iter;
        for (RawCluster::TowerConstIterator cltwr_iter = cltwrRange_ihcal.first; cltwr_iter != cltwrRange_ihcal.second; ++ cltwr_iter){
          
          RawTower *cltwr_ihcal = towers_ihcal->getTower(cltwr_iter->first);
          numDen++;  
          for (int itwr =0; itwr < (int)_b_towerjet4_ihcal_id.size(); itwr ++) { // for ihcal tower id
            if ( (int)_b_towerjet4_ihcal_id.at(itwr) == (int)cltwr_iter->first ) {
              esum += cltwr_ihcal->get_energy();
              numNum ++;  
            }
          }  
        }//EOL over towers in a cluster
        
        _b_towerjet4_clihcal_twrEsum.push_back( (float)esum );
        _b_towerjet4_clihcal_twrNfrac.push_back( (float)numNum/numDen );
        
        _b_towerjet4_clihcal_n++;
      }//EOL over clusters within a jet 
      
      ////////// 3) HCALOUT
      
      ////  a) loop tower info for constituents
      RawTowerContainer *towers_ohcal = findNode::getClass<RawTowerContainer>(topNode,"TOWER_CALIB_HCALOUT");
      if (!towers_ohcal) {
        std::cout << PHWHERE << "TOWER_CALIB_HCALOUT node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      RawTowerGeomContainer *towersgeom_ohcal = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
      if (!towers_ohcal) {
        std::cout << PHWHERE << "TOWER_CALIB_HCALOUT node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      RawTowerContainer::ConstRange range_ohcal = towers_ohcal->getTowers(); 
      
      _b_towerjet4_ohcal_n = 0;
      _b_towerjet4_ohcal_e.clear(); 
      _b_towerjet4_ohcal_p.clear(); 
      _b_towerjet4_ohcal_pt.clear(); 
      _b_towerjet4_ohcal_eta.clear(); 
      _b_towerjet4_ohcal_phi.clear(); 
      _b_towerjet4_ohcal_dR.clear(); 
      _b_towerjet4_ohcal_id.clear(); 
      
      for ( RawTowerContainer::ConstIterator tower_iter = range_ohcal.first; tower_iter != range_ohcal.second; ++tower_iter ) {
        
        RawTower *twr_ohcal = tower_iter->second; 
        RawTowerGeom *tgeom_ohcal = towersgeom_ohcal->get_tower_geometry(twr_ohcal -> get_key()); 
        //// loop jet components and get_id
        for (Jet::ConstIter comp_iter = this_jet->begin_comp(); comp_iter != this_jet->end_comp(); ++comp_iter) {
          
          if (twr_ohcal->get_id()>=0 
          && (unsigned int)twr_ohcal->get_id()==comp_iter->second){
    
            float r = tgeom_ohcal->get_center_radius();
            float phi = atan2(tgeom_ohcal->get_center_y(), tgeom_ohcal->get_center_x()); 
            float z0 = tgeom_ohcal->get_center_z();
            float z = z0 - _b_vtx_z;
            float eta = asinh(z/r); // eta after shift from vertex
            float energy = twr_ohcal->get_energy();
            float pt = energy / cosh(eta);
            float px = pt * cos(phi);
            float py = pt * sin(phi);
            float pz = pt * sinh(eta);
            float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
             
            float ohcal_e = energy;
            float ohcal_p = p;
            float ohcal_pt = pt;
            float ohcal_eta = eta;
            float ohcal_phi = phi;
            
            float ohcal_dphi = ohcal_phi-jet_phi;
            if (ohcal_dphi > +M_PI) { ohcal_dphi -= 2 * M_PI; }
            if (ohcal_dphi < -M_PI) { ohcal_dphi += 2 * M_PI; }
            float ohcal_dR = sqrt(pow(ohcal_eta-jet_eta,2)+pow(ohcal_dphi,2)); 
            int ohcal_id = (int)twr_ohcal->get_id(); 

            _b_towerjet4_ohcal_e.push_back( (float)ohcal_e );
            _b_towerjet4_ohcal_p.push_back( (float)ohcal_p );
            _b_towerjet4_ohcal_pt.push_back( (float)ohcal_pt );
            _b_towerjet4_ohcal_eta.push_back( (float)ohcal_eta );
            _b_towerjet4_ohcal_phi.push_back( (float)ohcal_phi );
            _b_towerjet4_ohcal_dR.push_back( (float)ohcal_dR );
            _b_towerjet4_ohcal_id.push_back( (float)ohcal_id );
            
            _b_towerjet4_ohcal_n++;
          }
        } // EOL over jet cons
      } // EOL over tower
      
      //// b)  loop cluster info
      RawClusterContainer *clusters_ohcal = findNode::getClass<RawClusterContainer>(topNode,"CLUSTER_HCALOUT");
        if (!clusters_ohcal) {
        std::cout << PHWHERE << "CLUSTER_HCALOUT node not found on node tree"<< std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      RawClusterContainer::ConstRange clrange_ohcal = clusters_ohcal->getClusters();

      _b_towerjet4_clohcal_n = 0;
      _b_towerjet4_clohcal_e.clear(); 
      _b_towerjet4_clohcal_p.clear(); 
      _b_towerjet4_clohcal_pt.clear(); 
      _b_towerjet4_clohcal_eta.clear(); 
      _b_towerjet4_clohcal_phi.clear(); 
      _b_towerjet4_clohcal_dR.clear(); 
      _b_towerjet4_clohcal_id.clear(); 
      _b_towerjet4_clohcal_ntwr.clear(); 
      _b_towerjet4_clohcal_twrEsum.clear(); 
      _b_towerjet4_clohcal_twrNfrac.clear(); 
      
      for ( RawClusterContainer::ConstIterator cluster_iter = clrange_ohcal.first; cluster_iter != clrange_ohcal.second; ++cluster_iter ) {
        
        RawCluster *cl_ohcal = cluster_iter->second; 
        float r = towersgeom_ohcal->get_radius();
        float phi = cl_ohcal->get_phi();
        float z0 = cl_ohcal->get_z();
        float z = z0 - _b_vtx_z;
        float eta = asinh(z/r); // eta after shift from vertex
        
        //// match with this_jet    
        float clohcal_dphi = phi-jet_phi;
        if (clohcal_dphi > +M_PI) { clohcal_dphi -= 2 * M_PI; }
        if (clohcal_dphi < -M_PI) { clohcal_dphi += 2 * M_PI; }
        float clohcal_dR = sqrt(pow(eta-jet_eta,2)+pow(clohcal_dphi,2)); 

        if (clohcal_dR > 0.4) continue;
        
        float energy = cl_ohcal->get_energy();
        float pt = energy / cosh(eta);
        float px = pt * cos(phi);
        float py = pt * sin(phi);
        float pz = pt * sinh(eta);
        float p = sqrt(pow(px,2)+pow(py,2)+pow(pz,2));
        int ntwr = cl_ohcal->getNTowers();
        int id = (int)cluster_iter->first;
      
        _b_towerjet4_clohcal_e.push_back( (float)energy );
        _b_towerjet4_clohcal_p.push_back( (float)p );
        _b_towerjet4_clohcal_pt.push_back( (float)pt );
        _b_towerjet4_clohcal_eta.push_back( (float)eta );
        _b_towerjet4_clohcal_phi.push_back( (float)phi );
        _b_towerjet4_clohcal_dR.push_back( (float)clohcal_dR );
        _b_towerjet4_clohcal_id.push_back( (float)id );
        _b_towerjet4_clohcal_ntwr.push_back( (float)ntwr );
     
        ////loop over towers in the cluster
        float esum = 0;
        int numDen = 0;
        int numNum = 0;

        RawCluster::TowerConstRange cltwrRange_ohcal = cl_ohcal->get_towers();
        //RawCluster::TowerConstIterator cl_twr_iter;
        for (RawCluster::TowerConstIterator cltwr_iter = cltwrRange_ohcal.first; cltwr_iter != cltwrRange_ohcal.second; ++ cltwr_iter){
          
          RawTower *cltwr_ohcal = towers_ohcal->getTower(cltwr_iter->first);
          numDen++;  
          for (int itwr =0; itwr < (int)_b_towerjet4_ohcal_id.size(); itwr ++) { // for ohcal tower id
            if ( (int)_b_towerjet4_ohcal_id.at(itwr) == (int)cltwr_iter->first ) {
              esum += cltwr_ohcal->get_energy();
              numNum ++;  
            }
          }  
        }//EOL over towers in a cluster
        
        _b_towerjet4_clohcal_twrEsum.push_back( (float)esum );
        _b_towerjet4_clohcal_twrNfrac.push_back( (float)numNum/numDen );
        
        _b_towerjet4_clohcal_n++;
      }//EOL over clusters within a jet 
      
      towerjet4_tree->Fill();
    }
  }    
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if (_fill_bhhit_tree){
    _b_bhhit_n = 0;
    _b_bhhit_edepTot = 0;
    _b_bhhit_edep.clear();
    _b_bhhit_p.clear();
    _b_bhhit_pt.clear();
    _b_bhhit_eta.clear();
    _b_bhhit_phi.clear();

    PHG4HitContainer* bh_hits = findNode::getClass<PHG4HitContainer> (topNode, "G4HIT_BH_1");
    if (!bh_hits) {
      std::cerr << PHWHERE << "G4HIT_BH_1 node not found on node tree"<< std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    PHG4HitContainer::ConstRange range_bh = bh_hits->getHits();   
    for (PHG4HitContainer::ConstIterator bh_iter = range_bh.first; bh_iter != range_bh.second; bh_iter++){
      PHG4Hit *this_hit =  bh_iter->second;
      float x = this_hit->get_x(0);
      float y = this_hit->get_y(0);
      float z = this_hit->get_z(0);
      float r = sqrt(pow(x,2)+pow(y,2));
      float px = this_hit->get_px(0);
      float py = this_hit->get_py(0);
      float pz = this_hit->get_pz(0);
      
      _b_bhhit_edepTot += this_hit->get_edep();
      _b_bhhit_edep.push_back( (float)this_hit->get_edep() );
      _b_bhhit_p.push_back( (float)sqrt(pow(px,2)+pow(py,2)+pow(pz,2)) );
      _b_bhhit_pt.push_back( (float)sqrt(pow(px,2)+pow(py,2)) );
      _b_bhhit_eta.push_back( (float)asinh(z/r) );
      _b_bhhit_phi.push_back( (float)atan2(y,x) );

      _b_bhhit_n ++;
    }
    bhhit_tree->Fill();
  }    
  
  return 0;
}
////////////////////////////////////////////////////////////////////////////
int KyoTreeMaker::End(PHCompositeNode *topNode)
{

  _f->Write();
  _f->Close();

  return 0;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

float KyoTreeMaker::PositionCorrectedNxNTowerE( PHCompositeNode *topNode, unsigned int cl_id, float eNN){

  // pull the clusters
  // Assumes a 1-1 correlation between the two lists!
  std::string clusternodename = "CLUSTER_CEMC";
  RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
  if (!clusters) {
    std::cout << PHWHERE << " ERROR: Can't find node " << clusternodename << std::endl;
    return 0.0;
  }    
  clusternodename = "CLUSTER_POS_COR_CEMC";
  RawClusterContainer *clusters_PC = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
  if (!clusters_PC) {
    std::cout << PHWHERE << " ERROR: Can't find node " << clusternodename << std::endl;
    return 0.0;
  }    

  if((cl_id>=0) && (cl_id<clusters->size())){
    RawCluster *cluster = clusters->getCluster(cl_id);
    RawCluster *cluster_PC = clusters_PC->getCluster(cl_id);
    return (cluster_PC->get_energy()/cluster->get_energy())*eNN; 
  }
  else
    return 0.0; 
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void KyoTreeMaker::getClusterByIndex( PHCompositeNode *topNode, std::string detName, unsigned int cl_id, float teta, float tphi){

  // pull the clusters
  std::string clusternodename = "CLUSTER_" + detName;
  RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
  if (!clusters) {
    std::cout << PHWHERE << " ERROR: Can't find node " << clusternodename << std::endl;
    return;
  }    
  std::string geonodename = "TOWERGEOM_" + detName;
  RawTowerGeomContainer *towersgeo = findNode::getClass<RawTowerGeomContainer>(topNode,geonodename.c_str());
  if (!towersgeo) {
    std::cout << PHWHERE << " ERROR: Can't find node " << geonodename << std::endl;
    return;
  }    
  
  float cl_e = -531.0;
  int cl_n = -531; 
  //float r = -531;
  //float phi = -531;
  //float z0 = -531;
  //float z = -531;
  //float eta = -531;
  //float deta = -531;
  //float dphi = -531;
  if((cl_id>=0) && (cl_id<clusters->size())){
    RawCluster *cluster = clusters->getCluster(cl_id);
    cl_e = cluster->get_energy();
    cl_n = cluster->getNTowers(); 
    //r = towersgeo->get_radius();
    //phi = cluster->get_phi();
    //z0 = cluster->get_z();
    //z = z0 - _b_vtx_z;
    //eta = asinh(z/r);
    //deta = eta-teta;
    //dphi = phi-tphi;
  }
  
  if(detName == "CEMC"){
    _b_cemc_clE.push_back( (float)cl_e );
    _b_cemc_clNtwr.push_back( (int)cl_n ); 
    //_b_cemc_cldeta.push_back( (float)deta );
    //_b_cemc_cldphi.push_back( (float)dphi );
  }
  else if(detName == "HCALIN"){
    _b_ihcal_clE.push_back( (float)cl_e );
    _b_ihcal_clNtwr.push_back( (int)cl_n ); 
    //_b_ihcal_cldeta.push_back( (float)deta );
    //_b_ihcal_cldphi.push_back( (float)dphi );
  }
  else if(detName == "HCALOUT"){
    _b_ohcal_clE.push_back( (float)cl_e );
    _b_ohcal_clNtwr.push_back( (int)cl_n ); 
    //_b_ohcal_cldeta.push_back( (float)deta );
    //_b_ohcal_cldphi.push_back( (float)dphi );
  }
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void KyoTreeMaker:: getClusterPCByIndex( PHCompositeNode *topNode, unsigned int cl_id){

  // only for CEMC
  // pull the clusters
  std::string clusternodename = "CLUSTER_POS_COR_CEMC";
  RawClusterContainer *clusters_PC = findNode::getClass<RawClusterContainer>(topNode,clusternodename.c_str());
  if (!clusters_PC) {
    std::cout << PHWHERE << " ERROR: Can't find node " << clusternodename << std::endl;
    return;
  }    

  float cl_e = 0.0;
  if((cl_id>=0) && (cl_id<clusters_PC->size())){
    RawCluster *cluster_PC = clusters_PC->getCluster(cl_id); 
    cl_e = cluster_PC->get_energy();
  }
  _b_cemc_clE_PC.push_back( (float)cl_e );
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void KyoTreeMaker::GetJetPrimaryContributors( PHCompositeNode *topNode, PHG4Particle* g4ptl){

  PHG4Particle* g4particle = g4ptl; 
  PHNodeIterator iter(topNode);
  
  //// *** node: g4truthinfo
  PHG4TruthInfoContainer* truthinfos = findNode::getClass<PHG4TruthInfoContainer>(topNode,"G4TruthInfo");
  if (!truthinfos) {
    std::cout << PHWHERE << "G4TruthInfo node not found on node tree"<< std::endl;
    //return Fun4AllReturnCodes::ABORTEVENT;
    return;
  }
  //// *** node: hepmcgenevt
  PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode,"PHHepMCGenEventMap"); 
  if (!genevtmap) {
    std::cout << PHWHERE << "PHHepMCGenEventMap node not found on node tree"<< std::endl;
    //return Fun4AllReturnCodes::ABORTEVENT;
    return;
  }

  PHHepMCGenEventMap::ConstReverseIter hepmciter = genevtmap->rbegin();
  //int testing = 0;
  if (hepmciter != genevtmap->rend())
  { 
    PHHepMCGenEvent * hepmc_evt = hepmciter->second;
    HepMC::GenEvent *evt = hepmc_evt->getEvent();
    HepMC::GenParticle *ptcle = NULL;
    // Track the g4particle back to its primary
    //std::cout << "Starting Particle ->" << g4particle->get_pid() << std::endl; 
    while( !truthinfos->is_primary(g4particle) ){
      // All particles
      PHG4Particle *prev_g4particle = g4particle; 
      PHG4TruthInfoContainer::ConstRange range = truthinfos->GetParticleRange();//all
      for (PHG4TruthInfoContainer::ConstIterator truth_itr = range.first;
	      truth_itr != range.second; ++truth_itr) {
        PHG4Particle* mother = truth_itr->second;
        if(g4particle->get_parent_id()==mother->get_track_id()){
	        g4particle = mother;  
	        break; 
        }
      }
      if(prev_g4particle==g4particle) return; // pure secondary?
      //std::cout << "->" << g4particle->get_pid() << std::endl; 
    }
  
    //// associate Geant4 primary particle with its HEPMC primary
    HepMC::GenEvent::particle_iterator pitr = evt->particles_begin();
    for( ; pitr!=evt->particles_end(); pitr++){
      ptcle = *pitr;        
      if (ptcle->barcode()  == g4particle->get_barcode() ) { break; } //match
    }
    if(!ptcle) {
      std::cout << " No HepMC primary found associated with G4 primary! " << std::endl;
      return;
    }
    //std::cout << "PARTICLE:" << std::endl;
    //ptcle->print();
    //std::cout << ""<<std::endl;
    //std::cout << "map iter...? (should be zero..)" <<testing<<std::endl;
    //testing++;
    
    // Recursively trace out the particle tree
    TraceHepMCParticle(ptcle);
  } 
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void KyoTreeMaker::TraceHepMCParticle( HepMC::GenParticle *ptcle){
  
  ////////// Important Pythia 8 Status Codes: 
  // 23: outgoing of hardest subprocess <-
  // 33: outgoing of subsequent subprocesses <-
  // 43: outgoing from an initial state branching <-
  // 51: outgoing from a final state branching 
  // 63: outgoing beam remnant <-
  ////////// Pythia6 Status Codes:
  // 1: Stable final state particle
  // 2/10902: Unstable particle
  // 3: Documentary particle <- 
  //////////

  if(ptcle->production_vertex()){
    for ( HepMC::GenVertex::particle_iterator mother 
	  = ptcle->production_vertex()->particles_begin(HepMC::parents);
	  mother != ptcle->production_vertex()->particles_end(HepMC::parents); 
	  ++mother ) {
      //std::cout << "->" << "particle = " << ptcle->pdg_id() << " status = " << ptcle->status() << " parent = " << (*mother)->pdg_id() << " parent status = " << (*mother)->status() << std::endl; 
      
      bool pythia8=true;
      bool pythia6=false;
      if(pythia8){
	      // These status codes terminate tracing 
	      if( ((*mother)->status()==23) || 
	          ((*mother)->status()==33) ||
	          ((*mother)->status()==43) ||
	          ((*mother)->status()==63) ){ 
          JetProgenitor jprog;
          jprog.pid = (*mother)->pdg_id();
          jprog.status = (*mother)->status();
          jprog.barcode = (*mother)->barcode();
          jprog.energy = 0.0;
          jprog.pg_energy = (*mother)->momentum().e();
          progenitors.push_back(jprog);
	      } else { TraceHepMCParticle(*mother);} 
      }//end of pythia8
      else if (pythia6){
	      if( ((*mother)->status()==3) ){
          JetProgenitor jprog;
          jprog.pid = (*mother)->pdg_id();
          jprog.status = (*mother)->status();
          jprog.barcode = (*mother)->barcode();
          jprog.energy = 0.0;
          jprog.pg_energy = (*mother)->momentum().e();
          progenitors.push_back(jprog);
        } else { TraceHepMCParticle(*mother);} 
      }//end of pythia6
      else {
      }
    }// end of mother iter
  }// end of production_vertex()
}

