// Rough copy paste from HGCalAnalysis
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/Math/interface/deltaPhi.h"
// track data formats
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/GeometrySurface/interface/PlaneBuilder.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CommonTools/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "FastSimulation/CaloGeometryTools/interface/Transform3DPJ.h"
#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/Event/interface/FSimVertex.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/VolumeGeometry/interface/MagVolumeOutsideValidity.h"
#include "TrackPropagation/RungeKutta/interface/defaultRKPropagator.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1F.h"
#include "TTree.h"

#include <map>
#include <set>
#include <string>
#include <vector>

using namespace std;

class PFAnalysis : public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::SharedResources> {
 public:
  typedef ROOT::Math::Transform3DPJ::Point Point;

  PFAnalysis();
  explicit PFAnalysis(const edm::ParameterSet &);
  ~PFAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  virtual void beginRun(edm::Run const &iEvent, edm::EventSetup const &) override;
  virtual void endRun(edm::Run const &iEvent, edm::EventSetup const &) override;

 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event &, const edm::EventSetup &) override;
  virtual void endJob() override;

  void clearVariables();

  // ---------parameters ----------------------------
  bool readGen_;
  bool storeMoreGenInfo_;
  bool storeGenParticleExtrapolation_;
  bool storePFCandidates_;
  bool storeGunParticles_;
  bool verbose_;

  // ----------member data ---------------------------

  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticles_;
  edm::EDGetTokenT<std::vector<SimTrack>> simTracks_;
  edm::EDGetTokenT<std::vector<SimVertex>> simVertices_;
  edm::EDGetTokenT<edm::HepMCProduct> hev_;
  edm::EDGetTokenT<std::vector<reco::Track>> tracks_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertices_;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandidates_;

  TTree *t_;

  ////////////////////
  // event
  //
  edm::RunNumber_t ev_run_;
  edm::LuminosityBlockNumber_t ev_lumi_;
  edm::EventNumber_t ev_event_;
  float vtx_x_;
  float vtx_y_;
  float vtx_z_;

  ////////////////////
  // GenParticles
  //
  std::vector<float> genpart_eta_;
  std::vector<float> genpart_phi_;
  std::vector<float> genpart_pt_;
  std::vector<float> genpart_energy_;
  std::vector<float> genpart_dvx_;
  std::vector<float> genpart_dvy_;
  std::vector<float> genpart_dvz_;
  std::vector<float> genpart_ovx_;
  std::vector<float> genpart_ovy_;
  std::vector<float> genpart_ovz_;
  std::vector<float> genpart_exx_;
  std::vector<float> genpart_exy_;
  std::vector<int> genpart_mother_;
  std::vector<float> genpart_exphi_;
  std::vector<float> genpart_exeta_;
  std::vector<float> genpart_fbrem_;
  std::vector<int> genpart_pid_;
  std::vector<int> genpart_gen_;
  std::vector<int> genpart_reachedEE_;
  std::vector<bool> genpart_fromBeamPipe_;
  std::vector<std::vector<float>> genpart_posx_;
  std::vector<std::vector<float>> genpart_posy_;
  std::vector<std::vector<float>> genpart_posz_;

  ////////////////////
  // reco::GenParticles
  //
  std::vector<float> gen_eta_;
  std::vector<float> gen_phi_;
  std::vector<float> gen_pt_;
  std::vector<float> gen_energy_;
  std::vector<int> gen_charge_;
  std::vector<int> gen_pdgid_;
  std::vector<int> gen_status_;
  std::vector<std::vector<int>> gen_daughters_;

  ////////////////////
  // high purity tracks
  //
  std::vector<float> track_eta_;
  std::vector<float> track_phi_;
  std::vector<float> track_pt_;
  std::vector<float> track_energy_;
  std::vector<int> track_charge_;
  std::vector<std::vector<float>> track_posx_;
  std::vector<std::vector<float>> track_posy_;
  std::vector<std::vector<float>> track_posz_;

  ////////////////////
  // PF candidates
  //
  std::vector<float> pfcandidate_eta_;
  std::vector<float> pfcandidate_phi_;
  std::vector<float> pfcandidate_pt_;
  std::vector<float> pfcandidate_energy_;
  std::vector<int> pfcandidate_pdgid_;

  ////////////////////
  // gun particles per vertex
  //
  std::vector<std::vector<int> > gunparticle_id_;
  std::vector<std::vector<float> > gunparticle_energy_;
  std::vector<std::vector<float> > gunparticle_pt_;
  std::vector<std::vector<float> > gunparticle_eta_;
  std::vector<std::vector<float> > gunparticle_phi_;


  // -------convenient tool to deal with simulated tracks
  FSimEvent *mySimEvent_;
  edm::ParameterSet particleFilter_;
  std::vector<float> layerPositions_;


  // and also the magnetic field
  MagneticField const *aField_;

};

PFAnalysis::PFAnalysis() { ; }

PFAnalysis::PFAnalysis(const edm::ParameterSet &iConfig)
    : readGen_(iConfig.getParameter<bool>("readGenParticles")),
      storeMoreGenInfo_(iConfig.getParameter<bool>("storeGenParticleOrigin")),
      storeGenParticleExtrapolation_(iConfig.getParameter<bool>("storeGenParticleExtrapolation")),
      storePFCandidates_(iConfig.getParameter<bool>("storePFCandidates")),
      storeGunParticles_(iConfig.getParameter<bool>("storeGunParticles")),
      verbose_(iConfig.getParameter<bool>("verbose")),
      particleFilter_(iConfig.getParameter<edm::ParameterSet>("TestParticleFilter")) {
  // now do what ever initialization is needed
  mySimEvent_ = new FSimEvent(particleFilter_);

  hev_ = consumes<edm::HepMCProduct>(edm::InputTag("generatorSmeared"));

  simTracks_ = consumes<std::vector<SimTrack>>(edm::InputTag("g4SimHits"));
  simVertices_ = consumes<std::vector<SimVertex>>(edm::InputTag("g4SimHits"));

  if (readGen_) {
    genParticles_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
  }

  if (storePFCandidates_) {
    pfCandidates_ = consumes<std::vector<reco::PFCandidate>>(edm::InputTag("particleFlow"));
  }

  tracks_ = consumes<std::vector<reco::Track>>(edm::InputTag("generalTracks"));
  vertices_ = consumes<std::vector<reco::Vertex>>(edm::InputTag("offlinePrimaryVertices"));

  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  fs->make<TH1F>("total", "total", 100, 0, 5.);

  t_ = fs->make<TTree>("pftree", "pftree");

  // event info
  t_->Branch("event", &ev_event_);
  t_->Branch("lumi", &ev_lumi_);
  t_->Branch("run", &ev_run_);
  t_->Branch("vtx_x", &vtx_x_);
  t_->Branch("vtx_y", &vtx_y_);
  t_->Branch("vtx_z", &vtx_z_);

  t_->Branch("genpart_eta", &genpart_eta_);
  t_->Branch("genpart_phi", &genpart_phi_);
  t_->Branch("genpart_pt", &genpart_pt_);
  t_->Branch("genpart_energy", &genpart_energy_);
  t_->Branch("genpart_dvx", &genpart_dvx_);
  t_->Branch("genpart_dvy", &genpart_dvy_);
  t_->Branch("genpart_dvz", &genpart_dvz_);
  if (storeMoreGenInfo_) {
    t_->Branch("genpart_ovx", &genpart_ovx_);
    t_->Branch("genpart_ovy", &genpart_ovy_);
    t_->Branch("genpart_ovz", &genpart_ovz_);
    t_->Branch("genpart_mother", &genpart_mother_);
  }
  if (storeGenParticleExtrapolation_) {
    t_->Branch("genpart_exphi", &genpart_exphi_);
    t_->Branch("genpart_exeta", &genpart_exeta_);
    t_->Branch("genpart_exx", &genpart_exx_);
    t_->Branch("genpart_exy", &genpart_exy_);
  }
  t_->Branch("genpart_fbrem", &genpart_fbrem_);
  t_->Branch("genpart_pid", &genpart_pid_);
  t_->Branch("genpart_gen", &genpart_gen_);
  t_->Branch("genpart_reachedEE", &genpart_reachedEE_);
  t_->Branch("genpart_fromBeamPipe", &genpart_fromBeamPipe_);
  t_->Branch("genpart_posx", &genpart_posx_);
  t_->Branch("genpart_posy", &genpart_posy_);
  t_->Branch("genpart_posz", &genpart_posz_);


  if (readGen_) {
    t_->Branch("gen_eta", &gen_eta_);
    t_->Branch("gen_phi", &gen_phi_);
    t_->Branch("gen_pt", &gen_pt_);
    t_->Branch("gen_energy", &gen_energy_);
    t_->Branch("gen_charge", &gen_charge_);
    t_->Branch("gen_pdgid", &gen_pdgid_);
    t_->Branch("gen_status", &gen_status_);
    t_->Branch("gen_daughters", &gen_daughters_);
  }

  ////////////////////
  // high purity trackstatep
  //
  t_->Branch("track_eta", &track_eta_);
  t_->Branch("track_phi", &track_phi_);
  t_->Branch("track_pt", &track_pt_);
  t_->Branch("track_energy", &track_energy_);
  t_->Branch("track_charge", &track_charge_);
  t_->Branch("track_posx", &track_posx_);
  t_->Branch("track_posy", &track_posy_);
  t_->Branch("track_posz", &track_posz_);

  ////////////////////
  // PF candidates
  //
  if (storePFCandidates_) {
    t_->Branch("pfcandidate_eta", &pfcandidate_eta_);
    t_->Branch("pfcandidate_phi", &pfcandidate_phi_);
    t_->Branch("pfcandidate_pt", &pfcandidate_pt_);
    t_->Branch("pfcandidate_energy", &pfcandidate_energy_);
    t_->Branch("pfcandidate_pdgid", &pfcandidate_pdgid_);
  }

  ////////////////////
  // gun particles
  //
  if (storeGunParticles_) {
    t_->Branch("gunparticle_id", &gunparticle_id_);
    t_->Branch("gunparticle_energy", &gunparticle_energy_);
    t_->Branch("gunparticle_pt", &gunparticle_pt_);
    t_->Branch("gunparticle_eta", &gunparticle_eta_);
    t_->Branch("gunparticle_phi", &gunparticle_phi_);
  }

}
PFAnalysis::~PFAnalysis() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//
void PFAnalysis::clearVariables() {
  ev_run_ = 0;
  ev_lumi_ = 0;
  ev_event_ = 0;
  vtx_x_ = 0;
  vtx_y_ = 0;
  vtx_z_ = 0;

  ////////////////////
  // GenParticles
  //
  genpart_eta_.clear();
  genpart_phi_.clear();
  genpart_pt_.clear();
  genpart_energy_.clear();
  genpart_dvx_.clear();
  genpart_dvy_.clear();
  genpart_dvz_.clear();
  genpart_ovx_.clear();
  genpart_ovy_.clear();
  genpart_ovz_.clear();
  genpart_exx_.clear();
  genpart_exy_.clear();
  genpart_mother_.clear();
  genpart_exphi_.clear();
  genpart_exeta_.clear();
  genpart_fbrem_.clear();
  genpart_pid_.clear();
  genpart_gen_.clear();
  genpart_reachedEE_.clear();
  genpart_fromBeamPipe_.clear();
  genpart_posx_.clear();
  genpart_posy_.clear();
  genpart_posz_.clear();

  ////////////////////
  // reco::GenParticles
  //
  gen_eta_.clear();
  gen_phi_.clear();
  gen_pt_.clear();
  gen_energy_.clear();
  gen_charge_.clear();
  gen_pdgid_.clear();
  gen_status_.clear();
  gen_daughters_.clear();

  ////////////////////
  // high purity tracks
  //
  track_eta_.clear();
  track_phi_.clear();
  track_pt_.clear();
  track_energy_.clear();
  track_charge_.clear();
  track_posx_.clear();
  track_posy_.clear();
  track_posz_.clear();

  ////////////////////
  // PF candidates
  //
  pfcandidate_eta_.clear();
  pfcandidate_phi_.clear();
  pfcandidate_pt_.clear();
  pfcandidate_energy_.clear();
  pfcandidate_pdgid_.clear();

  ////////////////////
  // gun particles
  //
  gunparticle_id_.clear();
  gunparticle_energy_.clear();
  gunparticle_pt_.clear();
  gunparticle_eta_.clear();
  gunparticle_phi_.clear();
}

void PFAnalysis::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using namespace edm;

  clearVariables();

  Handle<edm::HepMCProduct> hevH;
  Handle<std::vector<SimTrack>> simTracksHandle;
  Handle<std::vector<SimVertex>> simVerticesHandle;

  iEvent.getByToken(hev_, hevH);

  iEvent.getByToken(simTracks_, simTracksHandle);
  iEvent.getByToken(simVertices_, simVerticesHandle);
  mySimEvent_->fill(*simTracksHandle, *simVerticesHandle);

  Handle<std::vector<reco::Track>> trackHandle;

  iEvent.getByToken(tracks_, trackHandle);
  const std::vector<reco::Track> &tracks = *trackHandle;

  Handle<std::vector<reco::Vertex>> verticesHandle;
  iEvent.getByToken(vertices_, verticesHandle);
  auto const &vertices = *verticesHandle;

  HepMC::GenVertex *primaryVertex = *(hevH)->GetEvent()->vertices_begin();
  float vx_ = primaryVertex->position().x() / 10.;  // to put in official units
  float vy_ = primaryVertex->position().y() / 10.;
  float vz_ = primaryVertex->position().z() / 10.;
  Point sim_pv(vx_, vy_, vz_);

  // fill the gunparticles per vertex
  if (storeGunParticles_) {
    HepMC::GenEvent::vertex_const_iterator vertex_it;
    for (vertex_it = hevH->GetEvent()->vertices_begin();
         vertex_it != hevH->GetEvent()->vertices_end(); vertex_it++) {
      std::vector<int> gunparticle_id;
      std::vector<float> gunparticle_energy;
      std::vector<float> gunparticle_pt;
      std::vector<float> gunparticle_eta;
      std::vector<float> gunparticle_phi;

      HepMC::GenVertex::particles_out_const_iterator particle_it;
      for (particle_it = (*vertex_it)->particles_out_const_begin();
           particle_it != (*vertex_it)->particles_out_const_end(); particle_it++) {
        gunparticle_id.push_back((*particle_it)->pdg_id());
        gunparticle_energy.push_back((*particle_it)->momentum().e());
        gunparticle_pt.push_back((*particle_it)->momentum().perp());
        gunparticle_eta.push_back((*particle_it)->momentum().eta());
        gunparticle_phi.push_back((*particle_it)->momentum().phi());
      }

      gunparticle_id_.push_back(gunparticle_id);
      gunparticle_energy_.push_back(gunparticle_energy);
      gunparticle_pt_.push_back(gunparticle_pt);
      gunparticle_eta_.push_back(gunparticle_eta);
      gunparticle_phi_.push_back(gunparticle_phi);
    }
  } //storeGunParticles

  std::vector<FSimTrack *> allselectedgentracks;
  unsigned int npart = mySimEvent_->nTracks();
  cout << "mySimEvent->nTracks()=" << npart << endl;
  for (unsigned int i = 0; i < npart; ++i) {
    std::vector<float> xp, yp, zp;
    FSimTrack &myTrack(mySimEvent_->track(i));
    math::XYZTLorentzVectorD vtx(0, 0, 0, 0);

    int reachedEE = 0;  // compute the extrapolations for the particles reaching EE
                        // and for the gen particles
    if (!myTrack.noEndVertex()) {
      vtx = myTrack.endVertex().position();
    }
    auto orig_vtx = myTrack.vertex().position();

    allselectedgentracks.push_back(&mySimEvent_->track(i));
    // fill branches
    genpart_eta_.push_back(myTrack.momentum().eta());
    genpart_phi_.push_back(myTrack.momentum().phi());
    genpart_pt_.push_back(myTrack.momentum().pt());
    genpart_energy_.push_back(myTrack.momentum().energy());
    genpart_dvx_.push_back(vtx.x());
    genpart_dvy_.push_back(vtx.y());
    genpart_dvz_.push_back(vtx.z());

    genpart_ovx_.push_back(orig_vtx.x());
    genpart_ovy_.push_back(orig_vtx.y());
    genpart_ovz_.push_back(orig_vtx.z());

    genpart_pid_.push_back(myTrack.type());
    genpart_gen_.push_back(myTrack.genpartIndex());
    genpart_reachedEE_.push_back(reachedEE);
    genpart_fromBeamPipe_.push_back(true);

    genpart_posx_.push_back(xp);
    genpart_posy_.push_back(yp);
    genpart_posz_.push_back(zp);
  } //tracks from geant

  // associate gen particles to mothers
  genpart_mother_.resize(genpart_posz_.size(), -1);
  for (size_t i = 0; i < allselectedgentracks.size(); i++) {
    const auto tracki = allselectedgentracks.at(i);

    for (size_t j = i + 1; j < allselectedgentracks.size(); j++) {
      const auto trackj = allselectedgentracks.at(j);

      if (!tracki->noMother()) {
        if (&tracki->mother() == trackj) genpart_mother_.at(i) = j;
      }
      if (!trackj->noMother()) {
        if (&trackj->mother() == tracki) genpart_mother_.at(j) = i;
      }
    }
  }

  if (readGen_) {
    Handle<std::vector<reco::GenParticle>> genParticlesHandle;
    iEvent.getByToken(genParticles_, genParticlesHandle);
    for (std::vector<reco::GenParticle>::const_iterator it_p = genParticlesHandle->begin();
         it_p != genParticlesHandle->end(); ++it_p) {
      gen_eta_.push_back(it_p->eta());
      gen_phi_.push_back(it_p->phi());
      gen_pt_.push_back(it_p->pt());
      gen_energy_.push_back(it_p->energy());
      gen_charge_.push_back(it_p->charge());
      gen_pdgid_.push_back(it_p->pdgId());
      gen_status_.push_back(it_p->status());
      std::vector<int> daughters(it_p->daughterRefVector().size(), 0);
      for (unsigned j = 0; j < it_p->daughterRefVector().size(); ++j) {
        daughters[j] = static_cast<int>(it_p->daughterRefVector().at(j).key());
      }
      gen_daughters_.push_back(daughters);
    }
  }

  // loop over tracks

  // prepare for RK propagation
  defaultRKPropagator::Product prod(aField_, alongMomentum, 5.e-5);
  auto &RKProp = prod.propagator;

  for (std::vector<reco::Track>::const_iterator it_track = tracks.begin(); it_track != tracks.end();
       ++it_track) {
    if (!it_track->quality(reco::Track::highPurity)) continue;

    double energy = it_track->pt() * cosh(it_track->eta());

    // save info in tree
    track_pt_.push_back(it_track->pt());
    track_eta_.push_back(it_track->eta());
    track_phi_.push_back(it_track->phi());
    track_energy_.push_back(energy);
    track_charge_.push_back(it_track->charge());
  }  // end loop over tracks

  if (storePFCandidates_) {
    Handle<std::vector<reco::PFCandidate>> pfCandidatesHandle;
    iEvent.getByToken(pfCandidates_, pfCandidatesHandle);
    cout << "pfCandidates=" << pfCandidatesHandle->size() << endl;
    int npi = 0;
    for (std::vector<reco::PFCandidate>::const_iterator it_p = pfCandidatesHandle->begin();
         it_p != pfCandidatesHandle->end(); ++it_p) {
      pfcandidate_eta_.push_back(it_p->eta());
      pfcandidate_phi_.push_back(it_p->phi());
      pfcandidate_pt_.push_back(it_p->pt());
      pfcandidate_energy_.push_back(it_p->energy());
      pfcandidate_pdgid_.push_back(it_p->pdgId());
      if (abs(it_p->pdgId())==211 && it_p->trackRef().isNonnull()) {
          npi += 1;
      }
    }
    cout << "pfCandidates(pi)=" << npi << endl;
  }

  ev_event_ = iEvent.id().event();
  ev_lumi_ = iEvent.id().luminosityBlock();
  ev_run_ = iEvent.id().run();

  vtx_x_ = vx_;
  vtx_y_ = vy_;
  vtx_z_ = vz_;

  t_->Fill();
}

void PFAnalysis::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
  edm::ESHandle<HepPDT::ParticleDataTable> pdt;
  es.getData(pdt);
  mySimEvent_->initializePdt(&(*pdt));

  edm::ESHandle<MagneticField> magfield;
  es.get<IdealMagneticFieldRecord>().get(magfield);

  aField_ = &(*magfield);
}

void PFAnalysis::endRun(edm::Run const &iEvent, edm::EventSetup const &) {}

void PFAnalysis::beginJob() { ; }

// ------------ method called once each job just after ending the event loop
// ------------
void PFAnalysis::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the
// module  ------------

void PFAnalysis::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  // The following says we do not know what parameters are allowed so do no
  // validation
  // Please change this to state exactly what you do use, even if it is no
  // parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(PFAnalysis);
