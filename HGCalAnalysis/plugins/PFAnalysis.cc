// Based on RecoNtuple/HGCalAnalysis with modifications for PF
//
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
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

#include "DataFormats/Math/interface/deltaPhi.h"
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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TH1F.h"
#include "TTree.h"

#include <map>
#include <set>
#include <string>
#include <vector>
#include <map>

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

  // ----------member data ---------------------------

  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_;
  edm::EDGetTokenT<edm::View<TrackingParticle>> trackingParticles_;
  edm::EDGetTokenT<edm::View<reco::Track>> tracks_;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandidates_;
  edm::EDGetTokenT<reco::RecoToSimCollection> tracks_recotosim_;

  TTree *t_;

  edm::RunNumber_t ev_run_;
  edm::LuminosityBlockNumber_t ev_lumi_;
  edm::EventNumber_t ev_event_;

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
  std::vector<int> genpart_pid_;
  std::vector<int> genpart_gen_;

  std::vector<float> gen_eta_;
  std::vector<float> gen_phi_;
  std::vector<float> gen_pt_;
  std::vector<float> gen_energy_;
  std::vector<int> gen_charge_;
  std::vector<int> gen_pdgid_;
  std::vector<int> gen_status_;
  std::vector<std::vector<int>> gen_daughters_;

  std::vector<float> track_eta_;
  std::vector<float> track_phi_;
  std::vector<float> track_pt_;
  std::vector<float> track_energy_;
  std::vector<int> track_charge_;
  std::vector<std::vector<float>> track_posx_;
  std::vector<std::vector<float>> track_posy_;
  std::vector<std::vector<float>> track_posz_;

  std::vector<float> pfcandidate_eta_;
  std::vector<float> pfcandidate_phi_;
  std::vector<float> pfcandidate_pt_;
  std::vector<float> pfcandidate_energy_;
  std::vector<int> pfcandidate_pdgid_;
  
  // and also the magnetic field
  MagneticField const *aField_;
};

PFAnalysis::PFAnalysis() { ; }

PFAnalysis::PFAnalysis(const edm::ParameterSet &iConfig) {

  tracks_recotosim_ = consumes<reco::RecoToSimCollection>(edm::InputTag("trackingParticleRecoTrackAsssociation"));
  trackingParticles_ = consumes<edm::View<TrackingParticle>>(edm::InputTag("mix", "MergedTrackTruth"));
  genParticles_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
  pfCandidates_ = consumes<std::vector<reco::PFCandidate>>(edm::InputTag("particleFlow"));
  tracks_ = consumes<edm::View<reco::Track>>(edm::InputTag("generalTracks"));

  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  fs->make<TH1F>("total", "total", 100, 0, 5.);

  t_ = fs->make<TTree>("pftree", "pftree");

  // event info
  t_->Branch("event", &ev_event_);
  t_->Branch("lumi", &ev_lumi_);
  t_->Branch("run", &ev_run_);

  t_->Branch("genpart_eta", &genpart_eta_);
  t_->Branch("genpart_phi", &genpart_phi_);
  t_->Branch("genpart_pt", &genpart_pt_);
  t_->Branch("genpart_energy", &genpart_energy_);
  t_->Branch("genpart_dvx", &genpart_dvx_);
  t_->Branch("genpart_dvy", &genpart_dvy_);
  t_->Branch("genpart_dvz", &genpart_dvz_);
  t_->Branch("genpart_pid", &genpart_pid_);
  t_->Branch("genpart_gen", &genpart_gen_);

  t_->Branch("gen_eta", &gen_eta_);
  t_->Branch("gen_phi", &gen_phi_);
  t_->Branch("gen_pt", &gen_pt_);
  t_->Branch("gen_energy", &gen_energy_);
  t_->Branch("gen_charge", &gen_charge_);
  t_->Branch("gen_pdgid", &gen_pdgid_);
  t_->Branch("gen_status", &gen_status_);
  t_->Branch("gen_daughters", &gen_daughters_);

  t_->Branch("track_eta", &track_eta_);
  t_->Branch("track_phi", &track_phi_);
  t_->Branch("track_pt", &track_pt_);
  t_->Branch("track_energy", &track_energy_);
  t_->Branch("track_charge", &track_charge_);
  t_->Branch("track_posx", &track_posx_);
  t_->Branch("track_posy", &track_posy_);
  t_->Branch("track_posz", &track_posz_);

  t_->Branch("pfcandidate_eta", &pfcandidate_eta_);
  t_->Branch("pfcandidate_phi", &pfcandidate_phi_);
  t_->Branch("pfcandidate_pt", &pfcandidate_pt_);
  t_->Branch("pfcandidate_energy", &pfcandidate_energy_);
  t_->Branch("pfcandidate_pdgid", &pfcandidate_pdgid_);
} // constructor

PFAnalysis::~PFAnalysis() {
}

void PFAnalysis::clearVariables() {
  ev_run_ = 0;
  ev_lumi_ = 0;
  ev_event_ = 0;

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
  genpart_pid_.clear();
  genpart_gen_.clear();

  gen_eta_.clear();
  gen_phi_.clear();
  gen_pt_.clear();
  gen_energy_.clear();
  gen_charge_.clear();
  gen_pdgid_.clear();
  gen_status_.clear();
  gen_daughters_.clear();

  track_eta_.clear();
  track_phi_.clear();
  track_pt_.clear();
  track_energy_.clear();
  track_charge_.clear();
  track_posx_.clear();
  track_posy_.clear();
  track_posz_.clear();

  pfcandidate_eta_.clear();
  pfcandidate_phi_.clear();
  pfcandidate_pt_.clear();
  pfcandidate_energy_.clear();
  pfcandidate_pdgid_.clear();
} //clearVariables

void PFAnalysis::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  clearVariables();

  //Simulated tracks, cleaned up by TrackingTruthAccumulator
  edm::Handle<edm::View<TrackingParticle>> trackingParticlesHandle;
  iEvent.getByToken(trackingParticles_, trackingParticlesHandle);
  const edm::View<TrackingParticle> &trackingParticles = *trackingParticlesHandle;

  //Matches reco tracks to sim tracks (TrackingParticle)
  edm::Handle<reco::RecoToSimCollection> recotosimCollection;
  iEvent.getByToken(tracks_recotosim_, recotosimCollection);
  const auto recotosim = *recotosimCollection;

  edm::Handle<edm::View<reco::Track>> trackHandle;
  iEvent.getByToken(tracks_, trackHandle);
  const edm::View<reco::Track> &tracks = *trackHandle;

  //Keep track of sim tracks that were matched to a good, high-purity reco track
  std::map<long, bool> trackingparticle_matched_to_recotrack;

  // loop over high-purity reco tracks
  for (unsigned long ntrack = 0; ntrack < tracks.size(); ntrack++) {
    const auto &track = tracks.at(ntrack);
    if (!track.quality(reco::Track::highPurity)) {
      continue;
    }

    edm::RefToBase<reco::Track> trackref(trackHandle, ntrack);

    if (recotosim.find(trackref) != recotosim.end()) {
      const auto &tps = recotosim[trackref];
      for (const auto tp : tps) {
        edm::Ref<std::vector<TrackingParticle>> tpr = tp.first;
        //double associationQuality = tp.second;

        //make a note that this trackingparticle was matched to a high-purity reco track
        trackingparticle_matched_to_recotrack[tpr.key()] = true;
      }
    }

    double energy = track.pt() * cosh(track.eta());

    // save info in tree
    track_pt_.push_back(track.pt());
    track_eta_.push_back(track.eta());
    track_phi_.push_back(track.phi());
    track_energy_.push_back(energy);
    track_charge_.push_back(track.charge());
  }  // end loop over tracks

  for (unsigned long ntrackingparticle = 0; ntrackingparticle < trackingParticles.size(); ntrackingparticle++) {
    const auto &tp = trackingParticles.at(ntrackingparticle);
    edm::RefToBase<TrackingParticle> tpref(trackingParticlesHandle, ntrackingparticle);

    math::XYZTLorentzVectorD vtx(0, 0, 0, 0);

    //trackingparticle was not matched to a reco track, skip it
    if (trackingparticle_matched_to_recotrack.find(tpref.key()) == trackingparticle_matched_to_recotrack.end()) {
      continue;
    }

    if (tp.decayVertices().size() > 0) {
      vtx = tp.decayVertices().at(0)->position();
    }
    auto orig_vtx = tp.vertex();

    // fill branches
    genpart_eta_.push_back(tp.p4().eta());
    genpart_phi_.push_back(tp.p4().phi());
    genpart_pt_.push_back(tp.p4().pt());
    genpart_energy_.push_back(tp.p4().energy());
    genpart_dvx_.push_back(vtx.x());
    genpart_dvy_.push_back(vtx.y());
    genpart_dvz_.push_back(vtx.z());

    genpart_ovx_.push_back(orig_vtx.x());
    genpart_ovy_.push_back(orig_vtx.y());
    genpart_ovz_.push_back(orig_vtx.z());

    genpart_pid_.push_back(tp.pdgId());
    genpart_gen_.push_back(0);
  }

  edm::Handle<std::vector<reco::GenParticle>> genParticlesHandle;
  iEvent.getByToken(genParticles_, genParticlesHandle);
  for (std::vector<reco::GenParticle>::const_iterator it_p = genParticlesHandle->begin();
       it_p != genParticlesHandle->end();
       ++it_p) {
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

  edm::Handle<std::vector<reco::PFCandidate>> pfCandidatesHandle;
  iEvent.getByToken(pfCandidates_, pfCandidatesHandle);
  for (std::vector<reco::PFCandidate>::const_iterator it_p = pfCandidatesHandle->begin();
       it_p != pfCandidatesHandle->end();
       ++it_p) {
    pfcandidate_eta_.push_back(it_p->eta());
    pfcandidate_phi_.push_back(it_p->phi());
    pfcandidate_pt_.push_back(it_p->pt());
    pfcandidate_energy_.push_back(it_p->energy());
    pfcandidate_pdgid_.push_back(it_p->pdgId());
  }

  ev_event_ = iEvent.id().event();
  ev_lumi_ = iEvent.id().luminosityBlock();
  ev_run_ = iEvent.id().run();

  t_->Fill();
} //analyze

void PFAnalysis::beginRun(edm::Run const &iEvent, edm::EventSetup const &es) {
  edm::ESHandle<HepPDT::ParticleDataTable> pdt;
  es.getData(pdt);

  edm::ESHandle<MagneticField> magfield;
  es.get<IdealMagneticFieldRecord>().get(magfield);

  aField_ = &(*magfield);
}

void PFAnalysis::endRun(edm::Run const &iEvent, edm::EventSetup const &) {}

void PFAnalysis::beginJob() { ; }

void PFAnalysis::endJob() {}

void PFAnalysis::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PFAnalysis);
