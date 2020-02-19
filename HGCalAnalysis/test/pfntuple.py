import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process("Demo")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase1_2021_realistic", "")
from FastSimulation.Event.ParticleFilter_cfi import *

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file://step3_RAW2DIGI_L1Reco_RECO_RECOSIM_phase1.root'
    ),
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck")
)

process.ana = cms.EDAnalyzer('PFAnalysis',
    verbose = cms.bool(True),
    storeGenParticleOrigin = cms.bool(True),
    storeGenParticleExtrapolation = cms.bool(True),
    storePFCandidates = cms.bool(True),
    storeGunParticles = cms.bool(True),
    readGenParticles = cms.bool(True),
    TestParticleFilter = ParticleFilterBlock.ParticleFilter
)

process.ana.TestParticleFilter.protonEMin = cms.double(100000)
process.ana.TestParticleFilter.etaMax = cms.double(3.1)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("pfntuple.root")
)
process.p = cms.Path(
  process.tpClusterProducer*
  process.quickTrackAssociatorByHits*
  process.trackingParticleRecoTrackAsssociation*
  process.ana)
