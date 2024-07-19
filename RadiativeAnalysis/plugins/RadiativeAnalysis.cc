// -*- C++ -*-
// Package:    BsToMuMuGammaAnalysis/RadiativeAnalysis
// Class:      RadiativeAnalysis
// Original Author:  Alibordi Muhammad
//         Created:  Fri, 12 Jul 2024 10:27:36 GMT
//============================================================
#include "BsToMuMuGammaAnalysis/RadiativeAnalysis/interface/KinematicBMMFit.h"
#include "BsToMuMuGammaAnalysis/RadiativeAnalysis/interface/RadiativeRootTree.h"
#include "BsToMuMuGammaAnalysis/RadiativeAnalysis/interface/RadiativeAnalysis.h"


#include <memory>
#include <cstddef>
#include <cfloat>
#include <string>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TVector.h"
#include "TLorentzRotation.h"
#include <iostream>
#include <TMath.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "BHadron/JpsiTrkTrk/interface/JpsiTrkTrk.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Association.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTrackLinks.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/CaloMuon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/OverlapChecker.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerAlgorithm.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L2MuonTrajectorySeedCollection.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"
#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "MuonAnalysis/MuonAssociators/interface/L1MuonMatcherAlgo.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

using namespace reco;
using namespace edm;
using namespace std;
using namespace pat;


RadiativeAnalysis::RadiativeAnalysis(const edm::ParameterSet& iConfig): 
	theConfig_(iConfig),
	nominalJpsiMass( 3.096916 ),
	nominalPsiMass( 3.096916 ),
	nominalPhiMass(1.019 ),
	nominalElectronMass(0.00051099893),
	nominalMuonMass(0.1056583)
{
	isMCstudy_ = iConfig.getParameter<bool>("isMCstudy");
	genParticlesLabel                 = iConfig.getParameter<InputTag>("genParticlesLabel");
        genParticlesTok                   = consumes<edm::View<reco::GenParticle>>(genParticlesLabel);
        MuonTag                           = iConfig.getParameter<edm::InputTag>("MuonTag");
        MuonTagTok                        = consumes<edm::View<pat::Muon>>(MuonTag);
	JetTag                            = iConfig.getParameter<edm::InputTag>("JetTag");
        JetTagTok                         = consumes<edm::View<pat::Jet>>(JetTag);
	PhotonTag                         = iConfig.getParameter<edm::InputTag>("PhotonTag");
        PhotonTagTok                      = consumes<edm::View<pat::Photon>>(PhotonTag);
	OOTPhotonTag                      = iConfig.getParameter<edm::InputTag>("OOTPhotonTag");
        OOTPhotonTagTok                   = consumes<edm::View<pat::Photon>>(OOTPhotonTag);
	ElectronTag                       = iConfig.getParameter<edm::InputTag>("ElectronTag");
        ElectronTagTok                    = consumes<edm::View<pat::Electron>>(ElectronTag);
	SuperClusterTag                   = iConfig.getParameter<edm::InputTag>("SuperClusterTag");
	SuperClusterTagTok                = consumes<edm::View<reco::SuperCluster>>(SuperClusterTag);
	OOTSuperClusterTag                = iConfig.getParameter<edm::InputTag>("OOTSuperClusterTag");
        OOTSuperClusterTagTok             = consumes<edm::View<reco::SuperCluster>>(OOTSuperClusterTag);
	PUInfo                            = iConfig.getParameter<InputTag>("PUInfo");
        PUInfoTok                         = consumes<edm::View<PileupSummaryInfo>>(PUInfo);
        vertexBeamSpot                    = iConfig.getParameter<edm::InputTag>("vertexBeamSpot");
        vertexBeamSpotTok                 = consumes<reco::BeamSpot>(vertexBeamSpot);
	primaryvertex                     = iConfig.getParameter<edm::InputTag>("primaryvertex");
        primaryvertexTok                  = consumes<edm::View<reco::Vertex>>(primaryvertex);
	triggerresults                    = iConfig.getParameter<edm::InputTag>("triggerresults");
        triggerresultsTok                 = consumes<edm::TriggerResults>(triggerresults);
	pfCandTag                         = iConfig.getParameter<edm::InputTag>("pfCandTag");
        pfCandTagTok                      = consumes<edm::View<pat::PackedCandidate>>(pfCandTok);
	IsoTrackTag                       = iConfig.getParameter<edm::InputTag>("IsoTrackTag");
        IsoTrackTagTok                    = consumes<edm::View<pat::IsolatedTrack>>(IsoTrackTag);
        trackBuilderTok                   = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));


	StoreDeDxInfo_                    = iConfig.getParameter<bool>("StoreDeDxInfo");
        JpsiMassWindowBeforeFit_          = iConfig.getParameter<double>("JpsiMassWindowBeforeFit");
        JpsiMassWindowAfterFit_           = iConfig.getParameter<double>("JpsiMassWindowAfterFit");
        JpsiPtCut_                        = iConfig.getParameter<double>("JpsiPtCut");
	PsiMassWindowBeforeFit_           = iConfig.getParameter<double>("PsiMassWindowBeforeFit");
        PsiMassWindowAfterFit_            = iConfig.getParameter<double>("PsiMassWindowAfterFit");
	BsLowerMassCutBeforeFit_          = iConfig.getParameter<double>("BsLowerMassCutBeforeFit");
        BsUpperMassCutBeforeFit_          = iConfig.getParameter<double>("BsUpperMassCutBeforeFit");
        BsLowerMassCutAfterFit_           = iConfig.getParameter<double>("BsLowerMassCutAfterFit");
        BsUpperMassCutAfterFit_           = iConfig.getParameter<double>("BsUpperMassCutAfterFit");
	BdPDGMass_                        = iConfig.getParameter<double>("BdPDGMass");
        BpPDGMass_                        = iConfig.getParameter<double>("BpPDGMass");
        BsPDGMass_                        = iConfig.getParameter<double>("BsPDGMass");
	outputFile_                       = iConfig.getUntrackedParameter<std::string>("outputFile");
	verbose_                          = iConfig.getParameter<bool>("verbose");
        TestVerbose_                      = iConfig.getParameter<bool>("TestVerbose");
	
	event_counter_ = 0;
        elecounter_    = 0;
        muoncounter_   = 0;
        jetcounter_    = 0;
        tagmucounter_  = 0;
	photoncounter_ = 0;


}

RadiativeAnalysis::~RadiativeAnalysis() {}
void RadiativeAnalysis::beginJob() {
  bmmgRootTree_ = new RadiativeRootTree();
  bmmgRootTree_->createTree(outputFile_);
}

// ------------ method called for each event  ------------
void RadiativeAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}
void RadiativeAnalysis::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void RadiativeAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RadiativeAnalysis);
