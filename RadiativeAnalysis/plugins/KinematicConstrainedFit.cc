#include "BsToMuMuGammaAnalysis/RadiativeAnalysis/interface/KinematicConstrainedFit.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
using namespace reco;
using namespace edm;
using namespace std;
using namespace pat;
#include <TMath.h>
KinematicConstrainedFit::KinematicConstrainedFit(){

}

bool KinematicConstrainedFit::doFit(std::vector<reco::TransientTrack> t_tracks, const double muonMass, const double mass1, const double  mass2){
	reco::TransientTrack track_MuP = t_tracks[0];
	reco::TransientTrack track_MuM = t_tracks[1];
	//Creating a KinematicParticleFactory
	KinematicParticleFactoryFromTransientTrack pFactory;
	//The mass of a muon and the insignificant mass sigma to avoid singularities in the covariance matrix.
        float muon_sigma = 0.0000000001;
	//initial chi2 and ndf before kinematic fits. The chi2 of the reconstruction is not considered
	float chi = 0.;
	float ndf = 0.;
	std::vector<RefCountedKinematicParticle> allParticlesMu;
	allParticlesMu.push_back(pFactory.particle (track_MuP, muonMass, chi, ndf, muon_sigma));
	allParticlesMu.push_back(pFactory.particle (track_MuM, muonMass, chi, ndf, muon_sigma));
	KinematicParticleVertexFitter Fitter;
        std::vector<RefCountedKinematicParticle> allParticlesTrk;
	RefCountedKinematicTree BsMMTree = Fitter.fit(allParticlesMu);
	//if the fit fails, do not consider this as candidate
	if(BsMMTree->isEmpty()) return 0;
	KinematicParticleFitter constFitter;
	double nominalBsMass = 5.36689;
	double bsMSigma = 0.00019;
	KinematicConstraint * bsmm_const = new MassKinematicConstraint(nominalBsMass, bsMSigma);
        BsMMTree = constFitter.fit(bsmm_const,BsMMTree);
        myTree_BsMM= BsMMTree;
        if(BsMMTree->isEmpty()) {
		delete bsmm_const;
		return 0;
	} 
	BsMMTree->movePointerToTheTop();
	RefCountedKinematicParticle BsMM_branch = BsMMTree->currentParticle();
	allParticlesTrk.push_back(BsMM_branch);
	myTree_Bs = Fitter.fit(allParticlesTrk);
	if(myTree_Bs->isEmpty()) {
		delete bsmm_const;
		return 0;
       	}
	myTree_Bs->movePointerToTheTop();
	bsmmg = myTree_Bs->currentParticle();
        bsVertex = myTree_Bs->currentDecayVertex();
	if (!bsVertex->vertexIsValid()) {
        delete bsmm_const;
        return 0;
       }
	vtxprob_Bs = TMath::Prob(bsmmg->chiSquared(), (int)bsmmg->degreesOfFreedom());
	delete bsmm_const;
        return 1;
}
bool KinematicConstrainedFit::dobsphikkgFit(std::vector<reco::TransientTrack> t_tracks, const double mass1, const double  mass2, const pat::Photon &photon){
	reco::TransientTrack track_KP = t_tracks[0];
	reco::TransientTrack track_KM = t_tracks[1];
	KinematicParticleFactoryFromTransientTrack pFactory;
	float kaon_sigma = 0.0000000001;
	float chi = 0.;
	float ndf = 0.;
	std::vector<RefCountedKinematicParticle> allParticlesK;
	allParticlesK.push_back(pFactory.particle (track_KP, mass1, chi, ndf, kaon_sigma));
	allParticlesK.push_back(pFactory.particle (track_KM, mass1, chi, ndf, kaon_sigma));
	KinematicParticleVertexFitter Fitter;
	RefCountedKinematicTree PhiTree = Fitter.fit(allParticlesK);
	if(PhiTree->isEmpty()) return 0;
	KinematicParticleFitter constFitter;
	double nominalPhiMass =  1.019455;
	float phiMsigma = 0.00002;
	KinematicConstraint * phi_const = new MassKinematicConstraint( nominalPhiMass, phiMsigma);
	PhiTree = constFitter.fit(phi_const,PhiTree);
	renewed_BsConstrainedTree = PhiTree;
	if(renewed_BsConstrainedTree->isEmpty()) {
		delete phi_const;
		return 0;
	}
	renewed_BsConstrainedTree->movePointerToTheTop();
	RefCountedKinematicParticle phi = renewed_BsConstrainedTree->currentParticle();

	// Incorporate the photon's four-momentum
    TLorentzVector photonP4;
    photonP4.SetPtEtaPhiE(photon.pt(), photon.eta(), photon.phi(), photon.energy());
    // Create a dummy particle for the photon (mass = 0)
    KinematicParticleFactoryFromTransientTrack pFactoryPhotons;
    RefCountedKinematicParticle photonParticle = pFactoryPhotons.particle(photonP4, 0.0, chi, ndf, 0.0);
    
    // Combine the phi meson and the photon to form the Bs candidate
    std::vector<RefCountedKinematicParticle> allParticlesBs;
    allParticlesBs.push_back(phi);
    allParticlesBs.push_back(photonParticle);
	RefCountedKinematicTree BsTree = Fitter.fit(allParticlesBs);
    
    if (BsTree->isEmpty()) {
        delete phi_const;
        return false;
    }
    
    BsTree->movePointerToTheTop();
	bs = BsTree->currentParticle();
    bsVertex = BsTree->currentDecayVertex();
    
    if (!bsVertex->vertexIsValid()) {
        delete phi_const;
        return false;
    }
    
    vtxprob_Bs = TMath::Prob(bs->chiSquared(), (int)bs->degreesOfFreedom());
    
    delete phi_const;
    return true;
	
	
}

bool KinematicConstrainedFit::dobsphimmgFit(std::vector<reco::TransientTrack> t_tracks, const double muonMass){
	reco::TransientTrack track_MuP = t_tracks[0];
	reco::TransientTrack track_MuM = t_tracks[1];
	KinematicParticleFactoryFromTransientTrack pFactory;
	float muon_sigma =  0.0000000001;
	float chi =0.;
	float ndf =0.;
	std::vector<RefCountedKinematicParticle> allParticlesMu;
	allParticlesMu.push_back(pFactory.particle (track_MuP, muonMass, chi, ndf, muon_sigma));
	allParticlesMu.push_back(pFactory.particle (track_MuM, muonMass, chi, ndf, muon_sigma));
	KinematicParticleVertexFitter Fitter;
	RefCountedKinematicTree PhiTree = Fitter.fit(allParticlesMu);
	if(PhiTree->isEmpty()) return 0;
	KinematicParticleFitter constFitter;
	double nominalPhiMass =  1.019455;
	float phiMsigma = 0.00002;
	KinematicConstraint * phi_const = new MassKinematicConstraint( nominalPhiMass, phiMsigma);
	PhiTree = constFitter.fit(phi_const,PhiTree);
	renewed_BsConstrainedTree = PhiTree;
	if(renewed_BsConstrainedTree->isEmpty()) {
		delete phi_const;
		return 0;
	}
	renewed_BsConstrainedTree->movePointerToTheTop();
	bs = renewed_BsConstrainedTree->currentParticle();
	bsVertex = renewed_BsConstrainedTree->currentDecayVertex();
	if (!bsVertex->vertexIsValid()) {
		delete phi_const;
		return 0;
	}
	vtxprob_Bs = TMath::Prob(bs->chiSquared(), (int)bs->degreesOfFreedom());
	delete phi_const;
	return 1;
}
