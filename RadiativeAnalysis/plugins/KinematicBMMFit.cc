#include "BsToMuMuGammaAnalysis/RadiativeAnalysis/interface/KinematicBMMFit.h"

#include <TMath.h>
KinematicBMMFit::KinematicBMMFit(){

}

bool KinematicBMMFit::doFit(std::vector<reco::TransientTrack> t_tracks, const double muonMass, const double mass1, const double  mass2){
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
