#ifndef BsToMuMuGammaAnalysis_RadiativeAnalysis_KinematicConstrainedFit_h
#define BsToMuMuGammaAnalysis_RadiativeAnalysis_KinematicConstrainedFit_h


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"

#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

class KinematicConstrainedFit{
	public: 
		KinematicConstrainedFit();
		~KinematicConstrainedFit(){}
		bool doFit(std::vector<reco::TransientTrack> t_tracks, const double muonMass, const  double mass1, const double  mass2);
		bool dobsphikkgFit(std::vector<reco::TransientTrack> t_tracks, const  double mass1, const double  mass2);
		bool dobsphimmgFit(std::vector<reco::TransientTrack> t_tracks, const double muonMass);
		double getProb() { return vtxprob_Bs; }
		RefCountedKinematicParticle getParticle() {return bsmmg; }
                RefCountedKinematicVertex  getVertex()   { return bsVertex; }
                RefCountedKinematicTree    getTree()     { return myTree_Bs; }
                RefCountedKinematicTree    getBsMMTree()     { return myTree_BsMM; }

	private: 
		double vtxprob_Bs;
		RefCountedKinematicParticle bsmmg;
		RefCountedKinematicParticle bs;
		RefCountedKinematicVertex bsVertex;
		RefCountedKinematicTree myTree_Bs;
		RefCountedKinematicTree myTree_BsMM;
		RefCountedKinematicTree renewed_BsConstrainedTree;

};
#endif
