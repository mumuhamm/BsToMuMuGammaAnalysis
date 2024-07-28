#ifndef BsToMuMuGammaAnalysis_RadiativeAnalysis_RadiativeRootTree
#define BsToMuMuGammaAnalysis_RadiativeAnalysis_RadiativeRootTree

#include <string>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <vector>

class RadiativeRootTree {
	public: 
		RadiativeRootTree();
		~RadiativeRootTree();
		void resetEntries();
		void writeFile();
		void createTree(const std::string filename);
		void readTree(const std::string filename);
		void readTree(std::vector<std::string> filenames);
		void getAngles(const double aa, const double bb, const double cc, const double dd);
		void getVtx(const double aa, const double bb, const double cc, const double dd, const double ee, const double ff, const double gg, const double hh, const double ii);
		void getDeDx(const double f1, const double f2, const int f3);
		void fill();
		void setBranchAddresses();

	public: 
		 int    runNumber_;
                 int    PUinteraction_;
		 int    PUTrueinteraction_;
                 unsigned int eventNumber_;
                 int    lumiSection_;
		 int    isPV_;
		 int    NVerticesbeforecut_;
		 int    NVerticesaftercut_;
		 int    isBS_;
		 double BSx_ ;
		 double BSy_ ;
		 double BSz_ ;
		 double BSdx_ ;
		 double BSdy_ ;
		 double BSdz_ ;
		 double BSdydz_;
		 double BSdxdz_;
		 double BSsigmaZ_;
		 double BSdsigmaZ_;
		 double PVx_ ;
		 double PVy_ ;
		 double PVz_ ;
		 double PVerrx_ ;
		 double PVerry_ ;
		 double PVerrz_ ;
		 double PVndof_;
		 double PVrho_;
		 double dedxTrk_;
		 double errdedxTrk_;
		 int    numdedxTrk_;
		 int    photonMultiplicity_;
                 double photonPt_;
                 double photonEta_;
                 double photonPhi_;
		 double photonTrkIso_;
		 double photonEcalIso_;
		 double photonHcalIso_;
		 double photonCaloIso_;
		 double photonSigmaiEtaiEta_;
		 double photonSigmaiEtaiPhi_;
		 double photonSigmaiPhiiPhi_;
		 double photonSigmaEtaEta_;
		 double photone1x5_;
		 double photone2x5_;
		 double photone3x3_;
		 double photone5x5_;
		 double photonhcalDepth1OverEcal_;
		 double photonhcalDepth2OverEcal_;
		 double photonhcalDepth1OverEcalBc_;
		 double photonhcalDepth2OverEcalBc_;
		 int    electronMultiplicity_;
		 int    pfCandMultiplicity_;

		 double costheta_;
		 double phi_;
		 double cospsi_;
		 double AngleBsDecayLength_;



		 //arrays 
		 double PVcovariance_[9];


		TFile* bmmgFile_;
                TTree* bmmgTree_;
};
#endif


