#include "BsToMuMuGammaAnalysis/RadiativeAnalysis/interface/RadiativeRootTree.h"
using namespace std;

RadiativeRootTree::RadiativeRootTree()
	{
		resetEntries();
		bmmgTree_ = 0;
		bmmgFile_ = 0;
	}

void RadiativeRootTree::createTree(const std::string filename)
{
  // open root file
  bmmgFile_ = new TFile (filename.c_str(), "RECREATE" );
  int bufsize = 256000;
  // create tree structure
  bmmgTree_ = new TTree("BMMGTree","BMMGTree",bufsize);

  bmmgTree_->Branch("runNumber",&runNumber_ ,"runNumber/I");
  bmmgTree_->Branch("eventNumber",&eventNumber_,"eventNumber/i");
  bmmgTree_->Branch("lumiSection",&lumiSection_,"lumiSection/I");
  bmmgTree_->Branch("PUinteraction",&PUinteraction_,"PUinteraction/I");
  bmmgTree_->Branch("PUTrueinteraction",&PUTrueinteraction_,"PUTrueinteraction/I");
  bmmgTree_->Branch("isPV",&isPV_,"isPV/I");
  bmmgTree_->Branch("NVerticesbeforecut",&NVerticesbeforecut_,"NVerticesbeforecut/I");
  bmmgTree_->Branch("NVerticesaftercut",&NVerticesaftercut_,"NVerticesaftercut/I");
  bmmgTree_->Branch("PVx",&PVx_,"PVx/D");
  bmmgTree_->Branch("PVy",&PVy_,"PVy/D");
  bmmgTree_->Branch("PVz",&PVz_,"PVz/D");
  bmmgTree_->Branch("PVerrx",&PVerrx_,"PVerrx/D");
  bmmgTree_->Branch("PVerry",&PVerry_,"PVerry/D");
  bmmgTree_->Branch("PVerrz",&PVerrz_,"PVerrz/D");
  bmmgTree_->Branch("PVcovariance",&PVcovariance_, "PVcovariance[9]/D");
  bmmgTree_->Branch("PVndof",&PVndof_,"PVndof/D");
  bmmgTree_->Branch("PVrho",&PVrho_,"PVrho/D");
  bmmgTree_->Branch("isBS",&isBS_,"isBS/I");
  bmmgTree_->Branch("BSx",&BSx_,"BSx/D");
  bmmgTree_->Branch("BSy",&BSy_,"BSy/D");
  bmmgTree_->Branch("BSz",&BSz_,"BSz/D");
  bmmgTree_->Branch("BSdx",&BSdx_,"BSdx/D");
  bmmgTree_->Branch("BSdy",&BSdy_,"BSdy/D");
  bmmgTree_->Branch("BSdz",&BSdz_,"BSdz/D");
  bmmgTree_->Branch("BSdxdz",&BSdxdz_,"BSdxdz/D");
  bmmgTree_->Branch("BSdydz",&BSdydz_,"BSdydz/D");
  bmmgTree_->Branch("BSsigmaZ",&BSsigmaZ_,"BSsigmaZ/D");
  bmmgTree_->Branch("BSdsigmaZ",&BSdsigmaZ_,"BSdsigmaZ/D");
  bmmgTree_->Branch("dedxTrk",&dedxTrk_,"dedxTrk/D");
  bmmgTree_->Branch("errdedxTrk",&errdedxTrk_,"errdedxTrk/D");
  bmmgTree_->Branch("numdedxTrk",&numdedxTrk_,"numdedxTrk/I");
  bmmgTree_->Branch("photonMultiplicity", &photonMultiplicity_,"photonMultiplicity/I");
  bmmgTree_->Branch("photonPt",&photonPt_,"photonPt/D");
  bmmgTree_->Branch("photonEta",&photonEta_,"photonEta/D");
  bmmgTree_->Branch("photonPhi",&photonPhi_,"photonPhi/D");
  bmmgTree_->Branch("photonTrkIso",&photonTrkIso_,"photonTrkIso/D");
  bmmgTree_->Branch("photonEcalIso",&photonEcalIso_,"photonEcalIso/D");
  bmmgTree_->Branch("photonHcalIso",&photonHcalIso_,"photonHcalIso/D");
  bmmgTree_->Branch("photonCaloIso",&photonCaloIso_,"photonCaloIso/D");
  bmmgTree_->Branch("photonSigmaiEtaiEta",&photonSigmaiEtaiEta_,"photonSigmaiEtaiEta/D");
  bmmgTree_->Branch("photonSigmaiEtaiPhi",&photonSigmaiEtaiPhi_,"photonSigmaiEtaiPhi/D");
  bmmgTree_->Branch("photonSigmaiPhiiPhi",&photonSigmaiPhiiPhi_,"photonSigmaiPhiiPhi/D");
  bmmgTree_->Branch("photonSigmaEtaEta",&photonSigmaEtaEta_,"photonSigmaEtaEta/D");
  bmmgTree_->Branch("photone1x5",&photone1x5_,"photone1x5/D");
  bmmgTree_->Branch("photone2x5",&photone2x5_,"photone2x5/D");
  bmmgTree_->Branch("photone3x3",&photone3x3_,"photone3x3/D");
  bmmgTree_->Branch("photone5x5",&photone5x5_,"photone5x5/D");
  bmmgTree_->Branch("photonhcalDepth1OverEcal",&photonhcalDepth1OverEcal_,"photonhcalDepth1OverEcal/D");
  bmmgTree_->Branch("photonhcalDepth2OverEcal",&photonhcalDepth2OverEcal_,"photonhcalDepth2OverEcal/D");
  bmmgTree_->Branch("photonhcalDepth1OverEcalBc",&photonhcalDepth1OverEcalBc_,"photonhcalDepth1OverEcal/D");
  bmmgTree_->Branch("photonhcalDepth2OverEcalBc",&photonhcalDepth2OverEcalBc_,"photonhcalDepth2OverEcalBc/D");
  bmmgTree_->Branch("electronMultiplicity",&electronMultiplicity_,"electronMultiplicity/D");
  bmmgTree_->Branch("pfCandMultiplicity",&pfCandMultiplicity_,"pfCandMultiplicity/D");
  bmmgTree_->Branch("costheta",&costheta_,"costheta/D");
  bmmgTree_->Branch("phi",&phi_,"phi/D");
  bmmgTree_->Branch("cospsi",&cospsi_,"cospsi/D");
  bmmgTree_->Branch("AngleBsDecayLength",&AngleBsDecayLength_,"AngleBsDecayLength/D");

}

RadiativeRootTree::~RadiativeRootTree()
{}

void RadiativeRootTree::writeFile()
{
  bmmgFile_->Write();
  bmmgFile_->Close();

}
void RadiativeRootTree::resetEntries()
{
	runNumber_          = -9999999;
        eventNumber_        = -9999999;
        lumiSection_        = -9999999;
        PUinteraction_      = -9999999;
	PUTrueinteraction_  = -9999999;
	NVerticesbeforecut_ = -9999999;
	NVerticesaftercut_  = -9999999;
	BSx_                = -9999999;
	BSy_                = -9999999;
	BSz_                = -9999999;
	BSdx_               = -9999999;
       	BSdy_               = -9999999;
	BSdz_               = -9999999;
	BSdxdz_             = -9999999;
	BSdydz_             = -9999999;
	BSsigmaZ_           = -9999999;
	BSdsigmaZ_          = -9999999;
	PVx_                = -9999999;
	PVy_                = -9999999;
	PVz_                = -9999999;
	PVerrx_             = -9999999;
	PVerry_             = -9999999;
	PVerrz_             = -9999999;
	PVndof_             = -9999999;
	PVrho_              = -9999999;
	costheta_           = -9999999;
	phi_                = -9999999;
	cospsi_             = -9999999;
	AngleBsDecayLength_ = -9999999;
	isBS_               = -9999999;
	isPV_               = -9999999;
	dedxTrk_            = -9999999;
	errdedxTrk_         = -9999999;
	numdedxTrk_         = -9999999;
	photonMultiplicity_ = -9999999;
	photonPt_           = -9999999;
        photonEta_          = -9999999;
        photonPhi_          = -9999999;
	photonTrkIso_       = -9999999;
	photonEcalIso_        = -9999999;
	photonHcalIso_        = -9999999;
	photonCaloIso_        = -9999999;
	photonSigmaiEtaiEta_  = -9999999;
	photonSigmaiEtaiPhi_  = -9999999;
	photonSigmaiPhiiPhi_  = -9999999;
	photonSigmaEtaEta_    = -9999999;
	photone1x5_           = -9999999;
	photone2x5_           = -9999999;
	photone3x3_           = -9999999;
	photone5x5_           = -9999999;
	photonhcalDepth1OverEcal_ = -9999999;
	photonhcalDepth2OverEcal_ = -9999999;
	photonhcalDepth1OverEcalBc_ = -9999999;
	photonhcalDepth2OverEcalBc_ = -9999999;
	electronMultiplicity_ = -9999999;
	pfCandMultiplicity_   = -9999999;
	for(int i=0; i<9;++i){PVcovariance_[i] = -9999999;}
}

void RadiativeRootTree::getDeDx(const double f1, const double f2, const int f3)
{
  dedxTrk_ = f1;
  errdedxTrk_ = f2;
  numdedxTrk_ = f3;
}

void RadiativeRootTree::getVtx(const double aa, const double bb, const double cc, const double dd, const double ee, const double ff,
                                 const double gg, const double hh, const double ii)
{
  BSx_ = aa;
  BSy_ = bb;
  BSz_ = cc;
  PVx_ = dd;
  PVy_ = ee;
  PVz_ = ff;
  PVerrx_ = gg;
  PVerry_ = hh;
  PVerrz_ = ii;
}


void RadiativeRootTree::getAngles(const double aa, const double bb, const double cc, const double dd)
{
  costheta_ = aa;
  phi_ = bb;
  cospsi_ = cc;
  AngleBsDecayLength_ = dd;
}

void RadiativeRootTree::fill()
{
  bmmgTree_->Fill();
}
void RadiativeRootTree::readTree(const std::string filename)
{
  bmmgFile_ = new TFile (filename.c_str(), "READ" );
  bmmgTree_ =  (TTree*) bmmgFile_->Get("BMMGTree");
  setBranchAddresses();
}

void RadiativeRootTree::readTree(std::vector<std::string> filenames){
  TChain * myChain = new TChain("BMMGTree");
  for(std::vector<std::string>::iterator it = filenames.begin(); it != filenames.end(); it++){
    myChain->Add( (*it).c_str());
  }
  bmmgTree_ = myChain;
  setBranchAddresses();
}
void RadiativeRootTree::setBranchAddresses(){
  bmmgTree_->SetBranchAddress("runNumber", &runNumber_);
  bmmgTree_->SetBranchAddress("eventNumber", &eventNumber_);
  bmmgTree_->SetBranchAddress("lumiSection", &lumiSection_ );
  bmmgTree_->SetBranchAddress("PUinteraction", &PUinteraction_);
  bmmgTree_->SetBranchAddress("PUTrueinteraction", &PUTrueinteraction_);
  bmmgTree_->SetBranchAddress("NVerticesbeforecut",&NVerticesbeforecut_);
  bmmgTree_->SetBranchAddress("NVerticesaftercut",&NVerticesaftercut_);
  bmmgTree_->SetBranchAddress("BSx" , &BSx_ );
  bmmgTree_->SetBranchAddress("BSy" , &BSy_ );
  bmmgTree_->SetBranchAddress("BSz", &BSz_ );
  bmmgTree_->SetBranchAddress("BSdx", &BSdx_ );
  bmmgTree_->SetBranchAddress("BSdy"  , &BSdy_ );
  bmmgTree_->SetBranchAddress("BSdz" , &BSdz_ );
  bmmgTree_->SetBranchAddress("BSdxdz",&BSdxdz_);
  bmmgTree_->SetBranchAddress("BSdydz",&BSdydz_);
  bmmgTree_->SetBranchAddress("BSsigmaZ",&BSsigmaZ_);
  bmmgTree_->SetBranchAddress("BSdsigmaZ",&BSdsigmaZ_);
  bmmgTree_->SetBranchAddress("PVx", &PVx_);
  bmmgTree_->SetBranchAddress("PVy" , &PVy_);
  bmmgTree_->SetBranchAddress("PVz", &PVz_);
  bmmgTree_->SetBranchAddress("PVerrx", &PVerrx_);
  bmmgTree_->SetBranchAddress("PVerry", &PVerry_);
  bmmgTree_->SetBranchAddress("PVerrz", &PVerrz_);
  bmmgTree_->SetBranchAddress("PVrho", &PVrho_);
  bmmgTree_->SetBranchAddress("PVndof", &PVndof_);
  bmmgTree_->SetBranchAddress("PVcovariance", &PVcovariance_);
  bmmgTree_->SetBranchAddress("isPV", &isPV_);
  bmmgTree_->SetBranchAddress("isBS" , &isBS_);
  bmmgTree_->SetBranchAddress("dedxTrk", &dedxTrk_);
  bmmgTree_->SetBranchAddress("errdedxTrk", &errdedxTrk_);
  bmmgTree_->SetBranchAddress("numdedxTrk", &numdedxTrk_);
  bmmgTree_->SetBranchAddress("costheta", &costheta_);
  bmmgTree_->SetBranchAddress("photonMultiplicity", &photonMultiplicity_ );
  bmmgTree_->SetBranchAddress("photonPt", &photonPt_ );
  bmmgTree_->SetBranchAddress("photonEta", &photonEta_ );
  bmmgTree_->SetBranchAddress("photonPhi", &photonPhi_ );
  bmmgTree_->SetBranchAddress("photonTrkIso", &photonTrkIso_);
  bmmgTree_->SetBranchAddress("photonEcalIso", &photonEcalIso_);
  bmmgTree_->SetBranchAddress("photonHcalIso", &photonHcalIso_);
  bmmgTree_->SetBranchAddress("photonCaloIso", &photonCaloIso_);
  bmmgTree_->SetBranchAddress("photonSigmaiEtaiEta", &photonSigmaiEtaiEta_);
  bmmgTree_->SetBranchAddress("photonSigmaiEtaiPhi", &photonSigmaiEtaiPhi_);
  bmmgTree_->SetBranchAddress("photonSigmaiPhiiPhi", &photonSigmaiPhiiPhi_);
  bmmgTree_->SetBranchAddress("photonSigmaEtaEta", &photonSigmaEtaEta_);
  bmmgTree_->SetBranchAddress("photone1x5", &photone1x5_);
  bmmgTree_->SetBranchAddress("photone2x5", &photone2x5_);
  bmmgTree_->SetBranchAddress("photone3x3", &photone3x3_);
  bmmgTree_->SetBranchAddress("photone5x5", &photone5x5_);
  bmmgTree_->SetBranchAddress("photonhcalDepth1OverEcal", &photonhcalDepth1OverEcal_);
  bmmgTree_->SetBranchAddress("photonhcalDepth2OverEcal", &photonhcalDepth2OverEcal_);
  bmmgTree_->SetBranchAddress("photonhcalDepth1OverEcalBc", &photonhcalDepth1OverEcalBc_);
  bmmgTree_->SetBranchAddress("photonhcalDepth2OverEcal", &photonhcalDepth2OverEcalBc_);
  bmmgTree_->SetBranchAddress("electronMultiplicity", &electronMultiplicity_);
  bmmgTree_->SetBranchAddress("pfCandMultiplicity", &pfCandMultiplicity_);
}

