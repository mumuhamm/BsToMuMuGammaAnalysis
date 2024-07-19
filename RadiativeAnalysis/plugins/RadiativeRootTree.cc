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
  bmmgTree_->Branch("photonMultiplicity", &photonMultiplicity_, "photonMultiplicity/I");
  bmmgTree_->Branch("photonEta", &photonEta_, "photonEta/D");
  bmmgTree_->Branch("photonPhi", &photonPhi_, "photonPhi/D");
  bmmgTree_->Branch("phi", &phi_,"phi/D");
  bmmgTree_->Branch("cospsi",&cospsi_,"cospsi/D");
  bmmgTree_->Branch("AngleBsDecayLength",&AngleBsDecayLength_,"AngleBsDecayLength/D");
  bmmgTree_->Branch("BSx", &BSx_, "BSx/D");
  bmmgTree_->Branch("BSy" , &BSy_,  "BSy/D");
  bmmgTree_->Branch("BSz" , &BSz_, "BSz/D");
  bmmgTree_->Branch("BSdx" , &BSdx_,  "BSdx/D");
  bmmgTree_->Branch("BSdy"  , &BSdy_, "BSdy/D");
  bmmgTree_->Branch("BSdz" , &BSdz_, "BSdz/D");
  bmmgTree_->Branch("PVx" , &PVx_, "PVx/D");
  bmmgTree_->Branch("PVy", &PVy_,"PVy/D");
  bmmgTree_->Branch("PVz", &PVz_, "PVz/D");
  bmmgTree_->Branch("PVerrx", &PVerrx_,  "PVerrx/D");
  bmmgTree_->Branch("PVerry", &PVerry_, "PVerry/D");
  bmmgTree_->Branch("PVerrz", &PVerrz_, "PVerrz/D");
  bmmgTree_->Branch("dedxTrk", &dedxTrk_, "dedxTrk/D");
  bmmgTree_->Branch("errdedxTrk", &errdedxTrk_, "errdedxTrk/D");
  bmmgTree_->Branch("numdedxTrk", &numdedxTrk_, "numdedxTrk/I");
  bmmgTree_->Branch("isPV",&isPV_,"isPV/I");
  bmmgTree_->Branch("isBS",&isBS_,"isBS/I");
  bmmgTree_->Branch("costheta",&costheta_,"costheta/D");

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
	photonMultiplicity_ = -9999999;
	photonEta_  = -9999999;
	photonPhi_  = -9999999;
	BSx_ = -9999999;
  BSy_ = -9999999;
  BSz_ = -9999999;
  BSdx_ = -9999999;
  BSdy_ = -9999999;
  BSdz_ = -9999999;
  PVx_ = -9999999;
  PVy_ = -9999999;
  PVz_ = -9999999;
  PVerrx_ = -9999999;
  PVerry_ = -9999999;
  PVerrz_ = -9999999;
  costheta_ = -9999999;
  phi_ = -9999999;
  cospsi_ = -9999999;
  AngleBsDecayLength_ = -9999999;
  isBS_ = -9999999;
  isPV_ = -9999999;
  dedxTrk_ = -9999999;
  errdedxTrk_ = -9999999;
  numdedxTrk_ = -9999999;
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

  bmmgTree_->SetBranchAddress("photonMultiplicity", &photonMultiplicity_ );
  bmmgTree_->SetBranchAddress("photonEta", &photonEta_ );
  bmmgTree_->SetBranchAddress("photonEta", &photonEta_ );
  bmmgTree_->SetBranchAddress("BSx" , &BSx_ );
  bmmgTree_->SetBranchAddress("BSy" , &BSy_ );
  bmmgTree_->SetBranchAddress("BSz", &BSz_ );
  bmmgTree_->SetBranchAddress("BSdx", &BSdx_ );
  bmmgTree_->SetBranchAddress("BSdy"  , &BSdy_ );
  bmmgTree_->SetBranchAddress("BSdz" , &BSdz_ );
  bmmgTree_->SetBranchAddress("PVx"				  , &PVx_  );
  bmmgTree_->SetBranchAddress("PVy"				  , &PVy_  );
  bmmgTree_->SetBranchAddress("PVz"				  , &PVz_  );
  bmmgTree_->SetBranchAddress(  "PVerrx"			  , &PVerrx_  );
  bmmgTree_->SetBranchAddress(  "PVerry"			  , &PVerry_  );
  bmmgTree_->SetBranchAddress(  "PVerrz"			  , &PVerrz_  );
  bmmgTree_->SetBranchAddress(  "isPV"				  , &isPV_  );
  bmmgTree_->SetBranchAddress(  "isBS"				  , &isBS_  );
  bmmgTree_->SetBranchAddress(  "dedxTrk"			  , &dedxTrk_  );
  bmmgTree_->SetBranchAddress(  "errdedxTrk"			  , &errdedxTrk_  );
  bmmgTree_->SetBranchAddress(  "numdedxTrk"			  , &numdedxTrk_  );
  bmmgTree_->SetBranchAddress(  "costheta"			  , &costheta_  );
}

