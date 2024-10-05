## Ntuple Production for the  BsToMuMuGamma  analysis 
## Contributors
- ```Muhammad Alibordi```
(with update with full names soon)
### Generic steps : Clone & Run : 
```bash
cmsrel CMSSW_XXX
cd CMSSW_XXX/src
cmsenv
git clone -b amdevel_v1 https://github.com/mumuhamm/BsToMuMuGammaAnalysis.git 
scram b -j 8
cd BsToMuMuGammaAnalysis/RadiativeAnalysis/test
cmsRun makeTree_BsMuMuGamma_MC.py
```
### Description :

- ```RadiativeAnalysis.cc``` : Main code of doing calculations and fill the root tree
- ```RadiativeRootTree.cc``` : Root tree definition  
- ```KinematicConstrainedFit.cc``` : Kinematic fit to obtain the constrained vertex, as for the low momentum gamma we have to an significat amount of work ( Alibordi)   
- ```BsToMuMuGammaAnalysis/RadiativeAnalysis/test/Dynamic_Crab.py``` : We need a dynamic multi-crab submission script, in order to produce ntuples in one go
- ```/RadiativeAnalysis/test/fragmentRequest.sh``` : genFragment download script 

 
### For AOD MC production : 
-Work in Progress, please have a look in another repository : in a branch devel_GPExtrapolationPhaseI
```bash
for now : 
cmsrel CMSSW_14_0_15_patch1
cd CMSSW_14_0_15_patch1/src
git cms-init 
mkdir -p Configuration/GenProduction/
git clone git@github.com:cms-sw/genproductions.git Configuration/GenProduction
mv  Configuration/GenProduction/genfragments Configuration/GenProduction/python
rm -rf  Configuration/GenProduction/python/ThirteenTeV/DisappTrksAMSB/
rm -rf  Configuration/GenProduction/python/ThirteenTeV/DelayedJets/
rm -rf  Configuration/GenProduction/python/ThirteenTeV/DMSIMP_Extensions
rm -f   Configuration/GenProduction/python/EightTeV/Exotica_HSCP_SIM_cfi.py
scram b -j 4

git clone -b devel_GPExtrapolationPhaseI https://github.com/mumuhamm/PrivateMCProduction.git
git clone https://github.com/mumuhamm/GeneratorInterface-EvtGenInterface.git  : for DEC,dec,pdl 
```

-```GenFragments/privateCustomizations.py``` : Private customization in order to keep the sample size resonable
-```./submitJobs_BsToMuMuGamma.py```         : Produces AOD in one go. Though, an step wise generation is described elsewhere : ```https://twiki.cern.ch/twiki/bin/viewauth/CMS/GeneratorMain``` 
### Satements :  
Comming soon

