#list of sampes HTo2LongLivedTo4mu_MH-125_MFF-25_CTau-1500mm.list, Nu_E10.list

#cmsDriver command for neutrino
#without gen
#cmsDriver.py Nu_E10-pythia8-gun_l1Ntuple -s RAW2DIGI --python_filename=Nu_E10_mc.py -n 50 --era=Run3 --mc --conditions=112X_mcRun3_2021_realistic_v15 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulMCFromRAW --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleEMU --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2018_v1_3 --filein=/store/mc/Run3Winter20DRPremixMiniAOD/Nu_E10-pythia8-gun/GEN-SIM-RAW/SNB_110X_mcRun3_2021_realistic_v6-v1/10000/51FBC4DC-D5C7-824F-9AFA-3025F04F96FA.root --no_exec
#with gen
#cmsDriver.py Nu_E10-pythia8-gun_l1Ntuple -s RAW2DIGI --python_filename=Nu_E10_config/Nu_E10_mc_$2.py -n 114000 --era=Run3 --mc --conditions=112X_mcRun3_2021_realistic_v15 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulMCFromRAW --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMUGEN_MC --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2018_v1_3 --filein=$1 --no_exec

#cmsDriver command for signal
#cmsDriver.py HTo2LongLivedTo4mu_MH-125_MFF-25_CTau-1500mm_l1Ntuple -s RAW2DIGI --python_filename=HTo2LongLivedTo4mu_MH-125_MFF-25_CTau-1500mm_mc.py -n 50 --no_output --era=Run3 --mc --conditions=112X_mcRun3_2021_realistic_v15 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulMCFromRAW --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleEMU --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2018_v1_3 --filein=/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4mu_MH-125_MFF-25_CTau-1500mm_TuneCP5_14TeV_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v2/240000/C529FF36-B03F-E441-9E51-1B93A740240E.root --no_exec
cmsDriver.py HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_l1Ntuple -s RAW2DIGI --python_filename=signal_config/HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_mc_$2.py -n 3010 --no_output --era=Run3 --mc --conditions=112X_mcRun3_2021_realistic_v15 --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulMCFromRAW --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMUGEN_MC --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2018_v1_3 --filein=$1 --no_exec

#apply fix to neutrino
#cat fix_recipe.py >> Nu_E10_config/Nu_E10_mc_$2.py

#apply fix to signal
cat fix_recipe.py >> signal_config/HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_mc_$2.py 
