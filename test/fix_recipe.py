#fixme based on: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideL1TStage2Instructions#Workflows
import os                                                                                                                                                           
base = os.environ["CMSSW_BASE"]                                                                                                                                     
process.GlobalTag.toGet = cms.VPSet(                                                                                                                                
        cms.PSet(record = cms.string("GEMeMapRcd"),                                                                                                                 
                       tag = cms.string("GEMeMapDummy"),                                                                                                            
                       connect = cms.string("sqlite_file:" + base + "/src/L1Trigger/Configuration/test/GEMeMapDummy.db")                                            
                )                                                                                                                                                   
)                                                                                                                                                                   
process.muonGEMDigis.useDBEMap = True       

# fixme based on: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HcalPileupMitigation#PFA1_Filter
process.load("SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff")

process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHBQIE11 = 1
process.simHcalTriggerPrimitiveDigis.numberOfFilterPresamplesHEQIE11 = 1
process.simHcalTriggerPrimitiveDigis.weightsQIE11 = (
    "ieta1",  [0.0, 1.0],
    "ieta2",  [0.0, 1.0],
    "ieta3",  [0.0, 1.0],
    "ieta4",  [0.0, 1.0],
    "ieta5",  [0.0, 1.0],
    "ieta6",  [0.0, 1.0],
    "ieta7",  [0.0, 1.0],
    "ieta8",  [0.0, 1.0],
    "ieta9",  [0.0, 1.0],
    "ieta10", [0.0, 1.0],
    "ieta11", [0.0, 1.0],
    "ieta12", [0.0, 1.0],
    "ieta13", [0.0, 1.0],
    "ieta14", [0.0, 1.0],
    "ieta15", [0.0, 1.0],
    "ieta16", [0.0, 1.0],
    "ieta17", [0.0, 1.0],
    "ieta18", [0.0, 1.0],
    "ieta19", [0.0, 1.0],
    "ieta20", [0.0, 1.0],
    "ieta21", [0.0, 1.0],
    "ieta22", [0.0, 1.0],
    "ieta23", [0.0, 1.0],
    "ieta24", [0.0, 1.0],
    "ieta25", [0.0, 1.0],
    "ieta26", [0.0, 1.0],
    "ieta27", [0.0, 1.0],
    "ieta28", [0.0, 1.0]
)

process.HcalTPGCoderULUT.contain1TSHB = True
process.HcalTPGCoderULUT.contain1TSHE = True

# Pick one of the pairs of lines below based on the intended scenario for running
process.HcalTPGCoderULUT.containPhaseNSHB = 3.0 # For Run3 MC
process.HcalTPGCoderULUT.containPhaseNSHE = 3.0 # For Run3 MC

#process.HcalTPGCoderULUT.containPhaseNSHB = 0.0 # For Run2 2018 Data
#process.HcalTPGCoderULUT.containPhaseNSHE = 0.0 # For Run2 2018 Data


