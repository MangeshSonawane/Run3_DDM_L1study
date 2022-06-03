#! /usr/bin/env python

## **************************************************************** ##
##  Framework for rate estimation for Single and Double muon seeds  ##
##	Authors : Mangesh Sonawane (mangesh.sonawane@cern.ch)	    ##
##		  Alberto Escalante del Valle (aescalante@cern.ch)  ##
## **************************************************************** ##

import os

import ROOT as R
from math import * 
from array import array

R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

PRT_EVT  = 10000   ## Print every Nth event
MAX_EVT  = 10000 ## Number of events to process
VERBOSE  = False  ## Verbose print-out

style1 = R.TStyle("style1", "For Histograms")
style1.SetLineWidth(2)
style1.SetOptStat(0)

R.gROOT.SetStyle("style1")

def getPhi(globalPhiHw) :
    phi = globalPhiHw/287.5*pi if globalPhiHw<287.5 else \
                 (globalPhiHw-575.)/287.5*pi

    return phi 

def main():

    print '\nInside DisplacedMuons\n'

    evtclassid = 1				## A switch for choosing samples, refers to the index of evtclass and nevt

    evtclass = ["NuGun_rate", "private_NuGun_rate", "signal_3000_rate", "Ephemeral_zero_Bias_rate"]
    nevt = [2500, 88, 14, 47]			## Number of root files in the folder containing the sample ntuples

    inputdir = ['/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/bundocka/condor/reHcalTP_Nu_11_2_105p20p1_1623921599/', 
		'/eos/user/s/sonawane/temp/L1Ntuples/NuGun_tuples/ntuples_01_07_21/Nu_E10-pythia8-gun_01_07_21_1625166942/',
		'/eos/user/s/sonawane/temp/L1Ntuples/signal_tuples/ntuples_01_07_21/HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_11_2_X_1623847924/',
		'/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/bundocka/EphemeralZeroBias' ]

    workdir = '/afs/cern.ch/user/s/sonawane/L1T/L1studies/L1_scripts_Alberto/L1RunIII/macros/'		#working directory
    in_file_names = []
    if evtclassid == 3:											# Ephemeral zero Bias data ntuples are in 8 folders
        for s in 	['1/zbD1/210622_230922/0000/',
			 '2/zbD2/210628_215336/0000/',
			 '3/zbD3/210628_215518/0000/',
 			 '4/zbD4/210707_182431/0000/',
			 '5/zbD5/210707_193231/0000/',
			 '6/zbD6/210707_195834/0000/',
			 '7/zbD7/210707_200137/0000/',
			 '8/zbD8/210707_200259/0000/']:
                ntupledir = inputdir[evtclassid]+s
                for i in range(nevt[evtclassid]) :
                        path = ntupledir+"L1Ntuple_"+str(i)+".root"
                        if not os.path.exists(path): continue
                        in_file_names.append(path)

    else :
	    for i in range(nevt[evtclassid]):
       	 	path = inputdir[evtclassid]+str(i)+".root"
        	if not os.path.exists(path): continue
		in_file_names.append(path)
    
    scale = [2544*11246., 2544*11246., 1., 2544*11246.]				## scaling number of events to obtain pure rate

    if not os.path.exists(workdir+'plots'): os.makedirs(workdir+'plots')

    out_file_str = 'DisplacedMuons_'+evtclass[evtclassid]+"_added_rate_ER2"		## Output file name
    if MAX_EVT >= 1000000 : out_file_str += ('_%dM' % (MAX_EVT / 1000000))
    else : out_file_str += ('_%dk' % (MAX_EVT / 1000))
    out_file = R.TFile(workdir+'plots/'+out_file_str+'.root','recreate')

    chains = {}
    chains['Evt'] = []  ## Event info
   # chains['Emu'] = []  ## Emulated Kalman BMTF
    chains['EmuE'] = []  ## Emulated EMTF
    chains['uGT'] = []  ## Global Trigger
    chains['Evt'].append( R.TChain('l1EventTree/L1EventTree') )
    chains['EmuE'].append( R.TChain('l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree') )
    chains['uGT'].append( R.TChain('l1UpgradeEmuTree/L1UpgradeTree') )

    for i in range(len(in_file_names)):
        print 'Adding file %s' % in_file_names[i]
        chains['Evt'][0].Add( in_file_names[i] )
        chains['EmuE'][0].Add( in_file_names[i] )
        chains['uGT'][0].Add( in_file_names[i] )

    pt_bins = [40, 0, 40]
    dR_bins = [25, 0, 5]
    abs_eta_bins = [0, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5]

# Rate histograms
#    h_dimu_rate_ptVtx 	=	R.TH1F('h_dimu_rate_ptVtx', 	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]',	150, 0, 150)
#    h_dimu_rate_ptDisp 	=	R.TH1F('h_dimu_rate_ptDisp', 	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]', 	150, 0, 150)
#    h_dimu_rate_ptOr 	=	R.TH1F('h_dimu_rate_ptOr', 	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]', 	150, 0, 150)

###############

    ########################
    ### Seed event counter
    ########################

    L1_SingleMu7 	= 0
    L1_SingleMu22 	= 0
    L1_SingleMu25 	= 0
    L1_DoubleMu15_7	= 0
    L1_DoubleMu18_er2p1	= 0

    L1_SingleMu18er1p5			= 0
    L1_DoubleMu0_SQ			= 0
    L1_DoubleMu_15_5_SQ			= 0
    L1_DoubleMu0er1p5_SQ		= 0
    L1_DoubleMu0er1p5_SQ_OS		= 0
    L1_DoubleMu0er1p5_SQ_dR_Max1p4	= 0
    L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4	= 0
    L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7 = 0
    L1_DoubleMu4p5er2p0_SQ_OS		= 0
    L1_DoubleMu4p5_SQ_OS_dR_Max1p2 	= 0
    L1_DoubleMu4_SQ_OS			= 0
    L1_DoubleMu0_SQ_OS			= 0

    L1_SingleMu22	= 0
    L1_SingleMu22_BMTF	= 0
    L1_SingleMu22_OMTF	= 0
    L1_SingleMu22_EMTF	= 0

    L1_DoubleMu15_7_BMTF	= 0
    L1_DoubleMu15_7_BMTF_UPT	= 0
    L1_DoubleMu15_7_UPT		= 0
    L1_DoubleMu15_7_UPT_IP1 	= 0

    L1_DoubleMu_15_7_SQ		= 0
    L1_DoubleMu_15_7_UPT_SQ	= 0
    L1_DoubleMu15_7_UPT_DXY1	= 0 
    L1_DoubleMu0_IP1	= 0 

    Baseline_Run_2			= 0
    Extended_Run_2			= 0
    Extended_Run_2A			= 0
    Extended_Run_2B			= 0
    Baseline_Run_2_BMTF_UPT		= 0
    Baseline_Run_2_UPT			= 0
    Baseline_Run_2_UPT_DXY1		= 0
    Extended_Run_2_UPT_DXY1		= 0
    Extended_Run_2_UPT			= 0

######### new seeds ##########


##############################

    iEvt = 0 
    iEvtPU = 0 

    ###################

    MU_QLTY_SNGL = [12, 13, 14, 15]
    MU_QLTY_DBLE = [8, 9, 10, 11, 12, 13, 14, 15]


    ###################

    print '\nEntering loop over chains'
    for iCh in range(len(chains['uGT'])):

        if iEvt >= MAX_EVT: break

        ## Faster tecnhique, inspired by https://github.com/thomreis/l1tMuonTools/blob/master/L1Analysis.py
        Evt_br = R.L1Analysis.L1AnalysisEventDataFormat()
        Unp_br = R.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
        EmuE_br = R.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
        uGT_br = R.L1Analysis.L1AnalysisL1UpgradeDataFormat()
#        Kmt_br = R.L1Analysis.L1AnalysisBMTFOutputDataFormat()

        chains['Evt'][iCh].SetBranchAddress('Event',               R.AddressOf(Evt_br))
        chains['EmuE'][iCh].SetBranchAddress('L1UpgradeEmtfMuon',   R.AddressOf(EmuE_br))
        chains['uGT'][iCh].SetBranchAddress('L1Upgrade',   	   R.AddressOf(uGT_br))
#        chains['Emu'][iCh].SetBranchAddress('L1UpgradeBmtfOutput', R.AddressOf(Kmt_br))


        print '\nEntering loop over events for chain %d' % iCh
        for jEvt in range(chains['uGT'][iCh].GetEntries()):
	    
            if iEvt >= MAX_EVT: break
	    
            iEvt +=1

            if iEvt % PRT_EVT is 0: print '\nEvent # %d (%dth in chain)' % (iEvt, jEvt+1)

            chains['Evt'][iCh].GetEntry(jEvt)
            chains['EmuE'][iCh].GetEntry(jEvt)
            chains['uGT'][iCh].GetEntry(jEvt)

            # ## Use these lines if you don't explicitly define the DataFormat and then do SetBranchAddress above
            # Evt_br = chains['Evt'][iCh].Event
            # Unp_br = chains['Unp'][iCh].L1UpgradeBmtfMuon
            # Emu_br = chains['Emu'][iCh].L1UpgradeBmtfMuon

	    ev_pileup = int(Evt_br.nPV_True)

            if iEvt % PRT_EVT is 0: print '  * Run %d, LS %d, event %d' % (int(Evt_br.run), int(Evt_br.lumi), int(Evt_br.event))

	    if ev_pileup < 48 or ev_pileup > 58 : continue

	    iEvtPU +=1

            nEMTFMu = int(EmuE_br.nTfMuons)
            nuGTMu = int(uGT_br.nMuons)
#            nKmtMu = int(Kmt_br.nTrks)

#            if (VERBOSE and nUnpMu > 0 or nEmuMu > 0):
#                print 'Unpacked = %d, emulated = %d (internal %d) total muons in collection' % (nUnpMu, nEmuMu, nKmtMu)
#                print 'Unpacked = %d, emulated = %d total muons in collection' % (nUnpMu, nEmuMu)
        

	    uGT_mus = []
	    UPT_uGT_mus = []
	    EMTF_premus = []
	    EMTF_mus = []

###### L1 Seed trigger flags ###########

	    _12_L1_SingleMu7_flag 				= False
	    _19_L1_SingleMu22_flag 				= False
	    _23_L1_SingleMu25_flag 				= False

	    _33_L1_SingleMu18er1p5_flag				= False
	    _41_L1_DoubleMu0_SQ_flag				= False
 	    _42_L1_DoubleMu0_SQ_OS_flag				= False
    	    _47_L1_DoubleMu_15_5_SQ_flag			= False
	    _48_L1_DoubleMu15_7_flag 				= False
	    _49_L1_DoubleMu_15_7_SQ_flag			= False
	    _51_L1_DoubleMu18_er2p1_flag 			= False
	    _55_L1_DoubleMu0er1p5_SQ_flag			= False
	    _56_L1_DoubleMu0er1p5_SQ_OS_flag 			= False
	    _57_L1_DoubleMu0er1p5_SQ_dR_Max1p4_flag 		= False
	    _58_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_flag 		= False
 	    _60_L1_DoubleMu4_SQ_OS_flag				= False
	    _63_L1_DoubleMu4p5_SQ_OS_dR_Max1p2_flag		= False
	    _64_L1_DoubleMu4p5er2p0_SQ_OS_flag			= False
	    _65_L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_flag	= False

	    _20_L1_SingleMu22_BMTF_flag				= False
	    _21_L1_SingleMu22_OMTF_flag				= False
	    _22_L1_SingleMu22_EMTF_flag				= False

############ new seeds flags  #################

	    L1_DoubleMu15_7_BMTF_flag 		= False
	    L1_DoubleMu15_7_BMTF_UPT_flag 	= False
	    L1_DoubleMu15_7_UPT_flag		= False
	    L1_DoubleMu15_7_UPT_DXY1_flag	= False
	    L1_DoubleMu_15_7_UPT_SQ_flag	= False
	    L1_DoubleMu0_IP1_flag		= False

#	    L1_DoubleMu_sublead4_flags		= [False]*200
#	    L1_DoubleMu_sublead7_flags		= [False]*200
#	    L1_DoubleMu_sublead11_flags		= [False]*200

########################################

            #################################
            ###  Emulated (Global) muons  ###
            #################################
            for i in range(nuGTMu):
                BX      = int(uGT_br.muonBx[i])
		ptDisp  = float(uGT_br.muonEtUnconstrained[i])
		
#                if VERBOSE: print 'Emulated Global muon %d BX = %d, qual = %d, ptVtx = %.1f, ptDispl = %.1f, eta = %.2f' % (i, BX, qual, ptVtx, ptDispl, eta)
                
                if (BX  !=  0): continue

		uGT_mus.append(i)

	    for i in uGT_mus :
	
                qual 	= int(uGT_br.muonQual[i])
		ptVtx	= float(uGT_br.muonEt[i])
		eta 	= float(uGT_br.muonEtaAtVtx[i])
		phi1 	= float(uGT_br.muonPhiAtVtx[i]) 
		ptDisp	= float(uGT_br.muonEtUnconstrained[i])
		dxy 	= float(uGT_br.muonDxy[i])
			
		vec1 = R.TLorentzVector()
		vec1.SetPtEtaPhiM(ptVtx, eta, phi1, 1.05e-3)

		if ptVtx>=7. and qual in MU_QLTY_SNGL	: _12_L1_SingleMu7_flag = True
		if ptVtx>=22. and qual in MU_QLTY_SNGL	: _19_L1_SingleMu22_flag = True
		if ptVtx>=25. and qual in MU_QLTY_SNGL	: _23_L1_SingleMu25_flag = True

		if ptVtx>=22. and qual in MU_QLTY_SNGL :
			if abs(eta) <= 0.8				: _20_L1_SingleMu22_BMTF_flag 	= True
			if abs(eta) > 0.8 and abs(eta) <= 1.245 	: _21_L1_SingleMu22_OMTF_flag	= True
			if abs(eta) > 1.245 and abs(eta) < 2.450  	: _22_L1_SingleMu22_EMTF_flag	= True

		if ptVtx >= 18. and qual in MU_QLTY_SNGL and abs(eta) <= 1.5062 : _33_L1_SingleMu18er1p5_flag = True

		## Double muon loop	
		for j in uGT_mus :
			if j <= i : continue
			
			qual2    	= int(uGT_br.muonQual[j])
        	        ptVtx2   	= float(uGT_br.muonEt[j])
	                ptDisp2  	= float(uGT_br.muonEtUnconstrained[j])
			eta2	 	= float(uGT_br.muonEtaAtVtx[j])
			phi2 		= float(uGT_br.muonPhiAtVtx[j]) 

			vec2 = R.TLorentzVector()
			vec2.SetPtEtaPhiM(ptVtx2, eta2, phi2, 1.05e-3)

			pt1 = max(ptVtx, ptVtx2)
			pt2 = min(ptVtx, ptVtx2)

			ptdisp1 = max(ptDisp, ptDisp2)
			ptdisp2 = min(ptDisp, ptDisp2)

			ch1 = int(uGT_br.muonChg[i])
			ch2 = int(uGT_br.muonChg[j])

			dxy1 = -1
			dxy2 = -1

			if (ptDisp >= ptDisp2) :
				dxy1 = int(uGT_br.muonDxy[i])
				dxy2 = int(uGT_br.muonDxy[j])

			else :
				dxy1 = int(uGT_br.muonDxy[j])
				dxy2 = int(uGT_br.muonDxy[i])
				
			dR = vec1.DeltaR(vec2)
#			dPhi = abs(phi1-phi2)
#			if dPhi >= 3.14 : dPhi = 2*3.14 - dPhi
#			dR = sqrt((eta-eta2)**2 + dPhi**2)
			M0 = (vec1+vec2).M()

			if (qual in MU_QLTY_DBLE) and (qual2 in MU_QLTY_DBLE):
				if pt1 >= 15. and pt2 >= 7.	: _48_L1_DoubleMu15_7_flag = True
				if (pt1 >= 15. and abs(eta) < 0.8) and (pt2 >= 7. and abs(eta2) < 0.8)	: L1_DoubleMu15_7_BMTF_flag = True
 
				if (ptdisp1 >= 15. and abs(eta) < 0.8) and (ptdisp2 >= 7. and abs(eta2) < 0.8)	: L1_DoubleMu15_7_BMTF_UPT_flag = True
				if (ptdisp1 >= 15. and ptdisp2 >= 7. )						: L1_DoubleMu15_7_UPT_flag = True

				if (ptdisp1 >= 15.  and ptdisp2 >= 7. and dxy1 >=1 )		: L1_DoubleMu15_7_UPT_DXY1_flag = True
				if (dxy1 >= 1 )		: L1_DoubleMu0_IP1_flag = True


			if (qual in MU_QLTY_SNGL) and (qual2 in MU_QLTY_SNGL):

				if pt1 >= 15. and pt2 >= 5. : _47_L1_DoubleMu_15_5_SQ_flag = True

				if ptVtx >= 18. and ptVtx2 >= 18. :
					if abs(eta) <= 2.104 and abs(eta2) <= 2.104	: _51_L1_DoubleMu18_er2p1_flag = True

#				dPhi = vec1.DeltaPhi(vec2)

				if ptVtx >= 4. and ptVtx2 >= 4. and ch1*ch2 < 0. : _60_L1_DoubleMu4_SQ_OS_flag = True

				if ptVtx >= 0. and ptVtx2 >= 0. : _41_L1_DoubleMu0_SQ_flag	= True

				if ptVtx >= 0. and ptVtx2 >= 0. and abs(eta) <= 1.506 and abs(eta2) <= 1.506 : _55_L1_DoubleMu0er1p5_SQ_flag = True

				if ptVtx >= 0. and ptVtx2 >= 0. :
					if abs(eta) <= 1.506 and abs(eta2) <= 1.506 and ch1*ch2 < 0 : _56_L1_DoubleMu0er1p5_SQ_OS_flag = True

				if ptVtx >= 0. and ptVtx2 >= 0. and ch1*ch2 < 0 : _42_L1_DoubleMu0_SQ_OS_flag = True

				if (ptVtx >= 0. and ptVtx2 >= 0.) and (abs(eta) <= 1.506 and abs(eta2) <= 1.506):
					if (ch1*ch2 < 0 and dR <= 1.4)	: _58_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_flag = True

					if dR <= 1.4 : _57_L1_DoubleMu0er1p5_SQ_dR_Max1p4_flag = True

				if (ptVtx >= 4.5 and ptVtx2 >= 4.5) :
					if (ch1*ch2 < 0 and dR <= 1.2) : _63_L1_DoubleMu4p5_SQ_OS_dR_Max1p2_flag = True
					if (abs(eta) <= 2.006 and abs(eta2) <= 2.006) and (ch1*ch2 < 0) : 
						_64_L1_DoubleMu4p5er2p0_SQ_OS_flag = True
						if M0 >= 7. : _65_L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_flag = True

                                if pt1 >= 15. and pt2 >= 7. : _49_L1_DoubleMu_15_7_SQ_flag = True
                                if ptdisp1 >= 15. and ptdisp2 >= 7. : L1_DoubleMu_15_7_UPT_SQ_flag = True

	    if (_12_L1_SingleMu7_flag)					: L1_SingleMu7 				+=1
	    if (_19_L1_SingleMu22_flag)					: L1_SingleMu22 			+=1
	    if (_20_L1_SingleMu22_BMTF_flag)				: L1_SingleMu22_BMTF 			+=1
	    if (_21_L1_SingleMu22_OMTF_flag)				: L1_SingleMu22_OMTF 			+=1
	    if (_22_L1_SingleMu22_EMTF_flag)				: L1_SingleMu22_EMTF 			+=1
	    if (_23_L1_SingleMu25_flag)					: L1_SingleMu25 			+=1
	    if (_33_L1_SingleMu18er1p5_flag)				: L1_SingleMu18er1p5 			+=1
	    if (_41_L1_DoubleMu0_SQ_flag)				: L1_DoubleMu0_SQ 			+=1
	    if (_42_L1_DoubleMu0_SQ_OS_flag)				: L1_DoubleMu0_SQ_OS 			+=1
	    if (_47_L1_DoubleMu_15_5_SQ_flag)				: L1_DoubleMu_15_5_SQ 			+=1
	    if (_48_L1_DoubleMu15_7_flag)				: L1_DoubleMu15_7 			+=1
	    if (_49_L1_DoubleMu_15_7_SQ_flag)				: L1_DoubleMu_15_7_SQ			+=1
	    if (_51_L1_DoubleMu18_er2p1_flag)				: L1_DoubleMu18_er2p1 			+=1
	    if (_55_L1_DoubleMu0er1p5_SQ_flag)				: L1_DoubleMu0er1p5_SQ 			+=1
	    if (_56_L1_DoubleMu0er1p5_SQ_OS_flag)			: L1_DoubleMu0er1p5_SQ_OS 		+=1
	    if (_57_L1_DoubleMu0er1p5_SQ_dR_Max1p4_flag)		: L1_DoubleMu0er1p5_SQ_dR_Max1p4 	+=1
	    if (_58_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_flag)		: L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4 	+=1
	    if (_60_L1_DoubleMu4_SQ_OS_flag)				: L1_DoubleMu4_SQ_OS 			+=1
	    if (_63_L1_DoubleMu4p5_SQ_OS_dR_Max1p2_flag)		: L1_DoubleMu4p5_SQ_OS_dR_Max1p2 	+=1
	    if (_64_L1_DoubleMu4p5er2p0_SQ_OS_flag)			: L1_DoubleMu4p5er2p0_SQ_OS 		+=1
	    if (_65_L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7_flag)		: L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7 	+=1

	    if (L1_DoubleMu15_7_BMTF_flag)				: L1_DoubleMu15_7_BMTF 			+=1
	    if (L1_DoubleMu15_7_BMTF_UPT_flag)				: L1_DoubleMu15_7_BMTF_UPT 		+=1
	    if (L1_DoubleMu15_7_UPT_flag)				: L1_DoubleMu15_7_UPT	 		+=1
	    if (L1_DoubleMu15_7_UPT_DXY1_flag)				: L1_DoubleMu15_7_UPT_IP1		+=1
	    if (L1_DoubleMu_15_7_UPT_SQ_flag)				: L1_DoubleMu_15_7_UPT_SQ		+=1
	    if (L1_DoubleMu15_7_UPT_DXY1_flag)				: L1_DoubleMu15_7_UPT_DXY1		+=1 
	    if (L1_DoubleMu0_IP1_flag)					: L1_DoubleMu0_IP1			+=1

##############################

	    BR2_flag = False
	    ER2_flag = False

	    if (_48_L1_DoubleMu15_7_flag)	: 
		Baseline_Run_2	+=1
		BR2_flag = True
	    if (_48_L1_DoubleMu15_7_flag or _58_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_flag or _63_L1_DoubleMu4p5_SQ_OS_dR_Max1p2_flag) 	: 
		ER2_flag = True
		Extended_Run_2	+=1

	    if (_48_L1_DoubleMu15_7_flag or _58_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_flag) 						: Extended_Run_2A	+=1
	    if (_48_L1_DoubleMu15_7_flag or _63_L1_DoubleMu4p5_SQ_OS_dR_Max1p2_flag) 							: Extended_Run_2B	+=1
	    if (_48_L1_DoubleMu15_7_flag or L1_DoubleMu15_7_BMTF_UPT_flag ) 								: Baseline_Run_2_BMTF_UPT +=1
	    if (_48_L1_DoubleMu15_7_flag or L1_DoubleMu15_7_UPT_flag) 									: Baseline_Run_2_UPT 	+=1
	    if (_48_L1_DoubleMu15_7_flag or L1_DoubleMu15_7_UPT_DXY1_flag)								: Baseline_Run_2_UPT_DXY1 +=1
	    if (ER2_flag or L1_DoubleMu15_7_UPT_DXY1_flag) :
		Extended_Run_2_UPT_DXY1 	+=1
	    if (ER2_flag or L1_DoubleMu15_7_UPT_flag): 
		Extended_Run_2_UPT 		+=1

############################# Block to pick leading and subleading PT for Double Muon seed studies ###################

	    pt_lead 	= -10.
	    pt_sublead 	= -10.

	    leadid	= -1.
	    subleadid	= -1.

	    for f in uGT_mus :

		qual 	= int(uGT_br.muonQual[f])
		pt_temp = float(uGT_br.muonEtUnconstrained[f])

		if qual in MU_QLTY_DBLE :	

			if pt_temp >= pt_lead :	
				pt_sublead = pt_lead
				pt_lead = pt_temp
				subleadid = leadid
				leadid = f

			elif pt_temp >= pt_sublead :
				pt_sublead = pt_temp
				subleadid  = f

######################################################################################################################

        ## End loop: for jEvt in range(chains['Unp'][iCh].GetEntries()):
    ## End loop: for iCh in range(len(chains['Unp'])):

    ########Finish loop over chains/Events

    print '\nFinished loop over chains'

    out_file.cd()

#########################################################################

    print '\nWrote out file: plots/'+out_file_str+'.root'
    print "Events run over: ", iEvt
    print "Events run over in PU range [48, 58]: ", iEvtPU

    print "%3d. %-50s : %-7d \t rate : %.3f"%(12, "L1_SingleMu7", L1_SingleMu7, L1_SingleMu7*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(19, "L1_SingleMu22", L1_SingleMu22, L1_SingleMu22*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(20, "L1_SingleMu22_BMTF", L1_SingleMu22_BMTF, L1_SingleMu22_BMTF*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(21, "L1_SingleMu22_OMTF", L1_SingleMu22_OMTF, L1_SingleMu22_OMTF*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(22, "L1_SingleMu22_EMTF", L1_SingleMu22_EMTF, L1_SingleMu22_EMTF*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(23, "L1_SingleMu25", L1_SingleMu25, L1_SingleMu25*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(33, "L1_SingleMu18er1p5", L1_SingleMu18er1p5, L1_SingleMu18er1p5*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(41, "L1_DoubleMu0_SQ", L1_DoubleMu0_SQ, L1_DoubleMu0_SQ*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(42, "L1_DoubleMu0_SQ_OS", L1_DoubleMu0_SQ_OS, L1_DoubleMu0_SQ_OS*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(47, "L1_DoubleMu_15_5_SQ", L1_DoubleMu_15_5_SQ, L1_DoubleMu_15_5_SQ*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(48, "L1_DoubleMu15_7", L1_DoubleMu15_7, L1_DoubleMu15_7*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(49, "L1_DoubleMu15_7_SQ", L1_DoubleMu_15_7_SQ, L1_DoubleMu_15_7_SQ*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(51, "L1_DoubleMu18_er2p1", L1_DoubleMu18_er2p1, L1_DoubleMu18_er2p1*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(55, "L1_DoubleMu0er1p5_SQ", L1_DoubleMu0er1p5_SQ, L1_DoubleMu0er1p5_SQ*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(56, "L1_DoubleMu0er1p5_SQ_OS", L1_DoubleMu0er1p5_SQ_OS, L1_DoubleMu0er1p5_SQ_OS*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(57, "L1_DoubleMu0er1p5_SQ_dR_Max1p4", L1_DoubleMu0er1p5_SQ_dR_Max1p4, L1_DoubleMu0er1p5_SQ_dR_Max1p4*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(58, "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4, L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(60, "L1_DoubleMu4_SQ_OS", L1_DoubleMu4_SQ_OS, L1_DoubleMu4_SQ_OS*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(63, "L1_DoubleMu4p5_SQ_OS_dR_Max1p2", L1_DoubleMu4p5_SQ_OS_dR_Max1p2, L1_DoubleMu4p5_SQ_OS_dR_Max1p2*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(64, "L1_DoubleMu4p5er2p0_SQ_OS", L1_DoubleMu4p5er2p0_SQ_OS, L1_DoubleMu4p5er2p0_SQ_OS*scale[evtclassid]/iEvtPU)
    print "%3d. %-50s : %-7d \t rate : %.3f"%(65, "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7, L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7*scale[evtclassid]/iEvtPU)
#    print "%3d. %-50s : %-7d \t rate : %.3f"%(, , , *scale[evtclassid]/iEvtPU)

    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu15_7_BMTF ", L1_DoubleMu15_7_BMTF, L1_DoubleMu15_7_BMTF*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu15_7_BMTF_UPT ", L1_DoubleMu15_7_BMTF_UPT, L1_DoubleMu15_7_BMTF_UPT*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu15_7_UPT ", L1_DoubleMu15_7_UPT, L1_DoubleMu15_7_UPT*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu15_7_UPT_SQ ", L1_DoubleMu_15_7_UPT_SQ, L1_DoubleMu_15_7_UPT_SQ*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu15_7_UPT_IP1 ", L1_DoubleMu15_7_UPT_DXY1, L1_DoubleMu15_7_UPT_DXY1*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu0_IP1 ", L1_DoubleMu0_IP1, L1_DoubleMu0_IP1*scale[evtclassid]/iEvtPU)

    print " ######################################### Efficiencies ############################################"
    print "%-50s : %-7d \t rate : %.3f"%("Baseline Run 2 ", Baseline_Run_2, Baseline_Run_2*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("Extended Run 2 ", Extended_Run_2, Extended_Run_2*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("Extended Run 2A ", Extended_Run_2A, Extended_Run_2A*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("Extended Run 2B ", Extended_Run_2B, Extended_Run_2B*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("Baseline Run 2 + Unconstrained BMTF ", Baseline_Run_2_BMTF_UPT, Baseline_Run_2_BMTF_UPT*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("Baseline Run 2 + Unconstrained ", Baseline_Run_2_UPT, Baseline_Run_2_UPT*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("Baseline Run 2 + Unconstrained + DXY1 ", Baseline_Run_2_UPT_DXY1, Baseline_Run_2_UPT_DXY1*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("Extended Run 2 + Unconstrained + DXY1 ", Extended_Run_2_UPT_DXY1, Extended_Run_2_UPT_DXY1*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("Extended Run 2 + Unconstrained ", Extended_Run_2_UPT, Extended_Run_2_UPT*scale[evtclassid]/iEvtPU)

############################### Scale rate to kHz #############################

    Scale = scale[evtclassid]/(1e3*iEvtPU)

#########################################################################

    out_file.Close()
    del chains

if __name__ == '__main__':
    main()
