#! /usr/bin/env python

## **************************************************************** ##
##  Look at properties of displaced muons from Kalman algo in BMTF  ##
## **************************************************************** ##

import os

import ROOT as R
from math import * 
from array import array

R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

PRT_EVT  = 10000   ## Print every Nth event
#MAX_EVT  = 1000000 ## Number of events to process
MAX_EVT  = 40000 ## Number of events to process
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
    evtclass = ["NuGun_rate", "private_NuGun_rate", "signal_3000_rate", "Ephemeral_zero_Bias_rate"]
    nevt = [2500, 88, 14, 47]
    evtclassid = 2
    inputdir = ['/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/bundocka/condor/reHcalTP_Nu_11_2_105p20p1_1623921599/', 
		'/eos/user/s/sonawane/temp/L1Ntuples/NuGun_tuples/ntuples_01_07_21/Nu_E10-pythia8-gun_01_07_21_1625166942/',
		'/eos/user/s/sonawane/temp/L1Ntuples/signal_tuples/ntuples_01_07_21/HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_11_2_X_1623847924/',
		'/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/bundocka/EphemeralZeroBias' ]

    workdir = '/afs/cern.ch/user/s/sonawane/L1T/L1studies/L1_scripts_Alberto/L1RunIII/macros/'
    in_file_names = []
    if evtclassid == 3:
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
    
    scale = [2544*11246., 2544*11246., 1., 2544*11246.]

    if not os.path.exists(workdir+'plots'): os.makedirs(workdir+'plots')

    out_file_str = 'DisplacedMuons_'+evtclass[evtclassid]+"_added_rate_ER2"
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

#    h_dimu_rate_ptAnd		= 	R.TH1F('h_dimu_rate_ptAnd',	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]',	150, 0, 150)
#    h_dimu_rate_ptUniqueDisp	= 	R.TH1F('h_dimu_rate_ptUniqueDisp',	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]',	150, 0, 150)
#    h_dimu_rate_ptUniqueVtx	= 	R.TH1F('h_dimu_rate_ptUniqueVtx',	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]',	150, 0, 150)

#    h_kBMTF_singlemu_rate_ptVtx		= 	R.TH1F('h_kBMTF_singlemu_rate_ptVtx',	'Dimuon Pt rate ; L1 Mu Pt; Rate [Hz]',	150, 0, 150)
#    h_kBMTF_singlemu_rate_ptDisp	= 	R.TH1F('h_kBMTF_singlemu_rate_ptDisp',	'Dimuon Pt rate ; L1 Mu Pt; Rate [Hz]',	150, 0, 150)
#    h_kBMTF_singlemu_rate_ptOr		= 	R.TH1F('h_kBMTF_singlemu_rate_ptOr',	'Dimuon Pt rate ; L1 Mu Pt; Rate [Hz]',	150, 0, 150)
#    h_kBMTF_singlemu_rate_ptDisp_Lxy2	= 	R.TH1F('h_kBMTF_singlemu_rate_ptDisp_Lxy2',	'Dimuon Pt rate ; L1 Mu Pt; Rate [Hz]',	150, 0, 150)

#    h_EMTF_rate_vs_ptthreshold_L1dxy0 	= R.TH1F('h_EMTF_rate_vs_ptthreshold_L1dxy0', 'Rate vs L1 Pt Threshold ; L1 pT threshold [GeV]; Rate [kHz]', 300, 0, 300)
#    h_EMTF_rate_vs_ptthreshold_L1dxy1 	= R.TH1F('h_EMTF_rate_vs_ptthreshold_L1dxy1', 'Rate vs L1 Pt Threshold ; L1 pT threshold [GeV]; Rate [kHz]', 300, 0, 300)
#    h_EMTF_rate_vs_ptthreshold_L1dxy2 	= R.TH1F('h_EMTF_rate_vs_ptthreshold_L1dxy2', 'Rate vs L1 Pt Threshold ; L1 pT threshold [GeV]; Rate [kHz]', 300, 0, 300)
#    h_EMTF_rate_vs_ptthreshold_L1dxy3	= R.TH1F('h_EMTF_rate_vs_ptthreshold_L1dxy3', 'Rate vs L1 Pt Threshold ; L1 pT threshold [GeV]; Rate [kHz]', 300, 0, 300)

#    h_rate_vs_l1pt2threshold4		= R.TH1F('h_rate_vs_l1pt2threshold4', 'Rate vs L1 Pt Threshold ; L1 Lead pT threshold [GeV]; Rate [kHz]', 300, 0, 300)
#    h_rate_vs_l1pt2threshold7		= R.TH1F('h_rate_vs_l1pt2threshold7', 'Rate vs L1 Pt Threshold ; L1 Lead pT threshold [GeV]; Rate [kHz]', 300, 0, 300)
#    h_rate_vs_l1pt2threshold11		= R.TH1F('h_rate_vs_l1pt2threshold11', 'Rate vs L1 Pt Threshold ; L1 Lead pT threshold [GeV]; Rate [kHz]', 300, 0, 300)

    h_rate_l1pt1_vs_l1pt2_dxy0		= R.TH2F('h_rate_l1pt1_vs_l1pt2_dxy0', "Lead pt vs sublead pt, Dxy = (0,0), Added rate wrt ER2; Lead muon L1 pt [GeV]; Sublead Muon L1 pt [GeV]", \
					pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])

    h_rate_l1pt1_vs_l1pt2_dxy10		= R.TH2F('h_rate_l1pt1_vs_l1pt2_dxy10', "Lead pt vs sublead pt, Dxy = (1,0), Added rate wrt ER2; Lead muon L1 pt [GeV]; Sublead Muon L1 pt [GeV]", \
					pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])

    h_rate_l1pt1_vs_l1pt2_dxy11		= R.TH2F('h_rate_l1pt1_vs_l1pt2_dxy11', "Lead pt vs sublead pt, Dxy = (1,1), Added rate wrt ER2; Lead muon L1 pt [GeV]; Sublead Muon L1 pt [GeV]", \
					pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])

    h_rate_l1pt1_vs_l1pt2_dxy20		= R.TH2F('h_rate_l1pt1_vs_l1pt2_dxy20', "Lead pt vs sublead pt, Dxy = (2,0), Added rate wrt ER2; Lead muon L1 pt [GeV]; Sublead Muon L1 pt [GeV]", \
					pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])

###############
    h_rate_dR_vs_eta_dxy00	= R.TH2F('h_rate_dR_vs_eta_dxy00', "Added Rate wrt ER2, |#eta| vs #Delta R, Dxy (0,0); Lead muon |#eta|; #Delta R (#mu_{1}, #mu_{2})",\
					len(abs_eta_bins)-1, array('d', abs_eta_bins), dR_bins[0], dR_bins[1], dR_bins[2])
    h_rate_dR_vs_eta_dxy10	= R.TH2F('h_rate_dR_vs_eta_dxy10', "Added Rate wrt ER2, |#eta| vs #Delta R, Dxy (1,0); Lead muon |#eta|; #Delta R (#mu_{1}, #mu_{2})",\
					len(abs_eta_bins)-1, array('d', abs_eta_bins), dR_bins[0], dR_bins[1], dR_bins[2])

    h_dR_vs_leadeta_ER2		= R.TH2F('h_dR_vs_leadeta_ER2', "Differential #Delta R vs |#eta| for ER2; Lead muon |#eta|; #DeltaR (#mu_{1}, #mu_{2}))", \
					len(abs_eta_bins)-1, array('d', abs_eta_bins), dR_bins[0], dR_bins[1], dR_bins[2])
    h_dR_vs_subleadeta_ER2	= R.TH2F('h_dR_vs_subleadeta_ER2', "Differential #Delta R vs |#eta| for ER2; Sublead muon |#eta|; #DeltaR (#mu_{1}, #mu_{2})", \
					len(abs_eta_bins)-1, array('d', abs_eta_bins), dR_bins[0], dR_bins[1], dR_bins[2])
    h_lead_vs_subleadeta_ER2	= R.TH2F('h_lead_vs_subleadeta_ER2', "Differential |#eta_{1}| vs |#eta_{2}| for ER2; Lead muon |#eta|; #Sublead muon |#eta|", \
					len(abs_eta_bins)-1, array('d', abs_eta_bins), len(abs_eta_bins)-1, array('d', abs_eta_bins))
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
    L1_TripleMu_5_3_3			= 0

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
    L1_DoubleMu6_4_UPT_DXY1	= 0 
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

#	    if ev_pileup < 48 or ev_pileup > 58 : continue

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
	    _78_L1_TripleMu_5_3_3_flag 				= False

	    _20_L1_SingleMu22_BMTF_flag				= False
	    _21_L1_SingleMu22_OMTF_flag				= False
	    _22_L1_SingleMu22_EMTF_flag				= False

############ new seeds #################

	    L1_DoubleMu15_7_BMTF_flag 		= False
	    L1_DoubleMu15_7_BMTF_UPT_flag 	= False
	    L1_DoubleMu15_7_UPT_flag		= False
	    L1_DoubleMu15_7_UPT_DXY1_flag	= False
	    L1_DoubleMu6_4_UPT_DXY1_flag	= False
	    L1_DoubleMu_15_7_UPT_SQ_flag	= False
	    L1_DoubleMu0_IP1_flag		= False

#	    L1_DoubleMu_sublead4_flags		= [False]*200
#	    L1_DoubleMu_sublead7_flags		= [False]*200
#	    L1_DoubleMu_sublead11_flags		= [False]*200

########################################

	    #################################
            ###  Emulated (EMTF) muons  ###
            #################################

            for i in range(nEMTFMu):
                BX      = int(EmuE_br.tfMuonBx[i])
                
                if (BX  !=  0): continue

		EMTF_premus.append(i)

#	    if len(EMTF_premus) :	    print "Before cleaning EMTF : ", EMTF_premus

	    EMTF_eta_phi = []

	    for i in EMTF_premus :

		pt1  = float(EmuE_br.tfMuonHwPt[i]-1.)*0.5
		eta1 = float(EmuE_br.tfMuonHwEta[i])*0.010875
		phi1 = getPhi(float(EmuE_br.tfMuonGlobalPhi[i]))

		ep = [eta1, phi1]

		pt_i = [pt1, i]

		ptmax = pt1

		if ep in EMTF_eta_phi : continue

		EMTF_eta_phi.append(ep)

		for j in EMTF_premus :
			if j <= i : continue
			
			pt2 = float(EmuE_br.tfMuonHwPt[j]-1.)*0.5
			eta2 = float(EmuE_br.tfMuonHwEta[j])*0.010875
			phi2 = getPhi(float(EmuE_br.tfMuonGlobalPhi[j]))

			if eta1==eta2 and phi1==phi2 :
				if pt2 > ptmax : 
					ptmax = pt2
					pt_i = [pt2, j]

	    	EMTF_mus.append(pt_i[1])

#	    if len(EMTF_premus) :	    print "After cleaning EMTF : ", EMTF_mus, "\n"



            #################################
            ###  Emulated (Global) muons  ###
            #################################
            for i in range(nuGTMu):
                BX      = int(uGT_br.muonBx[i])
		ptDisp  = float(uGT_br.muonEtUnconstrained[i])
		
#                if VERBOSE: print 'Emulated Global muon %d BX = %d, qual = %d, ptVtx = %.1f, ptDispl = %.1f, eta = %.2f' % (i, BX, qual, ptVtx, ptDispl, eta)
                
                if (BX  !=  0): continue

		uGT_mus.append(i)
#		UPT_uGT_mus.append([ptDisp, i])

#	    UPT_uGT_mus.sort(reverse=True)

	    ## uGT EMTF matching
	
	    EMTF_dR_list = []
	
	    for i in uGT_mus :
		pt1  = float(uGT_br.muonEt[i])
		eta1 = float(uGT_br.muonEtaAtVtx[i])
		phi1 = float(uGT_br.muonPhiAtVtx[i])

		v1 = R.TLorentzVector()
		v1.SetPtEtaPhiM(pt1, eta1, phi1, 0.)

		for j in EMTF_mus :
			pt2  = float(EmuE_br.tfMuonHwPt[j]-1.)*0.5
			eta2 = float(EmuE_br.tfMuonHwEta[j])*0.010875
			phi2 = getPhi(float(EmuE_br.tfMuonGlobalPhi[j]))
			v2 = R.TLorentzVector()
			v2.SetPtEtaPhiM(pt2, eta2, phi2, 0.)

			dR = v1.DeltaR(v2)

			if dR < 0.4 : EMTF_dR_list.append([dR, i, j])

	    if len(EMTF_dR_list) :
		EMTF_dR_list.sort()

                k = True
                i = 0
                while k:
			k = False
                        idx1 = EMTF_dR_list[i][1]
                        idx2 = EMTF_dR_list[i][2]
                        j = i+1
                        while j < len(EMTF_dR_list):
                                idx_1  = EMTF_dR_list[j][1]
                                idx_2  = EMTF_dR_list[j][2]

                                if idx1 == idx_1 or idx2 == idx_2: del EMTF_dR_list[j]
                                else : j +=1
                        i +=1
                        if i < len(EMTF_dR_list):
                                k = True


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

  			for k in uGT_mus :
				if k <= j : continue

				ptVtx3 = float(uGT_br.muonEt[k])
		
				pt_list = [ptVtx, ptVtx2, ptVtx3]
				pt_list.sort(reverse=True)

				if pt_list[0] >= 5. and pt_list[1] >=3. and pt_list[2] >= 3 : _78_L1_TripleMu_5_3_3_flag = True

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

				if (ptdisp1 >= 6.  and ptdisp2 >= 4. and dxy1 >=1 )		: L1_DoubleMu15_7_UPT_DXY1_flag = True
				if (dxy1 >= 1 )		: L1_DoubleMu0_IP1_flag = True

	#			for ptthresh in range(200) :
	#				if ptdisp2 >= 4. and ptdisp1 >= ptthresh : L1_DoubleMu_sublead4_flags[ptthresh] = True
	#				if ptdisp2 >= 7. and ptdisp1 >= ptthresh : L1_DoubleMu_sublead7_flags[ptthresh] = True
	#				if ptdisp2 >= 11. and ptdisp1 >= ptthresh : L1_DoubleMu_sublead11_flags[ptthresh] = True


#				for l in range(150) :
#					ptthresh = float(l)
#					if pt1 >= ptthresh and pt2 >=7. : h_dimu_rate_ptVtx.Fill(ptthresh)
#					if ptdisp1 >= ptthresh and ptdisp2 >=7. : h_dimu_rate_ptDisp.Fill(ptthresh)
#					if (pt1 >= ptthresh and pt2 >=7.) or (ptdisp1 >= ptthresh and ptdisp2 >=7.) : h_dimu_rate_ptOr.Fill(ptthresh)
#					if (pt1 >= ptthresh and pt2 >=7.) and (ptdisp1 >= ptthresh and ptdisp2 >=7.) : h_dimu_rate_ptAnd.Fill(ptthresh)
#					if not (pt1 >= ptthresh and pt2 >=7.) and (ptdisp1 >= ptthresh and ptdisp2 >=7.) : h_dimu_rate_ptUniqueDisp.Fill(ptthresh)
#					if (pt1 >= ptthresh and pt2 >=7.) and not (ptdisp1 >= ptthresh and ptdisp2 >=7.) : h_dimu_rate_ptUniqueVtx.Fill(ptthresh)

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

#	    for el in EMTF_premus :
#		uGT_idx = el[1]		
#		EMTF_idx = el[2]		
	
#		ptVtx 	= float(EmuE_br.tfMuonHwPt[el]-1.)*0.5
#		ptDisp 	= float(EmuE_br.tfMuonHwPtUnconstrained[el]-1.)
#		L1_dxy 	= int(EmuE_br.tfMuonHwDxy[el])

#		qual = int(EmuE_br.tfMuonHwQual[el])

#		if qual not in [11,13,14,15] : continue

	#	for k in range(300):
	#		ptthresh = float(k)

	#		if ptVtx >= ptthresh :
         #                      h_EMTF_rate_vs_ptthreshold_L1dxy0.Fill(ptthresh)

          #                     if L1_dxy >= 1 : h_EMTF_rate_vs_ptthreshold_L1dxy1.Fill(ptthresh)
           #                    if L1_dxy >= 2 : h_EMTF_rate_vs_ptthreshold_L1dxy2.Fill(ptthresh)
            #                   if L1_dxy >= 3 : h_EMTF_rate_vs_ptthreshold_L1dxy3.Fill(ptthresh)

			
			

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
	    if (_78_L1_TripleMu_5_3_3_flag)				: L1_TripleMu_5_3_3 			+=1

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

	    if (_48_L1_DoubleMu15_7_flag or _47_L1_DoubleMu_15_5_SQ_flag)  or _78_L1_TripleMu_5_3_3_flag)	: 
		Baseline_Run_2	+=1
		BR2_flag = True
	    if (BR2_flag or _58_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_flag or _63_L1_DoubleMu4p5_SQ_OS_dR_Max1p2_flag) 	: 
		ER2_flag = True
		Extended_Run_2	+=1

	    if (BR2_flag or _58_L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4_flag) 						: Extended_Run_2A	+=1
	    if (BR2_flag or _63_L1_DoubleMu4p5_SQ_OS_dR_Max1p2_flag) 							: Extended_Run_2B	+=1
	    if (BR2_flag or L1_DoubleMu15_7_BMTF_UPT_flag ) 								: Baseline_Run_2_BMTF_UPT +=1
	    if (BR2_flag or L1_DoubleMu15_7_UPT_flag) 									: Baseline_Run_2_UPT 	+=1
	    if (BR2_flag or L1_DoubleMu15_7_UPT_DXY1_flag)								: Baseline_Run_2_UPT_DXY1 +=1
	    if (ER2_flag or L1_DoubleMu15_7_UPT_DXY1_flag) :
		Extended_Run_2_UPT_DXY1 	+=1
	    if (ER2_flag or L1_DoubleMu15_7_UPT_flag): 
		Extended_Run_2_UPT 		+=1

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

#	    if pt_lead > 40. : pt_lead = 40.
#	    if pt_sublead > 40.. : pt_sublead = 100.

	    if pt_lead >= 0 and pt_sublead >= 0 :

		    eta1 = float(uGT_br.muonEtaAtVtx[leadid])
	            eta2 = float(uGT_br.muonEtaAtVtx[subleadid])
		    phi1 = float(uGT_br.muonPhiAtVtx[leadid])
	            phi2 = float(uGT_br.muonPhiAtVtx[subleadid])
	
        	    dxy1 = float(uGT_br.muonDxy[leadid])

		    vec1 = R.TLorentzVector()
		    vec2 = R.TLorentzVector()

		    vec1.SetPtEtaPhiM(pt_lead, eta1, phi1, 0.)
		    vec2.SetPtEtaPhiM(pt_sublead, eta2, phi2, 0.)
		    dR = vec1.DeltaR(vec2)

		    if ER2_flag :
			h_dR_vs_leadeta_ER2.Fill(abs(eta), dR)
			h_dR_vs_subleadeta_ER2.Fill(abs(eta2), dR)
			h_lead_vs_subleadeta_ER2.Fill(abs(eta1), abs(eta2))

		    if not ER2_flag :
			if pt_lead >= 6 and pt_sublead >= 4.:
				for dRcount in range(26) :
					dRthresh = 0.2*dRcount
					for etathresh in abs_eta_bins :
		
						ch1 = int(uGT_br.muonChg[leadid])
						ch2 = int(uGT_br.muonChg[subleadid])

#					print "leadid : ", leadid, " eta1 : ", eta1, " phi1 : ", phi1, " ptlead : ", pt_lead 
#					print "subleadid : ", subleadid, " eta2 : ", eta2, " phi2 : ", phi2, " ptsublead : ", pt_sublead
#					print "delta R : ", dR, "\n"


						if etathresh > 0 and dRthresh > 0 :
							if abs(eta1) < etathresh and abs(eta2) < etathresh and dR < dRthresh :
								h_rate_dR_vs_eta_dxy00.Fill(etathresh-0.05, dRthresh-0.05)
								if dxy1 >= 1 : h_rate_dR_vs_eta_dxy10.Fill(etathresh-0.05, dRthresh-0.05)
				
	
			for ptsubleadthresh in range(0, 40) :
				for ptleadthresh in range(ptsubleadthresh, 40) :
					if pt_sublead >= ptsubleadthresh and pt_lead >= ptleadthresh and abs(eta1) < 2. and abs(eta2) < 2.:
						h_rate_l1pt1_vs_l1pt2_dxy0.Fill(ptleadthresh, ptsubleadthresh)
						dxy1 = float(uGT_br.muonDxy[leadid])
						dxy2 = float(uGT_br.muonDxy[subleadid])
						if dxy1 >= 1 : h_rate_l1pt1_vs_l1pt2_dxy10.Fill(ptleadthresh, ptsubleadthresh)
						if dxy1 >= 1 and dxy2 >= 1 : h_rate_l1pt1_vs_l1pt2_dxy11.Fill(ptleadthresh, ptsubleadthresh)
						if dxy2 >= 2 : h_rate_l1pt1_vs_l1pt2_dxy20.Fill(ptleadthresh, ptsubleadthresh)

        ## End loop: for jEvt in range(chains['Unp'][iCh].GetEntries()):
    ## End loop: for iCh in range(len(chains['Unp'])):

    ########Finish loop over chains/Events

    print '\nFinished loop over chains'

    out_file.cd()

#    h_dimu_rate_ptVtx.Scale(scale[evtclassid]/iEvt)
#    h_dimu_rate_ptDisp.Scale(scale[evtclassid]/iEvt)
#    h_dimu_rate_ptOr.Scale(scale[evtclassid]/iEvt)
#    h_dimu_rate_ptAnd.Scale(scale[evtclassid]/iEvt)
#    h_dimu_rate_ptUniqueDisp.Scale(scale[evtclassid]/iEvt)
#    h_dimu_rate_ptUniqueVtx.Scale(scale[evtclassid]/iEvt)

#    h_dimu_rate_ptVtx.SetLineColor(R.kBlack)
#    h_dimu_rate_ptVtx.Write()

#    h_dimu_rate_ptDisp.SetLineColor(R.kRed)
#    h_dimu_rate_ptDisp.Write()

#    h_dimu_rate_ptOr.SetLineColor(R.kBlue)
#    h_dimu_rate_ptOr.Write()

#    h_dimu_rate_ptAnd.SetLineColor(6)
#    h_dimu_rate_ptAnd.Write()

#    h_dimu_rate_ptUniqueDisp.SetLineColor(8)
#    h_dimu_rate_ptUniqueDisp.Write()

#    h_dimu_rate_ptUniqueVtx.SetLineColor(9)
#    h_dimu_rate_ptUniqueVtx.Write()

#########################################################################

#    L1_SingleMu22_BMTF = h_kBMTF_singlemu_rate_ptVtx.GetBinContent(23)
#    L1_SingleMu22_BMTF_UPT = h_kBMTF_singlemu_rate_ptDisp.GetBinContent(23)

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
    print "%3d. %-50s : %-7d \t rate : %.3f"%(78, "L1_TripleMu_5_3_3", L1_TripleMu_5_3_3, L1_TripleMu_5_3_3*scale[evtclassid]/iEvtPU)
#    print "%3d. %-50s : %-7d \t rate : %.3f"%(, , , *scale[evtclassid]/iEvtPU)

    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu15_7_BMTF ", L1_DoubleMu15_7_BMTF, L1_DoubleMu15_7_BMTF*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu15_7_BMTF_UPT ", L1_DoubleMu15_7_BMTF_UPT, L1_DoubleMu15_7_BMTF_UPT*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu15_7_UPT ", L1_DoubleMu15_7_UPT, L1_DoubleMu15_7_UPT*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu15_7_UPT_SQ ", L1_DoubleMu_15_7_UPT_SQ, L1_DoubleMu_15_7_UPT_SQ*scale[evtclassid]/iEvtPU)
    print "%-50s : %-7d \t rate : %.3f"%("L1_DoubleMu_upt6ipMin1_upt4 ", L1_DoubleMu15_7_UPT_DXY1, L1_DoubleMu15_7_UPT_DXY1*scale[evtclassid]/iEvtPU)
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

############################### Single Muon #############################

    Scale = scale[evtclassid]/(1e3*iEvtPU)

#    h_EMTF_rate_vs_ptthreshold_L1dxy0.Scale(Scale)
 #   h_EMTF_rate_vs_ptthreshold_L1dxy1.Scale(Scale)
  #  h_EMTF_rate_vs_ptthreshold_L1dxy2.Scale(Scale)
#    h_EMTF_rate_vs_ptthreshold_L1dxy3.Scale(Scale)

 #   h_EMTF_rate_vs_ptthreshold_L1dxy1.SetLineColor(2)
 #   h_EMTF_rate_vs_ptthreshold_L1dxy2.SetLineColor(4)
 #   h_EMTF_rate_vs_ptthreshold_L1dxy3.SetLineColor(6)

#    h_EMTF_rate_vs_ptthreshold_L1dxy0.SetLineWidth(2)
#    h_EMTF_rate_vs_ptthreshold_L1dxy1.SetLineWidth(2)
#    h_EMTF_rate_vs_ptthreshold_L1dxy2.SetLineWidth(2)
#    h_EMTF_rate_vs_ptthreshold_L1dxy3.SetLineWidth(2)

#    h_EMTF_rate_vs_ptthreshold_L1dxy0.Write()
#    h_EMTF_rate_vs_ptthreshold_L1dxy1.Write()
#    h_EMTF_rate_vs_ptthreshold_L1dxy2.Write()
#    h_EMTF_rate_vs_ptthreshold_L1dxy3.Write()

#    h_rate_vs_l1pt2threshold4.Scale(Scale)
#    h_rate_vs_l1pt2threshold4.Write()
#    h_rate_vs_l1pt2threshold7.Scale(Scale)
#    h_rate_vs_l1pt2threshold7.Write()
#    h_rate_vs_l1pt2threshold11.Scale(Scale)
#    h_rate_vs_l1pt2threshold11.Write()

    h_rate_l1pt1_vs_l1pt2_dxy0.Scale(Scale)
    h_rate_l1pt1_vs_l1pt2_dxy0.Write()

    h_rate_l1pt1_vs_l1pt2_dxy10.Scale(Scale)
    h_rate_l1pt1_vs_l1pt2_dxy10.Write()

    h_rate_l1pt1_vs_l1pt2_dxy11.Scale(Scale)
    h_rate_l1pt1_vs_l1pt2_dxy11.Write()

    h_rate_l1pt1_vs_l1pt2_dxy20.Scale(Scale)
    h_rate_l1pt1_vs_l1pt2_dxy20.Write()

    h_rate_dR_vs_eta_dxy00.Scale(Scale)
    h_rate_dR_vs_eta_dxy00.Write()
    h_rate_dR_vs_eta_dxy10.Scale(Scale)
    h_rate_dR_vs_eta_dxy10.Write()

##########################

    h_dR_vs_leadeta_ER2.Write()
    h_dR_vs_subleadeta_ER2.Write()
    h_lead_vs_subleadeta_ER2.Write()

##########################

#    h_kBMTF_singlemu_rate_ptVtx.Scale(scale[evtclassid]/iEvt)
#    h_kBMTF_singlemu_rate_ptDisp.Scale(scale[evtclassid]/iEvt)
#    h_kBMTF_singlemu_rate_ptOr.Scale(scale[evtclassid]/iEvt)
#    h_kBMTF_singlemu_rate_ptDisp_Lxy2.Scale(scale[evtclassid]/iEvt)

#    h_kBMTF_singlemu_rate_ptVtx.SetLineColor(R.kBlack)
#    h_kBMTF_singlemu_rate_ptVtx.Write()

#    h_kBMTF_singlemu_rate_ptDisp.SetLineColor(R.kBlue)
#    h_kBMTF_singlemu_rate_ptDisp.Write()

#    h_kBMTF_singlemu_rate_ptOr.SetLineColor(R.kRed)
#    h_kBMTF_singlemu_rate_ptOr.Write()

#    h_kBMTF_singlemu_rate_ptDisp_Lxy2.SetLineColor(8)
#    h_kBMTF_singlemu_rate_ptDisp_Lxy2.Write()

#########################################################################

    out_file.Close()
    del chains

#def BaselineRun(uGT_br, i, j):
#	ptVtx1 	= float(uGT_br.muonEt[i])
#	ptVtx2 	= float(uGT_br.muonEt[j])
	
#	pt1 = max(ptVtx1, ptVtx2)
#	pt2 = min(ptVtx1, ptVtx2)

if __name__ == '__main__':
    main()
