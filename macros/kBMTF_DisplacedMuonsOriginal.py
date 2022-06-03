#! /usr/bin/env python

## **************************************************************** ##
##  Look at properties of displaced muons from Kalman algo in BMTF  ##
## **************************************************************** ##

import os, sys
from math import *
from array import array

import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

style1 = R.TStyle("style1", "For Histograms")
style1.SetLineWidth(2)
style1.SetOptStat(0)

R.gROOT.SetStyle("style1")

PRT_EVT  = 10000	 	## Print every Nth event
MAX_EVT  = 40000	## Number of events to process
VERBOSE  = False	## Verbose print-out

def getPhi(globalPhiHw) :
    phi = globalPhiHw/287.5*pi if globalPhiHw<287.5 else \
                 (globalPhiHw-575.)/287.5*pi
    return phi

def getGenEtaPhi(Gen_br, i, mutype):
	if mutype == 1 : return getGenEtaPhiBarrel(Gen_br, i)
	elif mutype == 3 : return getGenEtaPhiEndcap(Gen_br, i)
	else  : return False, False

def getGenEtaPhiBarrel(Gen_br, i):

	eta 	= float(Gen_br.partEta[i])
	vx	= float(Gen_br.partVx[i])
	vy	= float(Gen_br.partVy[i])
	vz	= float(Gen_br.partVz[i])
	phi 	= float(Gen_br.partPhi[i])

	Lxy = sqrt(float(Gen_br.partVx[i])**2 + float(Gen_br.partVy[i])**2)	

	r = (512. - Lxy)
	z = r*sinh(eta)

	zStar = vz + z
	xStar = vx + r*cos(phi)
	yStar = vy + r*sin(phi)
 
	rStar = sqrt(xStar**2 + yStar**2)
	
	GenEta = asinh(zStar/512.)
	
	if (xStar >= 0) : GenPhi = atan(yStar/xStar)
      	elif (yStar >= 0 and xStar < 0) : GenPhi = pi + atan(yStar/xStar)
      	elif (yStar <= 0 and xStar < 0) : GenPhi = atan(yStar/xStar) - pi 	
	
	return GenEta, GenPhi

def IsinacceptBarrel(Gen_br, i):

	vx = float(Gen_br.partVx[i])
	vy = float(Gen_br.partVy[i])
	vz = float(Gen_br.partVz[i])

	eta = float(Gen_br.partEta[i])
	
#	if abs(vz) >= 650. : return False
	
	Lxy = sqrt(vx**2 + vy**2)
	if Lxy > 600. : return False

	if vz == 650. : theta1 = pi/2.
	else: theta1 = atan((700.-Lxy)/(650.-vz))
	if theta1 < 0 : theta1 = pi+theta1

	if vz == -650. : theta2 = pi/2.
	else: theta2 = pi-atan((700.-Lxy)/(650.+vz))
	if theta2 > pi : theta2 = theta2 - pi

	maxeta = -log(tan(0.5*theta1))
	mineta = -log(tan (0.5*theta2))

	if eta >= maxeta or eta <= mineta : return False
	
	return True

def getGenEtaPhiEndcap(Gen_br, i):				
	
	z_ME2 	= 800.
	eta 	= float(Gen_br.partEta[i])
	phi 	= float(Gen_br.partPhi[i])
	z	= float(Gen_br.partVz[i])	

	if (eta > 0) 	: r = abs(z_ME2 - z)/abs(sinh(eta))
	else	    	: r = abs(- z_ME2 - z)/abs(sinh(eta))

	xStar = float(Gen_br.partVx[i]) + r*cos(phi)
	yStar = float(Gen_br.partVy[i]) + r*sin(phi)
	rStar = sqrt(xStar**2 + yStar**2)

	GenEta	= asinh(z_ME2/rStar)*(eta/abs(eta))

	if (xStar >= 0) : GenPhi = atan(yStar/xStar)
      	elif (yStar >= 0 and xStar < 0) : GenPhi = pi + atan(yStar/xStar)
      	elif (yStar <= 0 and xStar < 0) : GenPhi = atan(yStar/xStar) - pi 	

	return GenEta, GenPhi

def IsinacceptEndcap(Gen_br, i):

	vx = float(Gen_br.partVx[i])
	vy = float(Gen_br.partVy[i])
	vz = float(Gen_br.partVz[i])

	if abs(vz) >=100. : return False

	Lxy = sqrt(vx**2 + vy**2)
	if Lxy > 490. : return False

	eta, phi = getGenEtaPhiEndcap(Gen_br, i)

	if abs(eta) <= 1.24 : return False
	if abs(eta) >= 2.5 : return False	

	return True

def getDxy(vx, vy, phi):
	return abs(vx*sin(phi)-vy*cos(phi))


##################################################################################################################


def main():

    print '\nInside DisplacedMuons\n'
    evtclass = ["signal","DisplacedMuGun"]
    evtclassid = 0
    inputdir = ['/eos/user/s/sonawane/temp/L1Ntuples/signal_tuples/ntuples_01_07_21/backup/',
		'/eos/user/s/sonawane/temp/L1Ntuples/Displaced_mu_gun_tuples/ntuples_16_06_21/']

    workdir = '/afs/cern.ch/user/s/sonawane/L1T/L1studies/L1_scripts_Alberto/L1RunIII/macros/'
    in_file_names = []

    if evtclassid == 1:
	for s in ['DisplacedMuGun_Pt2to10_11_2_X_1623847533/', 'DisplacedMuGun_Pt10to30_11_2_X_1623847663/', 'DisplacedMuGun_Pt30to100_11_2_X_1623847721/']:
#	for s in ['MuGun_Pt2to10_Nu_11_2_X_1623157277/']:
		ntupledir = inputdir[evtclassid]+s
		for i in range(23) :
        		path = ntupledir+str(i)+".root"
        		if not os.path.exists(path): continue
			in_file_names.append(path)

    elif evtclassid == 0:
        filestr = ['HTo2LongLivedTo4mu_MH-125_MFF-12_CTau-900mm_11_2_X_1625152032/', 'HTo2LongLivedTo4mu_MH-125_MFF-25_CTau-1500mm_11_2_X_1625152298/', 'HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_11_2_X_1623847924/']
        samplename = ["125_12_900", "125_25_1500", "125_50_3000"]
        for l in range(2,3) :
                s = filestr[l]
                sn = samplename[l]
                ntupledir = inputdir[evtclassid]+s
                for i in range(14) :
                        path = ntupledir+str(i)+".root"
                        if not os.path.exists(path): continue
                        in_file_names.append(path)			

    else:
	for i in range(14):
	        path = inputdir[evtclassid]+str(i)+".root"
       	 	if not os.path.exists(path): continue
		in_file_names.append(path)
    
    if not os.path.exists(workdir+'plots'): os.makedirs(workdir+'plots')

#	For testing
#    in_file_names=['/afs/cern.ch/user/s/sonawane/L1T/L1studies/L1_scripts_Alberto/L1RunIII/test/L1Ntuple.root']

    out_file_str = 'DisplacedMuons_'+evtclass[evtclassid]
    out_file_str += ('_%dk' % (MAX_EVT / 1000))
    out_file = R.TFile(workdir+'plots/'+out_file_str+'.root','recreate')

    chains = {}
    chains['Evt'] = []  ## Event info
    chains['Unp'] = []  ## Unpacked legacy BMTF
    chains['Emu'] = []  ## Emulated Kalman BMTF
    chains['EmuK'] = []  ## Emulated Kalman BMTF
    chains['EmuE'] = []  ## Emulated Kalman BMTF
    chains['uGT'] = []  ## Global Trigger
    chains['Gen'] = []  ## Generator information
    chains['Evt'].append( R.TChain('l1EventTree/L1EventTree') )
    chains['Unp'].append( R.TChain('l1UpgradeTfMuonTree/L1UpgradeTfMuonTree') )
    chains['Emu'].append( R.TChain('l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree') )
    chains['EmuK'].append( R.TChain('l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree') )
    chains['EmuE'].append( R.TChain('l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree') )
    chains['uGT'].append( R.TChain('l1UpgradeEmuTree/L1UpgradeTree') )
    chains['Gen'].append( R.TChain('l1GeneratorTree/L1GenTree') )

    for i in range(len(in_file_names)):
        print 'Adding file %s' % in_file_names[i]
        chains['Evt'][0].Add( in_file_names[i] )
        chains['Unp'][0].Add( in_file_names[i] )
        chains['Emu'][0].Add( in_file_names[i] )
        chains['EmuK'][0].Add( in_file_names[i] )
        chains['EmuE'][0].Add( in_file_names[i] )
        chains['uGT'][0].Add( in_file_names[i] )
        chains['Gen'][0].Add( in_file_names[i] )


    ###################
    ### Book histograms
    ###################

    pt_bins  = [150, 0, 300]
    q_pt_bins  = [100, -0.2, 0.2]
#    chi_bins = [100, 0, 100]
    dxy_bins = [10,   0, 150]
    vxyz_bins = [200,   -200, 200]
    lxy_bins = [30,   0, 150]
    phi_bins = [80, -4.,4.]
    eta_bins = [30, -3., 3.]
    dphi_bins = [40, 0., 4.]
    deta_bins = [30, 0., 3.]
#    Vx_bins  = [100,-50, 50]
#    phi_bins = [80, -4, 4]
    EC_eta_bins = [-3., -2.5, -2.1, -1.6, -1.2, 1.2, 1.6, 2.1, 2.5, 3.]

### Single Muon Efficiencies ########

    h_pt_blank = R.TH1F('h_pt_blank', 'Blank Gen Pt; Gen Muon Pt [GeV]; L1T Efficiency', pt_bins[0], pt_bins[1], pt_bins[2])
    h_dxy_blank = R.TH1F('h_dxy_blank', 'Blank Gen Dxy; Gen Muon Dxy [cm]; L1T Efficiency', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_EC_eta_blank = R.TH1F('h_EC_eta_blank', 'Blank Gen Eta; Gen Muon Eta; L1T Efficiency', len(EC_eta_bins)-1, array('d',EC_eta_bins))

    h_Eff_vs_dxy_l1pt0_den = R.TH1F('h_Eff_vs_dxy_l1pt0_den', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_Eff_vs_dxy_l1pt10_den = R.TH1F('h_Eff_vs_dxy_l1pt10_den', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_Eff_vs_dxy_l1pt20_den = R.TH1F('h_Eff_vs_dxy_l1pt20_den', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_Eff_vs_dxy_l1pt30_den = R.TH1F('h_Eff_vs_dxy_l1pt30_den', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])

    h_Eff_vs_dxy_l1pt0_num  = R.TH1F('h_Eff_vs_dxy_l1pt0_num', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_Eff_vs_dxy_l1pt10_num = R.TH1F('h_Eff_vs_dxy_l1pt10_num', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_Eff_vs_dxy_l1pt20_num = R.TH1F('h_Eff_vs_dxy_l1pt20_num', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_Eff_vs_dxy_l1pt30_num = R.TH1F('h_Eff_vs_dxy_l1pt30_num', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])

    genptplus = 10.

#Event counters

    GenMuevt_count 		= 0
    
    GenMuCount 			= 0
    GenMuinAcc			= 0

    EmuMuevt_count	 	= 0
    EmuMus_count		= 0
    EmuMus_unique_count 	= 0

    matched_gen_mu_count	= 0
    matched_gen_mu_evt		= 0

    mismatch_evt		= 0

    iEvt 			= 0 

    ##### Dimuon counters #####
   
    GenDimu_Evt = 0

    print '\nEntering loop over chains'
    for iCh in range(len(chains['Emu'])):

        if iEvt >= MAX_EVT: break

        ## Faster tecnhique, inspired by https://github.com/thomreis/l1tMuonTools/blob/master/L1Analysis.py
        Evt_br = R.L1Analysis.L1AnalysisEventDataFormat()
        Unp_br = R.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
        EmuK_br = R.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
        EmuE_br = R.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
        uGT_br = R.L1Analysis.L1AnalysisL1UpgradeDataFormat()
        Gen_br = R.L1Analysis.L1AnalysisGeneratorDataFormat()
#        Kmt_br = R.L1Analysis.L1AnalysisBMTFOutputDataFormat()

        chains['Evt'][iCh].SetBranchAddress('Event',               R.AddressOf(Evt_br))
        chains['Unp'][iCh].SetBranchAddress('L1UpgradeBmtfMuon',   R.AddressOf(Unp_br))
        chains['EmuK'][iCh].SetBranchAddress('L1UpgradeKBmtfMuon', R.AddressOf(EmuK_br))
        chains['EmuE'][iCh].SetBranchAddress('L1UpgradeEmtfMuon',  R.AddressOf(EmuE_br))
        chains['uGT'][iCh].SetBranchAddress('L1Upgrade',  	   R.AddressOf(uGT_br))
        chains['Gen'][iCh].SetBranchAddress('Generator',           R.AddressOf(Gen_br))
#        chains['Emu'][iCh].SetBranchAddress('L1UpgradeBmtfOutput', R.AddressOf(Kmt_br))

	

        print '\nEntering loop over events for chain %d' % iCh
        for jEvt in range(chains['Emu'][iCh].GetEntries()):

            if iEvt >= MAX_EVT: break

	    iEvt +=1

            if iEvt % PRT_EVT is 0: print '\nEvent # %d (%dth in chain)' % (iEvt, jEvt+1)

            chains['Evt'][iCh].GetEntry(jEvt)
            chains['Unp'][iCh].GetEntry(jEvt)
            chains['EmuK'][iCh].GetEntry(jEvt)
            chains['EmuE'][iCh].GetEntry(jEvt)
            chains['uGT'][iCh].GetEntry(jEvt)
            chains['Gen'][iCh].GetEntry(jEvt)

            if iEvt % PRT_EVT is 0: print '  * Run %d, LS %d, event %d' % (int(Evt_br.run), int(Evt_br.lumi), int(Evt_br.event))

            nUnpMu = int(Unp_br.nTfMuons)
            nEmuKMu = int(EmuK_br.nTfMuons)
            nEmuEMu = int(EmuE_br.nTfMuons)
            nuGTMu = int(uGT_br.nMuons)
	    nGenPart = int(Gen_br.nPart)
#            nKmtMu = int(Kmt_br.nTrks)
	    
	    ##########################
            ###  Index containers  ###
            ##########################
	    
  	    GenMus			=[]
	    GenMus_acc	 		=[]
	    EmuMus			=[]
	    EmuMus_unique		=[]
	    GenMuPt_acc			=[]

	    EmuKMus			=[]
	    EmuKMus_unique		=[]
	
	    #########################################
            ###  Generator information for muons  ###
            #########################################

	    for i in range(nGenPart):
	
		if (abs(Gen_br.partId[i])!=13): continue
		if (Gen_br.partStat[i]!=1): continue

		### weighting MuGun events ###

		GenMus.append(i)

		if IsinacceptBarrel(Gen_br, i) : MuType = 1 
		#if IsinacceptEndcap(Gen_br, i) : MuType = 3
		else : 	continue

		GenMus_acc.append([i, MuType])

                gen_eta, gen_phi = getGenEtaPhi(Gen_br, i, MuType)
		GenDxy  = getDxy(float(Gen_br.partVx[i]), float(Gen_br.partVy[i]), float(Gen_br.partPhi[i]))
                z0      = float(Gen_br.partVz[i])
		ptgen 	= float(Gen_br.partPt[i])

		if ptgen > 0. + genptplus : h_Eff_vs_dxy_l1pt0_den.Fill(GenDxy)
		if ptgen > 10. + genptplus : h_Eff_vs_dxy_l1pt10_den.Fill(GenDxy)
		if ptgen > 20. + genptplus : h_Eff_vs_dxy_l1pt20_den.Fill(GenDxy)
		if ptgen > 30. + genptplus : h_Eff_vs_dxy_l1pt30_den.Fill(GenDxy)

		GenMuPt_acc.append([Gen_br.partPt[i], i, MuType])

	    GenMuPt_acc.sort(reverse=True)

            #################################
            ###  Emulated uGT muons  ###
            #################################

            for i in range(nuGTMu):
                BX      = int(uGT_br.muonBx[i])
                                
                if (BX  !=  0): continue
#                if (qual < 11): continue

		EmuMus.append(i)

	    Eta_phi = []

            for i in EmuMus :

                pt1  = float(uGT_br.muonEtUnconstrained[i])
                eta1 = float(uGT_br.muonEta[i])
                phi1 = float(uGT_br.muonPhi[i])

                ep = [eta1, phi1]

                pt_i = [pt1, i]

                ptmax = pt1

                if ep in Eta_phi : continue

                Eta_phi.append(ep)

                for j in EmuMus :
                        if j <= i : continue

                        pt2 = float(uGT_br.muonEtUnconstrained[j])
                        eta2 = float(uGT_br.muonEta[j])
                        phi2 = float(uGT_br.muonPhi[j])

                        if eta1==eta2 and phi1==phi2 :
                                if pt2 > ptmax :
                                        ptmax = pt2
                                        pt_i = [pt2, j]

                EmuMus_unique.append(pt_i[1])


            #################################
            ###  Emulated EMTF muons  ###
            #################################

            for i in range(nEmuKMu):
                BX      = int(EmuK_br.tfMuonBx[i])
                                
                if (BX  !=  0): continue
#               if (qual < 11): continue

		EmuKMus.append(i)

	    Eta_phi = []

            for i in EmuKMus :

                pt1  = float(EmuK_br.tfMuonHwPtUnconstrained[i]-1.)
                eta1 = float(EmuK_br.tfMuonHwEta[i])
                phi1 = getPhi(float(EmuK_br.tfMuonGlobalPhi[i]))

                ep = [eta1, phi1]

                pt_i = [pt1, i]

                ptmax = pt1

                if ep in Eta_phi : continue

                Eta_phi.append(ep)

                for j in EmuKMus :
                        if j <= i : continue

                        pt2 = float(EmuK_br.tfMuonHwPtUnconstrained[j]-1.)
                        eta2 = float(EmuK_br.tfMuonHwEta[j])
                        phi2 = getPhi(float(EmuK_br.tfMuonGlobalPhi[j]))

                        if eta1==eta2 and phi1==phi2 :
                                if pt2 > ptmax :
                                        ptmax = pt2
                                        pt_i = [pt2, j]

                EmuKMus_unique.append(pt_i[1])


	    # Better matching subroutine

	    dR_list_corr=[]

	    EMTF_dR_list=[]
	    kBMTF_dR_list=[]

	    for i in EmuKMus:

		EmuMuvec = R.TLorentzVector()
		EmuMuvec.SetPtEtaPhiM(10., float(EmuK_br.tfMuonHwEta[i])*0.010875, getPhi(float(EmuK_br.tfMuonGlobalPhi[i])), 105.7e-3)

		for el in GenMus_acc :
			j = el[0]
			mutype = el[1]
		
			GenMuCorr = R.TLorentzVector()
			temp_eta, temp_phi = getGenEtaPhi(Gen_br, j, mutype)
			GenMuCorr.SetPtEtaPhiM(float(Gen_br.partPt[j]), temp_eta, temp_phi, 105.7e-3)
			eta1corr = GenMuCorr.Eta()
			phi1corr = GenMuCorr.Phi()
			dR = EmuMuvec.DeltaR(GenMuCorr)
			if dR < 0.5: dR_list_corr.append([dR, i, j])


	    dR_list = dR_list_corr

	    if len(dR_list):

		    dR_list.sort()

		    k = True 
		    i = 0
		    while k:
			k = False
			idx1 = dR_list[i][1]
			idx2 = dR_list[i][2]
			j = i+1
			while j < len(dR_list):
				idx_1  = dR_list[j][1]
				idx_2  = dR_list[j][2]
	
				if idx1 == idx_1 or idx2 == idx_2: del dR_list[j]
				else : j +=1
		  	i +=1
			if i < len(dR_list): 
				k = True

####### Dxy turnon for kBMTF, GMT and EMTF #######

####### Efficiencies for kBMTF ############
		
	    for i in range(len(dR_list)) :
		Emu_idx = dR_list[i][1]
		Gen_idx = dR_list[i][2]
		L1Dxy	= float(EmuK_br.tfMuonHwDxy[Emu_idx])
		GenDxy  = getDxy(float(Gen_br.partVx[Gen_idx]), float(Gen_br.partVy[Gen_idx]), float(Gen_br.partPhi[Gen_idx]))
		z0 	= float(Gen_br.partVz[Gen_idx])
		ptgen 	= float(Gen_br.partPt[Gen_idx])
#		dR	= dR_list[i][0] 
		ptVtx 	= float(EmuK_br.tfMuonHwPt[Emu_idx]-1.)*0.5
		ptDisp 	= float(EmuK_br.tfMuonHwPtUnconstrained[Emu_idx]-1.)
		gen_eta, gen_phi = getGenEtaPhi(Gen_br, Gen_idx, mutype)
		emu_eta = float(EmuK_br.tfMuonHwEta[Emu_idx])*0.010875

		pt = ptDisp


		if pt > 0. and ptgen > 0. + genptplus		: 
			h_Eff_vs_dxy_l1pt0_num.Fill(GenDxy)

		if pt > 10. and ptgen > 10. + genptplus :
			h_Eff_vs_dxy_l1pt10_num.Fill(GenDxy)

		if pt > 20. and ptgen > 20. + genptplus :
			h_Eff_vs_dxy_l1pt20_num.Fill(GenDxy)

		if pt > 30. and ptgen > 30. + genptplus :
			h_Eff_vs_dxy_l1pt30_num.Fill(GenDxy)

##################################################

	    ##End of sub loop over Gen, Emu or Unp muons

#	    print "GenMus: ", len(GenMus)
#	    print "GenMusAcc: ", len(GenMus_acc)
#	    print "EmuMus: ", len(EmuMus)
#	    print "EmuMus unique: ", len(EmuMus_unique)
	
	    ###########################
            ###  Updating counters  ###
            ###########################
 	    
	    if len(GenMus) : GenMuevt_count +=1
	    GenMuCount += len(GenMus)
	    GenMuinAcc += len(GenMus_acc)

	    if len(EmuMus) : EmuMuevt_count +=1
	    EmuMus_count += len(EmuMus)

	    EmuMus_unique_count +=len(EmuMus_unique)

	    if len(dR_list) : matched_gen_mu_evt +=1

        ## End loop: for jEvt in range(chains['Unp'][iCh].GetEntries()):
    ## End loop: for iCh in range(len(chains['Unp'])):

    print '\nFinished loop over chains'

    out_file.cd()

    h_pt_blank.Write()
    h_EC_eta_blank.Write()
    h_dxy_blank.Write()

    h_Teff_vs_dxy_l1pt0 = R.TEfficiency(h_Eff_vs_dxy_l1pt0_num, h_Eff_vs_dxy_l1pt0_den)
    h_Teff_vs_dxy_l1pt10 = R.TEfficiency(h_Eff_vs_dxy_l1pt10_num, h_Eff_vs_dxy_l1pt10_den)
    h_Teff_vs_dxy_l1pt20 = R.TEfficiency(h_Eff_vs_dxy_l1pt20_num, h_Eff_vs_dxy_l1pt20_den)
    h_Teff_vs_dxy_l1pt30 = R.TEfficiency(h_Eff_vs_dxy_l1pt30_num, h_Eff_vs_dxy_l1pt30_den)

    h_Teff_vs_dxy_l1pt0.SetName('h_Teff_vs_dxy_l1pt0')
    h_Teff_vs_dxy_l1pt10.SetName('h_Teff_vs_dxy_l1pt10')
    h_Teff_vs_dxy_l1pt20.SetName('h_Teff_vs_dxy_l1pt20')
    h_Teff_vs_dxy_l1pt30.SetName('h_Teff_vs_dxy_l1pt30')


    h_Teff_vs_dxy_l1pt10.SetLineColor(R.kBlue)
    h_Teff_vs_dxy_l1pt20.SetLineColor(R.kRed)
    h_Teff_vs_dxy_l1pt30.SetLineColor(R.kGreen)

    h_Teff_vs_dxy_l1pt0.Write()
    h_Teff_vs_dxy_l1pt10.Write()
    h_Teff_vs_dxy_l1pt20.Write()
    h_Teff_vs_dxy_l1pt30.Write()

    out_file.Close()
    del chains

    print '\nWrote out file: plots/'+out_file_str+'.root'
    print '\n Events run over: ', iEvt

    print 'Event counter'
    print 'Total GenMus: ', GenMuCount
    print 'GenMus in Acceptance: ', GenMuinAcc
    print 'EmuMus passing qual: ', EmuMus_count
    print 'EmuMus passing qual and unique: ', EmuMus_unique_count
    
    print 'Mismatch between corr and uncorr Gen Muons matched to qual emus: ', mismatch_evt

    print 'Matched gen mus : ', matched_gen_mu_count
    print 'Events with matched gen Mus: ', matched_gen_mu_evt

    print 'Events with Gen Dimuons in acceptance: ', GenDimu_Evt

if __name__ == '__main__':
    main()

