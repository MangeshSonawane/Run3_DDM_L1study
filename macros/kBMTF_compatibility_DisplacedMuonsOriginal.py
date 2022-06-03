#! /usr/bin/env python

## **************************************************************** ##
##  Look at properties of displaced muons from Kalman algo in BMTF  ##
## **************************************************************** ##

import os

import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

PRT_EVT  = 10000   ## Print every Nth event
MAX_EVT  = 1000000 ## Number of events to process
VERBOSE  = False  ## Verbose print-out

style1 = R.TStyle("style1", "For Histograms")
style1.SetLineWidth(2)
style1.SetOptStat(0)

R.gROOT.SetStyle("style1")

def main():

    print '\nInside DisplacedMuons\n'
    evtclass = ["NuGun_rate", "OldNuGun_rate", "signal_3000_rate"]
    evtclassid = 1
    inputdir = ['/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/elfontan/condor/reMu_reHcalTP_PFA1p_v15_LUTGenTrue_kBMTFghostbusting_uptScales_newNuGun_tagv16/', '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/stempl/condor/menu_Nu_11_0_X_1614189426/','/eos/user/s/sonawane/temp/L1Ntuples/signal_3000_tuples/ntuples_15_06_21/HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_11_2_X_1623671915/']
    workdir = '/afs/cern.ch/user/s/sonawane/L1T/L1studies/L1_scripts_Alberto/L1RunIII/macros/'
    in_file_names = []
    for i in range(88):
	if i == 44 : continue
        path = inputdir[evtclassid]+str(i)+".root"
        if not os.path.exists(path): continue
	in_file_names.append(path)
    
    scale = [2544*11246., 2544*11246., 1.]

    if not os.path.exists(workdir+'plots'): os.makedirs(workdir+'plots')

    out_file_str = 'DisplacedMuons_'+evtclass[evtclassid]
    out_file_str += ('_%dk' % (MAX_EVT / 1000))
    out_file = R.TFile(workdir+'plots/'+out_file_str+'.root','recreate')

    chains = {}
    chains['Evt'] = []  ## Event info
    chains['Emu'] = []  ## Emulated Kalman BMTF
    chains['uGT'] = []  ## Global Trigger
    for i in range(len(in_file_names)):
        print 'Adding file %s' % in_file_names[i]
        chains['Evt'].append( R.TChain('l1EventTree/L1EventTree') )
        chains['Emu'].append( R.TChain('l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree') )
        chains['uGT'].append( R.TChain('l1UpgradeEmuTree/L1UpgradeTree') )
        chains['Evt'][i].Add( in_file_names[i] )
        chains['Emu'][i].Add( in_file_names[i] )
        chains['uGT'][i].Add( in_file_names[i] )


# Rate histograms
    h_dimu_rate_ptVtx 	=	R.TH1F('h_dimu_rate_ptVtx', 	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]',	150, 0, 150)
    h_dimu_rate_ptDisp 	=	R.TH1F('h_dimu_rate_ptDisp', 	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]', 	150, 0, 150)
    h_dimu_rate_ptOr 	=	R.TH1F('h_dimu_rate_ptOr', 	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]', 	150, 0, 150)

    h_dimu_rate_ptAnd		= 	R.TH1F('h_dimu_rate_ptAnd',	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]',	150, 0, 150)
    h_dimu_rate_ptUniqueDisp	= 	R.TH1F('h_dimu_rate_ptUniqueDisp',	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]',	150, 0, 150)
    h_dimu_rate_ptUniqueVtx	= 	R.TH1F('h_dimu_rate_ptUniqueVtx',	'Dimuon Pt rate ; L1 lead Mu Pt; Rate [Hz]',	150, 0, 150)

    h_kBMTF_singlemu_rate_ptVtx		= 	R.TH1F('h_kBMTF_singlemu_rate_ptVtx',	'Dimuon Pt rate ; L1 Mu Pt; Rate [Hz]',	150, 0, 150)
    h_kBMTF_singlemu_rate_ptDisp	= 	R.TH1F('h_kBMTF_singlemu_rate_ptDisp',	'Dimuon Pt rate ; L1 Mu Pt; Rate [Hz]',	150, 0, 150)
    h_kBMTF_singlemu_rate_ptOr		= 	R.TH1F('h_kBMTF_singlemu_rate_ptOr',	'Dimuon Pt rate ; L1 Mu Pt; Rate [Hz]',	150, 0, 150)
    h_kBMTF_singlemu_rate_ptDisp_Lxy2	= 	R.TH1F('h_kBMTF_singlemu_rate_ptDisp_Lxy2',	'Dimuon Pt rate ; L1 Mu Pt; Rate [Hz]',	150, 0, 150)



    ########################
    ### Seed event counter
    ########################

    L1_SingleMu7 	= 0
    L1_SingleMu22 	= 0
    L1_SingleMu25 	= 0
    L1_DoubleMu15_7	= 0
    L1_DoubleMu18_er2p1	= 0

    L1_DoubleMu15_7_BMTF	= 0
    L1_DoubleMu15_7_BMTF_UPT	= 0

    iEvt = 0 

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
        Emu_br = R.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
        uGT_br = R.L1Analysis.L1AnalysisL1UpgradeDataFormat()
#        Kmt_br = R.L1Analysis.L1AnalysisBMTFOutputDataFormat()

        chains['Evt'][iCh].SetBranchAddress('Event',               R.AddressOf(Evt_br))
        chains['Emu'][iCh].SetBranchAddress('L1UpgradeKBmtfMuon',   R.AddressOf(Emu_br))
        chains['uGT'][iCh].SetBranchAddress('L1Upgrade',   	   R.AddressOf(uGT_br))
#        chains['Emu'][iCh].SetBranchAddress('L1UpgradeBmtfOutput', R.AddressOf(Kmt_br))


        print '\nEntering loop over events for chain %d' % iCh
        for jEvt in range(chains['Emu'][iCh].GetEntries()):
	    
            if iEvt >= MAX_EVT: break
	    
            iEvt +=1

            if iEvt % PRT_EVT is 0: print '\nEvent # %d (%dth in chain)' % (iEvt, jEvt+1)

            chains['Evt'][iCh].GetEntry(jEvt)
            chains['Emu'][iCh].GetEntry(jEvt)
            chains['uGT'][iCh].GetEntry(jEvt)

            # ## Use these lines if you don't explicitly define the DataFormat and then do SetBranchAddress above
            # Evt_br = chains['Evt'][iCh].Event
            # Unp_br = chains['Unp'][iCh].L1UpgradeBmtfMuon
            # Emu_br = chains['Emu'][iCh].L1UpgradeBmtfMuon

            if iEvt % PRT_EVT is 0: print '  * Run %d, LS %d, event %d' % (int(Evt_br.run), int(Evt_br.lumi), int(Evt_br.event))

            nEmuMu = int(Emu_br.nTfMuons)
            nuGTMu = int(uGT_br.nMuons)
#            nKmtMu = int(Kmt_br.nTrks)

#            if (VERBOSE and nUnpMu > 0 or nEmuMu > 0):
#                print 'Unpacked = %d, emulated = %d (internal %d) total muons in collection' % (nUnpMu, nEmuMu, nKmtMu)
#                print 'Unpacked = %d, emulated = %d total muons in collection' % (nUnpMu, nEmuMu)
        

	    uGT_mus = []

	    L1_SingleMu7_flag = False
	    L1_SingleMu22_flag = False
	    L1_SingleMu25_flag = False
	    L1_DoubleMu15_7_flag = False
	    L1_DoubleMu18_er2p1_flag = False

############ new seeds #################

	    L1_DoubleMu15_7_BMTF_flag = False
	    L1_DoubleMu15_7_BMTF_UPT_flag = False

########################################

            #################################
            ###  Emulated (Global) muons  ###
            #################################
            for i in range(nEmuMu):
                BX      = int(Emu_br.tfMuonBx[i])
                ptVtx   = float(Emu_br.tfMuonHwPt[i]-1.)*0.5  
                ptDisp  = float(Emu_br.tfMuonHwPtUnconstrained[i]-1.)  
                eta     = float(Emu_br.tfMuonHwEta[i])*0.010875
#                if VERBOSE: print 'Emulated Global muon %d BX = %d, qual = %d, ptVtx = %.1f, ptDispl = %.1f, eta = %.2f' % (i, BX, qual, ptVtx, ptDispl, eta)
                
                if (BX  !=  0): continue

		uGT_mus.append(i)
		
	    pt1 = -10. 
	    pt2 = -10.

	    for i in uGT_mus :
	
                qual 	= int(Emu_br.tfMuonHwQual[i])
		ptVtx	= float(Emu_br.tfMuonHwPt[i]-1.)*0.5
		eta 	= float(Emu_br.tfMuonHwEta[i])*0.010875
		ptDisp	= float(Emu_br.tfMuonHwPtUnconstrained[i]-1.)
		Lxy 	= float(Emu_br.tfMuonHwDxy[i])

		if ptVtx>=7. and qual in MU_QLTY_SNGL	: L1_SingleMu7_flag = True
		if ptVtx>=22. and qual in MU_QLTY_SNGL	: L1_SingleMu22_flag = True
		if ptVtx>=25. and qual in MU_QLTY_SNGL	: L1_SingleMu25_flag = True

		if abs(eta) < 0.8 :
			for j in range(150):
				ptthresh = float(j)
				if ptVtx >= ptthresh and qual in MU_QLTY_SNGL : h_kBMTF_singlemu_rate_ptVtx.Fill(ptthresh+0.1)
				if ptDisp >= ptthresh and qual in MU_QLTY_SNGL : h_kBMTF_singlemu_rate_ptDisp.Fill(ptthresh+0.1)
				if (ptVtx >= ptthresh or ptDisp ptthresh) and qual in MU_QLTY_SNGL	: h_kBMTF_singlemu_rate_ptOr.Fill(ptthresh+0.1)
				if ptDisp >= ptthresh and Lxy >=2 and qual in MU_QLTY_SNGL	: h_kBMTF_singlemu_rate_ptDisp_Lxy2.Fill(ptthresh+0.1) 
	
		for j in uGT_mus :
			if i == j : continue
			qual2    = int(Emu_br.tfMuonHwQual[j])
        	        ptVtx2   = float(Emu_br.tfMuonHwPt[j]-1.)*0.5
	                ptDisp2  = float(Emu_br.tfMuonHwPtUnconstrained[j]-1.)
			eta2	 = float(Emu_br.tfMuonHwEta[j])*0.010875 

			if (qual in MU_QLTY_DBLE) and (qual2 in MU_QLTY_DBLE):
				pt1 = max(ptVtx, ptVtx2)
				pt2 = min(ptVtx, ptVtx2)
				if pt1 >= 15. and pt2 >= 7.	: L1_DoubleMu15_7_flag = True
				if (pt1 >= 15. and abs(eta) < 0.8) and (pt2 >= 7. and abs(eta2) < 0.8)	: L1_DoubleMu15_7_BMTF_flag = True
 
				ptdisp1 = max(ptDisp, ptDisp2)
				ptdisp2 = min(ptDisp, ptDisp2)

				if (ptdisp1 >= 15. and abs(eta) < 0.8) and (ptdisp2 >= 7. and abs(eta2) < 0.8)	: L1_DoubleMu15_7_BMTF_UPT_flag = True

				for l in range(150) :
					ptthresh = float(l)
					if pt1 >= ptthresh and pt2 >=7. : h_dimu_rate_ptVtx.Fill(ptthresh)
					if ptdisp1 >= ptthresh and ptdisp2 >=7. : h_dimu_rate_ptDisp.Fill(ptthresh)
					if (pt1 >= ptthresh and pt2 >=7.) or (ptdisp1 >= ptthresh and ptdisp2 >=7.) : h_dimu_rate_ptOr.Fill(ptthresh)
					if (pt1 >= ptthresh and pt2 >=7.) and (ptdisp1 >= ptthresh and ptdisp2 >=7.) : h_dimu_rate_ptAnd.Fill(ptthresh)
					if not (pt1 >= ptthresh and pt2 >=7.) and (ptdisp1 >= ptthresh and ptdisp2 >=7.) : h_dimu_rate_ptUniqueDisp.Fill(ptthresh)
					if (pt1 >= ptthresh and pt2 >=7.) and not (ptdisp1 >= ptthresh and ptdisp2 >=7.) : h_dimu_rate_ptUniqueVtx.Fill(ptthresh)

			if (qual in MU_QLTY_SNGL) and (qual2 in MU_QLTY_SNGL):
				if ptVtx >= 18. and ptVtx2 >= 18. :
					eta1 = float(Emu_br.tfMuonHwEta[i])*0.010875
					eta2 = float(Emu_br.tfMuonHwEta[j])*0.010875
					if abs(eta1) <= 2.1 and abs(eta2) <= 2.1	: L1_DoubleMu18_er2p1_flag = True

	    if (L1_SingleMu7_flag): L1_SingleMu7 +=1
	    if (L1_SingleMu22_flag): L1_SingleMu22 +=1
	    if (L1_SingleMu25_flag): L1_SingleMu25 +=1
	    if (L1_DoubleMu15_7_flag): L1_DoubleMu15_7 +=1
	    if (L1_DoubleMu18_er2p1_flag): L1_DoubleMu18_er2p1 +=1

	    if (L1_DoubleMu15_7_BMTF_flag): L1_DoubleMu15_7_BMTF +=1
	    if (L1_DoubleMu15_7_BMTF_UPT_flag): L1_DoubleMu15_7_BMTF_UPT +=1


            ######################################
            ###  Extra info from Kalman muons  ###
            ######################################
            #Alberto, is this deprecated in CMSW11?
            #for i in range(nKmtMu):
            #    BX      = int(Kmt_br.bx[i])
            #    qual    = int(Kmt_br.quality[i])
            #    ptVtx   = -1
            #    ptDispl = float(Kmt_br.ptUnconstrained[i])  ## Is there an offset by 1 for displaced muons? - AWB 2019.05.29
            #    eta     = float(Kmt_br.coarseEta[i])*0.010875
            #    chi2    = float(Kmt_br.approxChi2[i])
            #    if VERBOSE: print 'Internal muon %d BX = %d, qual = %d, ptVtx = %.1f, ptDispl = %.1f, eta = %.2f' % (i, BX, qual, ptVtx, ptDispl, eta)
                
            #    if (BX  !=  0): continue
                # if (qual < 12): continue  ## Quality assignment not the same as uGMT quality

            #    h_pt_displ_kmt.Fill( min( max( ptDispl, pt_bins[1]+0.01), pt_bins[2]-0.01) )
            #    h_chi2_kmt    .Fill( min( max( chi2,   chi_bins[1]+0.01), chi_bins[2]-0.01) )

        ## End loop: for jEvt in range(chains['Unp'][iCh].GetEntries()):
    ## End loop: for iCh in range(len(chains['Unp'])):

    print '\nFinished loop over chains'

    out_file.cd()

    h_dimu_rate_ptVtx.Scale(scale[evtclassid]/iEvt)
    h_dimu_rate_ptDisp.Scale(scale[evtclassid]/iEvt)
    h_dimu_rate_ptOr.Scale(scale[evtclassid]/iEvt)
    h_dimu_rate_ptAnd.Scale(scale[evtclassid]/iEvt)
    h_dimu_rate_ptUniqueDisp.Scale(scale[evtclassid]/iEvt)
    h_dimu_rate_ptUniqueVtx.Scale(scale[evtclassid]/iEvt)

    h_dimu_rate_ptVtx.SetLineColor(R.kBlack)
    h_dimu_rate_ptVtx.Write()

    h_dimu_rate_ptDisp.SetLineColor(R.kRed)
    h_dimu_rate_ptDisp.Write()

    h_dimu_rate_ptOr.SetLineColor(R.kBlue)
    h_dimu_rate_ptOr.Write()

    h_dimu_rate_ptAnd.SetLineColor(6)
    h_dimu_rate_ptAnd.Write()

    h_dimu_rate_ptUniqueDisp.SetLineColor(8)
    h_dimu_rate_ptUniqueDisp.Write()

    h_dimu_rate_ptUniqueVtx.SetLineColor(9)
    h_dimu_rate_ptUniqueVtx.Write()

#########################################################################

    L1_SingleMu22_BMTF = h_kBMTF_singlemu_rate_ptVtx.GetBinContent(23)
    L1_SingleMu22_BMTF_UPT = h_kBMTF_singlemu_rate_ptDisp.GetBinContent(23)

    print '\nWrote out file: plots/'+out_file_str+'.root'
    print "Events run over: ", iEvt

    print "L1_SingleMu7 : \t", L1_SingleMu7/200, "\t rate: ", (L1_SingleMu7/200)*scale[evtclassid]/iEvt
    print "L1_SingleMu22 : \t", L1_SingleMu22, "\t rate: ", L1_SingleMu22*scale[evtclassid]/iEvt
    print "L1_SingleMu25 : \t", L1_SingleMu25, "\t rate: ", L1_SingleMu25*scale[evtclassid]/iEvt
    print "L1_DoubleMu15_7 : \t", L1_DoubleMu15_7, "\t rate: ", L1_DoubleMu15_7*scale[evtclassid]/iEvt
    print "L1_DoubleMu18_er2p1 : \t", L1_DoubleMu18_er2p1, "\t rate: ", L1_DoubleMu18_er2p1*scale[evtclassid]/iEvt

    print "L1_DoubleMu15_7_BMTF : \t", L1_DoubleMu15_7_BMTF, "\t rate: ", L1_DoubleMu15_7_BMTF*scale[evtclassid]/iEvt
    print "L1_DoubleMu15_7_BMTF_UPT : \t", L1_DoubleMu15_7_BMTF_UPT, "\t rate: ", L1_DoubleMu15_7_BMTF_UPT*scale[evtclassid]/iEvt

    print "L1_SingleMu22_BMTF : \t", L1_SingleMu22_BMTF, "\t rate: ", L1_SingleMu22_BMTF*scale[evtclassid]/iEvt
    print "L1_SingleMu22_BMTF_UPT : \t", L1_SingleMu22_BMTF_UPT, "\t rate: ", L1_SingleMu22_BMTF_UPT*scale[evtclassid]/iEvt

############################### Single Muon #############################

    h_kBMTF_singlemu_rate_ptVtx.Scale(scale[evtclassid]/iEvt)
    h_kBMTF_singlemu_rate_ptDisp.Scale(scale[evtclassid]/iEvt)
    h_kBMTF_singlemu_rate_ptOr.Scale(scale[evtclassid]/iEvt)
    h_kBMTF_singlemu_rate_ptDisp_Lxy2.Scale(scale[evtclassid]/iEvt)

    h_kBMTF_singlemu_rate_ptVtx.SetLineColor(R.kBlack)
    h_kBMTF_singlemu_rate_ptVtx.Write()

    h_kBMTF_singlemu_rate_ptDisp.SetLineColor(R.kBlue)
    h_kBMTF_singlemu_rate_ptDisp.Write()

    h_kBMTF_singlemu_rate_ptOr.SetLineColor(R.kRed)
    h_kBMTF_singlemu_rate_ptOr.Write()

    h_kBMTF_singlemu_rate_ptDisp_Lxy2.SetLineColor(8)
    h_kBMTF_singlemu_rate_ptDisp_Lxy2.Write()

#########################################################################

    out_file.Close()
    del chains

if __name__ == '__main__':
    main()
