#! /usr/bin/env python

## **************************************************************** ##
##  Look at properties of displaced muons from Kalman algo in BMTF  ##
## **************************************************************** ##

import os, sys
import commands
import math

import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

PRT_EVT  =  1000   ## Print every Nth event
MAX_EVT  = 30000  ## Number of events to process
VERBOSE  = True  ## Verbose print-out

def main(inputTag):

    print '\nInside DisplacedMuons\n'
    sampleTag = "_HTo2XTo4Mu_"+inputTag
    
    #Please add if you what to process new files add the DPM folder where thet are contained inside the ConfigureInputFiles
    in_file_names = ConfigureInputFiles(inputTag)

    if not os.path.exists('plots'): os.makedirs('plots')

    out_file_str = 'DisplacedMuons'
    out_file_str += ('_%dk' % (MAX_EVT / 1000))
#    out_file = R.TFile('plots/'+out_file_str+sampleTag+'.root','recreate'
    tag = "_Dec2019"
    out_file = R.TFile('plots/'+out_file_str+sampleTag+tag+'.root','recreate')
    
    chains = {}
    chains['Evt'] = []  ## Event info
    chains['Gen'] = []  ## Gen info
    chains['Unp'] = []  ## Unpacked legacy BMTF
    chains['Emu'] = []  ## Emulated Kalman BMTF
    for i in range(len(in_file_names)):
        print 'Adding file %s' % in_file_names[i]
        chains['Evt'].append( R.TChain('l1EventTree/L1EventTree') )
        chains['Gen'].append( R.TChain('l1GeneratorTree/L1GenTree') )
        chains['Unp'].append( R.TChain('l1UpgradeTfMuonTree/L1UpgradeTfMuonTree') )
        chains['Emu'].append( R.TChain('l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree') )
        chains['Evt'][i].Add( in_file_names[i] )
        chains['Gen'][i].Add( in_file_names[i] )
        chains['Unp'][i].Add( in_file_names[i] )
        chains['Emu'][i].Add( in_file_names[i] )


    ###################
    ### Book histograms
    ###################

    pt_bins  = [200, 0, 200]
    dxy_bins = [500, 0, 500]
    chi_bins = [200, 0, 200]
    dxy_bins = [500, 0, 500]
    dxy_emu_bins = [3, 0, 3]
#    dxy_gen_bins = [250, 0, 500]
#    lxy_gen_bins = [250, 0, 500]
    dxy_gen_bins = [500, 0, 500]
    lxy_gen_bins = [500, 0, 500]
    nmu_gen_bins = [5, 0, 5]

    res_bins = [100, -1.0, 5.0]

    dif_bins = [100, -100.0, 100.0]

    h_pt_vtx_unp   = R.TH1F('h_pt_vtx_unp'+sampleTag,   'Legacy BMTF vertex-constrained pT spectrum;p_{T};Events',         pt_bins[0], pt_bins[1], pt_bins[2])
    h_pt_vtx_emu   = R.TH1F('h_pt_vtx_emu'+sampleTag,   'Kalman BMTF vertex-constrained pT spectrum;p_{T};Events',         pt_bins[0], pt_bins[1], pt_bins[2])
    h_pt_displ_emu = R.TH1F('h_pt_displ_emu'+sampleTag, 'Kalman BMTF non-vertex-constrained pT spectrum;p_{T};Events',     pt_bins[0], pt_bins[1], pt_bins[2])
    h_pt_displ_kmt = R.TH1F('h_pt_displ_kmt'+sampleTag, 'Internal Kalman non-vertex-constrained pT spectrum;p_{T};Events', pt_bins[0], pt_bins[1], pt_bins[2])

    h_chi2_kmt = R.TH1F('h_chi2_kmt'+sampleTag, 'internal Kalman BMTF track #chi^{2} distribution; #chi^{2};Events', chi_bins[0], chi_bins[1], chi_bins[2])
    h_dxy_emu  = R.TH1F('h_dxy_emu'+sampleTag, 'Kalman BMTF track dxy distribution; d_{xy}[cm]; Events', dxy_emu_bins[0], dxy_emu_bins[1], dxy_emu_bins[2])
    h_dxy_kmt = R.TH1F('h_dxy_kmt'+sampleTag, 'internal Kalman BMTF track dxy distribution; |d_{xy}|[cm]; Events', dxy_bins[0], dxy_bins[1], dxy_bins[2])

    h_pt_gen       = R.TH1F('h_pt_gen'+sampleTag,        'Gen level momementum if #eta < 0.9;p_{T};Events'                , pt_bins[0], pt_bins[1], pt_bins[2])
    h_dxy_gen      = R.TH1F('h_dxy_gen'+sampleTag,       'Gen level dxy if #eta < 0.9;d_{xy}[cm];Events'                , dxy_gen_bins[0], dxy_gen_bins[1], dxy_gen_bins[2])
    h_lxy_gen      = R.TH1F('h_lxy_gen'+sampleTag,       'Gen level Lxy if #eta < 0.9;L_{xy}[cm];Events'                , lxy_gen_bins[0], lxy_gen_bins[1], lxy_gen_bins[2])

    h_pt_gen_eventDisplaced       = R.TH1F('h_pt_gen_eventDisplaced'+sampleTag,        'Gen level momementum if #eta < 0.9, 30cm<l_{xy}<250cm;p_{T};Events',                     pt_bins[0], pt_bins[1], pt_bins[2])
    h_dxy_gen_eventDisplaced      = R.TH1F('h_dxy_gen_eventDisplaced'+sampleTag,       'Gen level d_{xy} if #eta < 0.9, 30cm<l_{xy}<250cm;d_{xy}[cm];Events',                    dxy_gen_bins[0], dxy_gen_bins[1], dxy_gen_bins[2])
    h_lxy_gen_eventDisplaced      = R.TH1F('h_lxy_gen_eventDisplaced'+sampleTag,       'Gen level L_{xy} if #eta < 0.9, 30cm<l_{xy}<250cm;L_{xy}[cm];Events',                    lxy_gen_bins[0], lxy_gen_bins[1], lxy_gen_bins[2])
    h_nmu_gen_eventDisplaced      = R.TH1F('h_nmu_gen_eventDisplaced'+sampleTag,       'Number of gen if #eta < 0.9, 30cm<l_{xy}<250cm;N_{#mu};Events',                          nmu_gen_bins[0], nmu_gen_bins[1], nmu_gen_bins[2])

    h_pt_vtx_unp_eventDisplaced   = R.TH1F('h_pt_vtx_unp_eventDisplaced'+sampleTag,   'Legacy BMTF vertex-constrained pT spectrum, 30cm<l_{xy}<250cm; p_{T};Events',         pt_bins[0], pt_bins[1], pt_bins[2])
    h_pt_vtx_emu_eventDisplaced   = R.TH1F('h_pt_vtx_emu_eventDisplaced'+sampleTag,   'Kalman BMTF vertex-constrained pT spectrum, 30cm<l_{xy}<250cm; p_{T};Events',         pt_bins[0], pt_bins[1], pt_bins[2])
    h_pt_displ_emu_eventDisplaced = R.TH1F('h_pt_displ_emu_eventDisplaced'+sampleTag, 'Kalman BMTF non-vertex-constrained pT spectrum, 30cm<l_{xy}<250cm;p_{T};Events',     pt_bins[0], pt_bins[1], pt_bins[2])
    h_pt_kbmtf_max_eventDisplaced = R.TH1F('h_pt_kbmtf_max_eventDisplaced'+sampleTag, 'max Kalman BMTF (vertex-constrained, non-vertex-constrained) spectrum pT, 30cm<l_{xy}<250cm;p_{T};Events',     pt_bins[0], pt_bins[1], pt_bins[2])

    h_res_vtx_unp_eventDisplaced   = R.TH1F('h_res_vtx_unp_eventDisplaced'+sampleTag,   'Legacy BMTF ; (p_{T}-p^{gen}_{T})/p^{gen}_{T};Events',                            res_bins[0], res_bins[1], res_bins[2])
    h_res_vtx_emu_eventDisplaced   = R.TH1F('h_res_vtx_emu_eventDisplaced'+sampleTag,   'Kalman BMTF vertex-constrained ; (p_{T}-p^{gen}_{T})/p^{gen}_{T};Events',         res_bins[0], res_bins[1], res_bins[2])
    h_res_displ_emu_eventDisplaced = R.TH1F('h_res_displ_emu_eventDisplaced'+sampleTag, 'Kalman BMTF non-vertex-constrained ; (p_{T}-p^{gen}_{T})/p^{gen}_{T});Events',    res_bins[0], res_bins[1], res_bins[2])
    h_res_displ_kmt_eventDisplaced = R.TH1F('h_res_displ_kmt_eventDisplaced'+sampleTag, 'internal Kalman BMTF non-vertex-constrained ; (p_{T}-p^{gen}_{T})/p^{gen}_{T});Events',    res_bins[0], res_bins[1], res_bins[2])
    h_res_kbmtf_max_eventDisplaced = R.TH1F('h_res_kbmtf_max_eventDisplaced'+sampleTag, 'Kalman BMTF max ; (p_{T}-p^{gen}_{T})/p^{gen}_{T};Events',                        res_bins[0], res_bins[1], res_bins[2])

    h_resc_kbmtf_eventDisplaced = R.TH1F('h_resc_kbmtf_eventDisplaced'+sampleTag, 'BMTF vs kBMTF vertex constrained ; (p^{kBMTF(disp)}_{T}-p^{kBMTF}_{T})/p^{gen}_{T};Events',       res_bins[0], res_bins[1], res_bins[2])
    h_resdxy_displ_emu_eventDisplaced = R.TH1F('h_resdxy_displ_emu_eventDisplaced'+sampleTag, 'internal Kalman BMTF non-vertex-constrained ; (d_{xy}-d^{gen}_{xy})/d^{gen}_{xy};Events',    res_bins[0], res_bins[1], res_bins[2])

    #TOBECOMPLETED.
    h_difdxy_displ_emu_0         = R.TH1F('h_difdxy_displ_emu_0'+sampleTag, 'internal Kalman BMTF non-vertex-constrained ; (d_{xy}-d^{gen}_{xy});Events',    dif_bins[0], dif_bins[1], dif_bins[2])
    h_difdxy_displ_emu_1         = R.TH1F('h_difdxy_displ_emu_1'+sampleTag, 'internal Kalman BMTF non-vertex-constrained ; (d_{xy}-d^{gen}_{xy});Events',    dif_bins[0], dif_bins[1], dif_bins[2])
    h_difdxy_displ_emu_2         = R.TH1F('h_difdxy_displ_emu_2'+sampleTag, 'internal Kalman BMTF non-vertex-constrained ; (d_{xy}-d^{gen}_{xy});Events',    dif_bins[0], dif_bins[1], dif_bins[2])
    h_difdxy_displ_emu_3         = R.TH1F('h_difdxy_displ_emu_3'+sampleTag, 'internal Kalman BMTF non-vertex-constrained ; (d_{xy}-d^{gen}_{xy});Events',    dif_bins[0], dif_bins[1], dif_bins[2])
    h_difdxy_displ_emu_4         = R.TH1F('h_difdxy_displ_emu_4'+sampleTag, 'internal Kalman BMTF non-vertex-constrained ; (d_{xy}-d^{gen}_{xy});Events',    dif_bins[0], dif_bins[1], dif_bins[2])

    iEvt = -1
    print '\nEntering loop over chains'
    for iCh in range(len(chains['Unp'])):

        if iEvt > MAX_EVT: break

        ## Faster tecnhique, inspired by https://github.com/thomreis/l1tMuonTools/blob/master/L1Analysis.py
        Evt_br = R.L1Analysis.L1AnalysisEventDataFormat()
        Gen_br = R.L1Analysis.L1AnalysisGeneratorDataFormat()
        Unp_br = R.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
        Emu_br = R.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
        Kmt_br = R.L1Analysis.L1AnalysisBMTFOutputDataFormat()

        chains['Evt'][iCh].SetBranchAddress('Event',               R.AddressOf(Evt_br))
        chains['Gen'][iCh].SetBranchAddress('Generator',           R.AddressOf(Gen_br))
        chains['Unp'][iCh].SetBranchAddress('L1UpgradeBmtfMuon',   R.AddressOf(Unp_br))
        chains['Emu'][iCh].SetBranchAddress('L1UpgradeBmtfMuon',   R.AddressOf(Emu_br))
        chains['Emu'][iCh].SetBranchAddress('L1UpgradeBmtfOutput', R.AddressOf(Kmt_br))


        print '\nEntering loop over events for chain %d' % iCh
        for jEvt in range(chains['Unp'][iCh].GetEntries()):

            if iEvt > MAX_EVT: break
            iEvt += 1
            if iEvt % PRT_EVT is 0: print '\nEvent # %d (%dth in chain)' % (iEvt, jEvt)

            chains['Evt'][iCh].GetEntry(jEvt)
            chains['Gen'][iCh].GetEntry(jEvt)
            chains['Unp'][iCh].GetEntry(jEvt)
            chains['Emu'][iCh].GetEntry(jEvt)

            # ## Use these lines if you don't explicitly define the DataFormat and then do SetBranchAddress above
            # Evt_br = chains['Evt'][iCh].Event
            # Unp_br = chains['Unp'][iCh].L1UpgradeBmtfMuon
            # Emu_br = chains['Emu'][iCh].L1UpgradeBmtfMuon
            
            if iEvt % PRT_EVT is 0: print '  * Run %d, LS %d, event %d' % (int(Evt_br.run), int(Evt_br.lumi), int(Evt_br.event))

            nUnpMu = int(Unp_br.nTfMuons)
            nEmuMu = int(Emu_br.nTfMuons)
            nKmtMu = int(Kmt_br.nTrks)
            
            #Count gen level muons:
            nGenMu = 0
            for kGenParId in Gen_br.partId:
                if abs(kGenParId) == 13: nGenMu = nGenMu+1
            if VERBOSE: print 'nGen muons %i' % nGenMu

            if (VERBOSE == True and (nUnpMu > 0 or nEmuMu > 0)):
                print 'Unpacked = %d, emulated = %d (internal %d) total muons in collection' % (nUnpMu, nEmuMu, nKmtMu)
        
            #################################
            ###  Generator level  muons   ###
            #################################
            eventLxy = []
            eventLz  = []
            eventDxy = []
            eventDz  = []
            eventPt  = []
            eventEta = []

            eventDisplacedPt = []

            isDisplaced = True
            for i, pdgId in enumerate(Gen_br.partId):
                pt      = float(Gen_br.partPt[i])     ## Gen level pt momentum [GeV]
                eta     = float(Gen_br.partEta[i])    ## Gen level eta
                phi     = float(Gen_br.partPhi[i])    ## Gen level phi
                parent  = int(Gen_br.partParent[i])   ## Gen level parent   [GeV]
                x       = float(Gen_br.partVx[i])   ## Gen level X coordinate
                y       = float(Gen_br.partVy[i])   ## Gen level Y coordinate
                z       = float(Gen_br.partVz[i])   ## Gen level Z coordinate
                if (pdgId == 6000111 or pdgId == 6000113): # particle ID is a long lived X
                    pv_x = float(Gen_br.partVx[i])
                    pv_y = float(Gen_br.partVx[i])
                    pv_z = float(Gen_br.partVx[i])

                if abs(pdgId) != 13: continue
                if abs(parent) != 6000113: continue                                  ## pdgID for X in H->XX model
                px      = pt * math.cos(phi)                                         ## mu Gen level px momentum [GeV]
                py      = pt * math.sin(phi)                                         ## mu Gen level py momentum [GeV]
                pz      = pt * math.sinh(eta)                                        ## mu Gen level pz momentum [GeV]
                lxy = abs(math.sqrt((x - pv_x)*(x - pv_x)+(y - pv_y)*(y - pv_y)))    ## Transverse decay length
                lz  = abs(math.sqrt((z - pv_z)*(z - pv_z)))                          ## Longitudinal decay length
                dxy = abs((x*py+y*px)/pt)
                dz  = abs(z-pz/pt*(x*px+y*py)/pt)
                if VERBOSE: print 'gen muon %d pt = %.2f, eta = %.2f, (x:%.2f, y:%.2f , z:%.2f), lxy = %.2f, lz = %.2f, dxy = %.2f, dz = %.2f, parent = %i' % (i, pt, eta, x, y, z, lxy, lz, dxy, dz, parent)
                if abs(eta) < 0.9:                 
                    h_pt_gen.Fill( pt )
                    h_lxy_gen.Fill( lxy )
                    h_dxy_gen.Fill( dxy )
                    eventDxy.append(dxy)
                    eventDz.append(dz)
                    eventLxy.append(lxy)
                    eventLz.append(lz)
                    eventPt.append(pt)
                    eventEta.append(eta)
                    eventDisplacedPt.append(pt)

            #displaced category?
            if len(eventLxy) >0:
                for kEventLxy, kEventLz, kEventDxy in zip(eventLxy, eventLz, eventDxy):
                    if kEventLxy < 30 or kEventLxy > 250 or kEventLz > 500 or kEventDxy < 30:
                        isDisplaced = False # both genLxy should be larger than 30 but smaller than 250 cm than something.
            else:
                isDisplaced = False

            if isDisplaced == True:
                #print len(eventDisplacedPt), len(eventLxy)
                h_nmu_gen_eventDisplaced.Fill( len(eventDisplacedPt) )
                for i in range(len(eventDisplacedPt)):
                    h_pt_gen_eventDisplaced.Fill( eventDisplacedPt[i] )
                    h_dxy_gen_eventDisplaced.Fill( eventDxy[i] )
                    h_lxy_gen_eventDisplaced.Fill( eventLxy[i] )

            if VERBOSE and isDisplaced == True: 
                print ' displaced event lxy, lz: ', eventLxy, eventLz

            #################################
            ###  Unpacked (legacy) muons  ###
            #################################
            for i in range(nUnpMu):
                BX      = int(Unp_br.tfMuonBx[i])
                qual    = int(Unp_br.tfMuonHwQual[i])
                ptVtx   = float(Unp_br.tfMuonHwPt[i] - 1)*0.5   ## Vertex-constrained (standard) pT is stored in 0.5 GeV steps
                ptDispl = float(Unp_br.tfMuonHwPtDispl[i] - 1)  ## Is there an offset by 1 for displaced muons? - AWB 2019.05.29
                eta     = float(Unp_br.tfMuonHwEta[i])*0.010875
                if VERBOSE: print 'Unpacked muon %d BX = %d, qual = %d, ptVtx = %.1f, ptDispl = %.1f, eta = %.2f' % (i, BX, qual, ptVtx, ptDispl, eta)
                if VERBOSE and ptVtx > 138: print ' ------- > Unpacked muon %d BX = %d, qual = %d, ptVtx = %.1f, ptDispl = %.1f, eta = %.2f' % (i, BX, qual, ptVtx, ptDispl, eta) #overflow 139.5
                
                if (BX  !=  0): continue
                if (qual < 12): continue

                h_pt_vtx_unp.Fill( min( max( ptVtx+0.01, pt_bins[1]+0.01), pt_bins[2]-0.01) )
                if isDisplaced == True: 
                    h_pt_vtx_unp_eventDisplaced.Fill( min( max( ptVtx+0.01, pt_bins[1]+0.01), pt_bins[2]-0.01) )
                    if len(eventDisplacedPt) == 1 and ptVtx > 0:              # there is only one displaced gen muon in barrel
                        res_bmtf =  (ptVtx - eventDisplacedPt[0])/eventDisplacedPt[0]
                        h_res_vtx_unp_eventDisplaced.Fill( min(res_bmtf, res_bins[2]-0.01 ))
                    
            #################################
            ###  Emulated (Kalman) muons  ###
            #################################
            for i in range(nEmuMu):
                BX      = int(Emu_br.tfMuonBx[i])
                qual    = int(Emu_br.tfMuonHwQual[i])
                ptVtx   = float(Emu_br.tfMuonHwPt[i] - 1)*0.5  ## Vertex-constrained (standard) pT is stored in 0.5 GeV steps
                ptDispl = float(Emu_br.tfMuonHwPtDispl[i])     ## Is there an offset by 1 for displaced muons? - AWB 2019.05.29
                eta     = float(Emu_br.tfMuonHwEta[i])*0.010875
                dxy     = float(Emu_br.tfMuonHwDXY[i])         
                if VERBOSE: print 'Emulated muon %d BX = %d, qual = %d, ptVtx = %.1f, ptDispl = %.1f, eta = %.2f, pTRaw = %.2f, pTDispRaw = %.2f ' % (i, BX, qual, ptVtx, ptDispl, eta, Emu_br.tfMuonHwPt[i], Emu_br.tfMuonHwPtDispl[i])
                
                if (BX  !=  0): continue
                if (qual < 12): continue

                h_pt_vtx_emu  .Fill( min( max( ptVtx+0.01,   pt_bins[1]+0.01), pt_bins[2]-0.01) )
                h_dxy_emu  .Fill( dxy )
                h_pt_displ_emu.Fill( min( max( ptDispl+0.01, pt_bins[1]+0.01), pt_bins[2]-0.01) )

                if isDisplaced == True:
                    h_pt_vtx_emu_eventDisplaced  .Fill( min( max( ptVtx+0.01,   pt_bins[1]+0.01), pt_bins[2]-0.01) )
                    h_pt_displ_emu_eventDisplaced.Fill( min( max( ptDispl+0.01, pt_bins[1]+0.01), pt_bins[2]-0.01) )
                    pt_constrained = min( max( ptVtx+0.01,   pt_bins[1]+0.01), pt_bins[2]-0.01)
                    pt_unconstrained = min( max( ptDispl+0.01, pt_bins[1]+0.01), pt_bins[2]-0.01)
                    h_pt_kbmtf_max_eventDisplaced.Fill( max(pt_constrained,pt_unconstrained) )
                    if len(eventDisplacedPt) == 1:  # there is only one displaced gen muon in barrel 
                        if ptVtx > 0:             
                            res_kbmtf =  (ptVtx - eventDisplacedPt[0])/eventDisplacedPt[0]             # resolution kbmtf vertex constrained
                            h_res_vtx_emu_eventDisplaced.Fill( min(res_kbmtf, res_bins[2]-0.01) ) 
                        if ptDispl > 0:
                            res_kbmtf_disp =  (ptDispl - eventDisplacedPt[0])/eventDisplacedPt[0]      # resolution kbmtf vertex un-constrained
                            h_res_displ_emu_eventDisplaced.Fill( min(res_kbmtf_disp, res_bins[2]-0.01) )
                                
                        if ptVtx > 0 or ptDispl > 0:
                            res_kbmtf_max = (max(ptVtx, ptDispl) - eventDisplacedPt[0])/eventDisplacedPt[0]    # resolution kbmtf maximal
                            h_res_kbmtf_max_eventDisplaced.Fill( min(res_kbmtf_max, res_bins[2]-0.01))
                            
                            resc_kbmtf = (ptDispl - ptVtx)/eventDisplacedPt[0]                         #pT comparison between different kBMTF measurements.
                            h_resc_kbmtf_eventDisplaced.Fill( min(resc_kbmtf, res_bins[2]-0.01) )

                if len(eventPt) == 1:
                    if VERBOSE: print  '  o??o > ONLY one GEN MU'
                    for i in range(len(eventPt)):
                        if abs(eta-eventEta[i]) < 0.1:
                            if VERBOSE: print  '  oooo > MATCHED-HW mu genEta = %f, kEta = %f, genDxy = %.1f, kDxy = %.1f, diff %.2f' % ( eventEta[i], eta, eventDxy[i], dxy, abs(dxy) - eventDxy[i])

            ######################################
            ###  Extra info from Kalman muons  ###
            ######################################
            for i in range(nKmtMu):
                BX      = int(Kmt_br.bx[i])
                qual    = int(Kmt_br.quality[i])
                ptVtx   = -1
                ptDispl = float(Kmt_br.ptUnconstrained[i])  ## Is there an offset by 1 for displaced muons? - AWB 2019.05.29
                eta     = float(Kmt_br.coarseEta[i])*0.010875
                chi2    = float(Kmt_br.approxChi2[i])
                dxy     = float(Kmt_br.dxy[i])/10 ## Which units? I think they are in mm that I converted to cm here.
                if VERBOSE: print 'Internal muon %d BX = %d, qual = %d, ptVtx = %.1f, ptDispl = %.1f, eta = %.2f, dxy = %.2f' % (i, BX, qual, ptVtx, ptDispl, eta, dxy)
                if VERBOSE and ptDispl > 140: print ' ---- > Internal muon %d BX = %d, qual = %d, ptVtx = %.1f, ptDispl = %.1f, eta = %.2f, dxy = %.2f' % (i, BX, qual, ptVtx, ptDispl, eta, dxy) #overflow at 142.5
                
                if (BX  !=  0): continue
                # if (qual < 12): continue  ## Quality assignment not the same as uGMT quality

                h_pt_displ_kmt.Fill( min( max( ptDispl, pt_bins[1]+0.01), pt_bins[2]-0.01) )
                h_chi2_kmt    .Fill( min( max( chi2,   chi_bins[1]+0.01), chi_bins[2]-0.01) )
                h_dxy_kmt     .Fill( abs(dxy) )

                if len(eventPt) < 4:
                    for i in range(len(eventPt)):
                        if abs(eta-eventEta[i]) < 0.05:
                            if VERBOSE: print  '  oooo > MATCHED-INTERNAL mu genEta = %f, kEta = %f, genDxy = %.1f, kDxy = %.1f, diff %.2f, res = %.2f' % ( eventEta[i], eta, eventDxy[i], dxy, abs(dxy) - eventDxy[i], (abs(dxy) - eventDxy[i])/eventDxy[i])
                            
                        if eventDxy[i] < 30:
                            h_difdxy_displ_emu_0.Fill(abs(dxy) - eventDxy[i])
                        if 30 < eventDxy[i] < 60:
                            h_difdxy_displ_emu_1.Fill(abs(dxy) - eventDxy[i])
                        if 60 < eventDxy[i] < 90:
                            h_difdxy_displ_emu_2.Fill(abs(dxy) - eventDxy[i])
                        if 90 < eventDxy[i] < 300:
                            h_difdxy_displ_emu_3.Fill(abs(dxy) - eventDxy[i])
                        if 90 < eventDxy[i] < 300 and ptDispl > 10:
                            h_difdxy_displ_emu_4.Fill(abs(dxy) - eventDxy[i])

                if isDisplaced == True:
                    if len(eventDisplacedPt) == 1:                                                              # there is only one displaced gen muon in barrel 
                        res_kbmtf_disp =  (ptDispl - eventDisplacedPt[0])/eventDisplacedPt[0]                   # resolution kbmtf vertex un-constrained
                        h_res_displ_kmt_eventDisplaced.Fill(min(res_kbmtf_disp, res_bins[2]-0.01))

                        res_dxy_kbmtf_disp = (abs(dxy) - eventDxy[0])/eventDxy[0]                                # resolution kbmtf vertex un-constrained
                        h_resdxy_displ_emu_eventDisplaced.Fill(min(res_dxy_kbmtf_disp, res_bins[2]-0.01))
                        if VERBOSE : print ' xxxxx > matched muon ptGen = %.1f, ptDispl = %.1f, dxyGen = %.2f,  dxyDisp = %.2f' % (eventDisplacedPt[0], ptDispl, eventDxy[0], dxy) #overflow at 142.5

        ## End loop: for jEvt in range(chains['Unp'][iCh].GetEntries()):
    ## End loop: for iCh in range(len(chains['Unp'])):

    print '\nFinished loop over chains'

    out_file.cd()

    h_pt_vtx_unp.SetLineWidth(2)
    h_pt_vtx_unp.SetLineColor(R.kBlack)
    h_pt_vtx_unp.Write()

    h_pt_vtx_emu.SetLineWidth(2)
    h_pt_vtx_emu.SetLineColor(R.kBlue)
    h_pt_vtx_emu.Write()

    h_dxy_emu.SetLineWidth(2)
    h_dxy_emu.SetLineColor(R.kBlue)
    h_dxy_emu.Write()

    h_pt_displ_emu.SetLineWidth(2)
    h_pt_displ_emu.SetLineColor(R.kRed)
    h_pt_displ_emu.Write()
    
    h_pt_displ_kmt.SetLineWidth(2)
    h_pt_displ_kmt.SetLineColor(R.kMagenta)
    h_pt_displ_kmt.Write()
    
    h_pt_gen.SetLineWidth(2)
    h_pt_gen.SetLineColor(R.kGray+2)
    h_pt_gen.Write()

    h_dxy_gen.SetLineWidth(2)
    h_dxy_gen.SetLineColor(R.kGray+2)
    h_dxy_gen.Write()

    h_lxy_gen.SetLineWidth(2)
    h_lxy_gen.SetLineColor(R.kGray+2)
    h_lxy_gen.Write()

    h_chi2_kmt.SetLineWidth(2)
    h_chi2_kmt.SetLineColor(R.kBlack)
    h_chi2_kmt.Write()

    h_dxy_kmt.SetLineWidth(2)
    h_dxy_kmt.SetLineColor(R.kGray)
    h_dxy_kmt.Write()

    # displaced Histograms
    h_pt_vtx_unp_eventDisplaced.SetLineWidth(2)
    h_pt_vtx_unp_eventDisplaced.SetLineColor(R.kBlack)
    h_pt_vtx_unp_eventDisplaced.Write()

    h_pt_vtx_emu_eventDisplaced.SetLineWidth(2)
    h_pt_vtx_emu_eventDisplaced.SetLineColor(R.kBlue)
    h_pt_vtx_emu_eventDisplaced.Write()

    h_pt_displ_emu_eventDisplaced.SetLineWidth(2)
    h_pt_displ_emu_eventDisplaced.SetLineColor(R.kRed)
    h_pt_displ_emu_eventDisplaced.Write()

    h_pt_kbmtf_max_eventDisplaced.SetLineWidth(2)
    h_pt_kbmtf_max_eventDisplaced.SetLineColor(R.kRed)
    h_pt_kbmtf_max_eventDisplaced.Write()        

    h_pt_gen_eventDisplaced.SetLineWidth(2)
    h_pt_gen_eventDisplaced.SetLineColor(R.kGray+2)
    h_pt_gen_eventDisplaced.Write()

    h_nmu_gen_eventDisplaced.SetLineWidth(2)
    h_nmu_gen_eventDisplaced.SetLineColor(R.kGray+2)
    h_nmu_gen_eventDisplaced.Write()

    h_dxy_gen_eventDisplaced.SetLineWidth(2)
    h_dxy_gen_eventDisplaced.SetLineColor(R.kGray+2)
    h_dxy_gen_eventDisplaced.Write()

    h_lxy_gen_eventDisplaced.SetLineWidth(2)
    h_lxy_gen_eventDisplaced.SetLineColor(R.kGray+2)
    h_lxy_gen_eventDisplaced.Write()
    
    h_res_vtx_unp_eventDisplaced.SetLineWidth(2)
    h_res_vtx_unp_eventDisplaced.SetLineColor(R.kRed)
    h_res_vtx_unp_eventDisplaced.Write()

    h_res_vtx_emu_eventDisplaced.SetLineWidth(2)
    h_res_vtx_emu_eventDisplaced.SetLineColor(R.kRed)
    h_res_vtx_emu_eventDisplaced.Write()

    h_res_displ_emu_eventDisplaced.SetLineWidth(2)
    h_res_displ_emu_eventDisplaced.SetLineColor(R.kRed)
    h_res_displ_emu_eventDisplaced.Write()

    h_res_displ_kmt_eventDisplaced.SetLineWidth(2)
    h_res_displ_kmt_eventDisplaced.SetLineColor(R.kRed)
    h_res_displ_kmt_eventDisplaced.Write()

    h_res_kbmtf_max_eventDisplaced.SetLineWidth(2)
    h_res_kbmtf_max_eventDisplaced.SetLineColor(R.kRed)
    h_res_kbmtf_max_eventDisplaced.Write()

    h_resc_kbmtf_eventDisplaced.SetLineWidth(2)
    h_resc_kbmtf_eventDisplaced.SetLineColor(R.kRed)
    h_resc_kbmtf_eventDisplaced.Write()

    h_resdxy_displ_emu_eventDisplaced.SetLineWidth(2)
    h_resdxy_displ_emu_eventDisplaced.SetLineColor(R.kRed)
    h_resdxy_displ_emu_eventDisplaced.Write()

    h_difdxy_displ_emu_0.SetLineWidth(2)
    h_difdxy_displ_emu_0.SetLineColor(R.kRed)
    h_difdxy_displ_emu_0.Write()

    h_difdxy_displ_emu_1.SetLineWidth(2)
    h_difdxy_displ_emu_1.SetLineColor(R.kRed)
    h_difdxy_displ_emu_1.Write()

    h_difdxy_displ_emu_2.SetLineWidth(2)
    h_difdxy_displ_emu_2.SetLineColor(R.kRed)
    h_difdxy_displ_emu_2.Write()

    h_difdxy_displ_emu_3.SetLineWidth(2)
    h_difdxy_displ_emu_3.SetLineColor(R.kRed)
    h_difdxy_displ_emu_3.Write()

    h_difdxy_displ_emu_4.SetLineWidth(2)
    h_difdxy_displ_emu_4.SetLineColor(R.kRed)
    h_difdxy_displ_emu_4.Write()

    out_file.Close()
    del chains

    print '\nWrote out file: plots/'+out_file_str+sampleTag+tag+'.root'


def ConfigureInputFiles(inputTag):
    # Simple function that finds the location of a given set of NTuples that will be process by the analyzer.
    mySE = 'srm://hephyse.oeaw.ac.at:8446/srm/ymanagerv2?SFN=/dpm/oeaw.ac.at/home/cms/store/user/escalant/'
    preFix = 'root://cmsxrootd.fnal.gov//store/user/escalant/'


    whatToProcessinputFiles = {
        #'200_50_200':'/HTo2LongLivedTo2mu2jets_MH-200_MFF-50_CTau-200mm_TuneCP5_13TeV_pythia8/crab_HTo2LongLivedTo2mu2jets_MH-200_MFF-50_CTau-200mm_TuneCP5_13TeV_pythia8_testL1Ntuples/190606_195230/0000/', 
        #'200_50_20'  :'/HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-20mm_TuneCP5_13TeV_pythia8/crab_HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-20mm_TuneCP5_13TeV_pythia8_testL1Ntuples_withDisplacedGen_v1/190607_112212/0000/', 
        #'200_50_200' :'/HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-200mm_TuneCP5_13TeV_pythia8/crab_HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-200mm_TuneCP5_13TeV_pythia8_testL1Ntuples_withDisplacedGen_v1/190607_112059/0000/', 
        #'200_50_2000':'/HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-2000mm_TuneCP5_13TeV_pythia8/crab_HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-2000mm_TuneCP5_13TeV_pythia8_testL1Ntuples_withDisplacedGen_v1/190607_111945/0000/', 
        '200_50_20'  :'/HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-20mm_TuneCP5_13TeV_pythia8/crab_HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-20mm_TuneCP5_13TeV_pythia8_testL1Ntuples_withDisplacedGen_v3/190610_212836/0000/', 
        '200_50_200' :'/HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-200mm_TuneCP5_13TeV_pythia8/crab_HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-200mm_TuneCP5_13TeV_pythia8_testL1Ntuples_withDisplacedGen_v3/190610_212720/0000/', 
        '200_50_2000':'/HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-2000mm_TuneCP5_13TeV_pythia8/crab_HTo2LongLivedTo4mu_MH-200_MFF-50_CTau-2000mm_TuneCP5_13TeV_pythia8_testL1Ntuples_withDisplacedGen_v3/190610_212601/0000/', 
        }



    filesToProcess = list(commands.getstatusoutput('gfal-ls '+mySE+whatToProcessinputFiles[inputTag]))[1].split()
    for i in range(len(filesToProcess)):
        filesToProcess[i] = preFix + whatToProcessinputFiles[inputTag]+filesToProcess[i]
    return filesToProcess

if __name__ == '__main__':

    if len(sys.argv) == 2:
        print sys.argv[1], sys.argv[0]
        PlotTags = sys.argv[1]

    main(PlotTags)
    #main('200_50_2000')
    #main('200_50_200')
    #main('200_50_20')
