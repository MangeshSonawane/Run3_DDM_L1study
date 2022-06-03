#! /usr/bin/env python

## **************************************************************** ##
##  Look at properties of displaced muons from Kalman algo in BMTF  ##
## **************************************************************** ##

import os, sys
from math import *
from array import array

import ROOT as R
R.gROOT.SetBatch(False)  ## Don't print histograms to screen while processing

PRT_EVT  = 10000	 	## Print every Nth event
MAX_EVT  = 300000	## Number of events to process
VERBOSE  = False	## Verbose print-out

scale = [1., 1., 2544*11246., 1.]

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

	r = (490. - Lxy)
	z = r*sinh(eta)

	zStar = vz + z
	xStar = vx + r*cos(phi)
	yStar = vy + r*sin(phi)
 
	rStar = sqrt(xStar**2 + yStar**2)
	
	GenEta = asinh(zStar/490.)
	
	if (xStar >= 0) : GenPhi = atan(yStar/xStar)
      	elif (yStar >= 0 and xStar < 0) : GenPhi = pi + atan(yStar/xStar)
      	elif (yStar <= 0 and xStar < 0) : GenPhi = atan(yStar/xStar) - pi 	
	
	return GenEta, GenPhi

def IsinacceptBarrel(Gen_br, i):

	vx = float(Gen_br.partVx[i])
	vy = float(Gen_br.partVy[i])
	vz = float(Gen_br.partVz[i])

	eta = float(Gen_br.partEta[i])
	
	if abs(vz) >= 650. : return False
	
	Lxy = sqrt(vx**2 + vy**2)
	if Lxy > 490. : return False

	maxeta = -log(tan(0.5*atan((700.-Lxy)/(650.-vz))))
	mineta = -log(tan (0.5*(pi-atan((700.-Lxy)/(650.+vz)))))

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

	if abs(eta) <= 1.2 : return False
	if abs(eta) >= 2.5 : return False	

	return True

def getDxy(vx, vy, phi):
	return abs(vx*sin(phi)-vy*cos(phi))


##################################################################################################################


def main():

    print '\nInside DisplacedMuons\n'
    evtclass = ["signal_1500", "signal_3000","NuGun","MuGun"]
    evtclassid = 3
    inputdir = ['/eos/user/s/sonawane/temp/L1Ntuples/signal_1500_tuples/',
		'/eos/user/s/sonawane/temp/L1Ntuples/signal_3000_tuples/ntuples_15_06_21/HTo2LongLivedTo4mu_MH-125_MFF-50_CTau-3000mm_11_2_X_1623671915/',
		'/eos/user/s/sonawane/temp/L1Ntuples/NuGunmod_jobs_12k/ntuples/',
		'/eos/user/s/sonawane/temp/L1Ntuples/Displaced_mu_gun_tuples/ntuples_16_06_21/']

    workdir = '/afs/cern.ch/user/s/sonawane/L1T/L1studies/L1_scripts_Alberto/L1RunIII/macros/'
    in_file_names = []

    if evtclassid == 3:
	for s in ['DisplacedMuGun_Pt2to10_11_2_X_1623847533/','DisplacedMuGun_Pt10to30_11_2_X_1623847663/', 'DisplacedMuGun_Pt30to100_11_2_X_1623847721/']:
		ntupledir = inputdir[3]+s
		for i in range(23) :
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
    chains['uGT'] = []  ## Global Trigger
    chains['Gen'] = []  ## Generator information
    for i in range(len(in_file_names)):
        print 'Adding file %s' % in_file_names[i]
        chains['Evt'].append( R.TChain('l1EventTree/L1EventTree') )
        chains['Unp'].append( R.TChain('l1UpgradeTfMuonTree/L1UpgradeTfMuonTree') )
        chains['Emu'].append( R.TChain('l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree') )
        chains['uGT'].append( R.TChain('l1UpgradeEmuTree/L1UpgradeTree') )
        chains['Gen'].append( R.TChain('l1GeneratorTree/L1GenTree') )
        chains['Evt'][i].Add( in_file_names[i] )
        chains['Unp'][i].Add( in_file_names[i] )
        chains['Emu'][i].Add( in_file_names[i] )
        chains['uGT'][i].Add( in_file_names[i] )
        chains['Gen'][i].Add( in_file_names[i] )


    ###################
    ### Book histograms
    ###################

    pt_bins  = [150, 0, 300]
    q_pt_bins  = [100, -0.2, 0.2]
#    chi_bins = [100, 0, 100]
    dxy_bins = [20,   0, 100]
    vxyz_bins = [200,   -200, 200]
    lxy_bins = [30,   0, 150]
    phi_bins = [80, -4.,4.]
    eta_bins = [30, -3., 3.]
    dphi_bins = [40, 0., 4.]
    deta_bins = [30, 0., 3.]
#    Vx_bins  = [100,-50, 50]
#    phi_bins = [80, -4, 4]
    EC_eta_bins = [-3., -2.5, -2.1, -1.6, -1.2, 1.2, 1.6, 2.1, 2.5, 3.]

    h_dxy_blank = R.TH1F('h_dxy_blank', 'Blank Gen Dxy for TEfficiency plots', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_pt_blank = R.TH1F('h_pt_blank', 'Blank Gen Pt for TEfficiency plot', pt_bins[0], pt_bins[1], pt_bins[2])
    h_eta_blank = R.TH1F('h_eta_blank', 'Blank Endcap Gen Eta for TEfficiency plot', len(EC_eta_bins)-1, array('d',EC_eta_bins))

    h_eff_gen_dxy   = R.TH1F('h_eff_gen_dxy',   'Eff numerator for Gen Dxy',	dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_eff_gen_dxy_pt10   = R.TH1F('h_eff_gen_dxy_pt10',   'Eff numerator for Gen Dxy, pt > 10 GeV',	dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_eff_gen_dxy_pt20   = R.TH1F('h_eff_gen_dxy_pt20',   'Eff numerator for Gen Dxy, pt > 20 GeV',	dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_eff_gen_dxy_pt30   = R.TH1F('h_eff_gen_dxy_pt30',   'Eff numerator for Gen Dxy, pt > 30 GeV',	dxy_bins[0], dxy_bins[1], dxy_bins[2])

    h_gen_dxy   = R.TH1F('h_gen_dxy',   'Denominator Gen Dxy',	dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_gen_dxy_pt10 = R.TH1F('h_gen_dxy_pt_10',   'Denominator Gen Dxy for pt Thresh 10',	dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_gen_dxy_pt20 = R.TH1F('h_gen_dxy_pt_20',   'Denominator Gen Dxy for pt Thresh 20',	dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_gen_dxy_pt30 = R.TH1F('h_gen_dxy_pt_30',   'Denominator Gen Dxy for pt Thresh 30',	dxy_bins[0], dxy_bins[1], dxy_bins[2])


    h_gen_pt   = R.TH1F('h_gen_pt',   'Denominator Gen Pt',	pt_bins[0], pt_bins[1], pt_bins[2])
    h_gen_pt_pt10 = R.TH1F('h_gen_pt_pt_10',   'Denominator Gen Pt for pt Thresh 10',	pt_bins[0], pt_bins[1], pt_bins[2])
    h_gen_pt_pt20 = R.TH1F('h_gen_pt_pt_20',   'Denominator Gen Pt for pt Thresh 20',	pt_bins[0], pt_bins[1], pt_bins[2])

    h_eff_gen_pt   = R.TH1F('h_eff_pt_dxy',   'Eff numerator for Gen Pt',	pt_bins[0], pt_bins[1], pt_bins[2])
    h_eff_gen_pt_pt10   = R.TH1F('h_eff_gen_pt_pt10',   'Eff numerator for Gen Pt, pt > 10 GeV',	pt_bins[0], pt_bins[1], pt_bins[2])
    h_eff_gen_pt_pt20   = R.TH1F('h_eff_gen_pt_pt20',   'Eff numerator for Gen Pt, pt > 20 GeV',	pt_bins[0], pt_bins[1], pt_bins[2])
    
    #2D histos 

    h_emupt_gendxy      = R.TH2F('h_emupt_gendxy', 'p_{T} vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])

    h_emu_dxy_gen_dxy   = R.TH2F('h_emu_dxy_gen_dxy', 'Emu dxy vs Gen dxy', 4, 0., 4., dxy_bins[0], dxy_bins[1], dxy_bins[2])

#    h_uGTMu_pt_pt	= R.TH2F('h_uGTMu_pt_pt', 'uGT mu pt vs pt', pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])
#    h_uGTMu_ptvtx_ptvtx	= R.TH2F('h_uGTMu_ptvtx_ptvtx', 'uGT mu ptvtx vs ptvtx', pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])
#    h_uGTMu_eta_eta	= R.TH2F('h_uGTMu_eta_eta', 'uGT mu eta vs eta', eta_bins[0], eta_bins[1], eta_bins[2], eta_bins[0], eta_bins[1], eta_bins[2])
#    h_uGTMu_phi_phi	= R.TH2F('h_uGTMu_phi_phi', 'uGT mu phi vs phi', phi_bins[0], phi_bins[1], phi_bins[2], phi_bins[0], phi_bins[1], phi_bins[2])

#    h_uGTMu_pt_diff	= R.TH1F('h_uGTMu_pt_diff', 'uGT Mu relative pt difference', pt_bins[0], pt_bins[1], pt_bins[2])
#    h_uGTMu_ptvtx_diff	= R.TH1F('h_uGTMu_ptvtx_diff', 'uGT Mu relative pt difference', pt_bins[0], pt_bins[1], pt_bins[2])
#    h_uGTMu_phi_diff	= R.TH1F('h_uGTMu_phi_diff', 'uGT Mu relative phi difference', phi_bins[0], phi_bins[1], phi_bins[2])
#    h_uGTMu_eta_diff	= R.TH1F('h_uGTMu_eta_diff', 'uGT Mu relative eta difference', eta_bins[0], eta_bins[1], eta_bins[2])

###############

#For Pt debug

#    h_uGTMu_ptVtx_genpt       = R.TH2F('h_uGTMu_ptVtx_genpt', 'uGT ptVtx vs Gen pt', pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])
#    h_uGTMu_ptIVtx_genpt       = R.TH2F('h_uGTMu_ptIVtx_genpt', 'uGT ptIVtx vs Gen pt', pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])
#    h_uGTMu_ptDisp_genpt       = R.TH2F('h_uGTMu_ptDisp_genpt', 'uGT ptDisp vs Gen pt', pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])
#    h_uGTMu_ptIDisp_genpt       = R.TH2F('h_uGTMu_ptIDisp_genpt', 'uGT ptIDisp vs Gen pt', pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])


#Single mu histograms, comparing constrained vs unconstrained pt, vs dxy

    h_gen_dxy_den	= R.TH1F('h_gen_dxy_den', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])

    h_ptVtx_dxy_0 	= R.TH1F('h_ptVtx_dxy_0', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptDisp_dxy_0 	= R.TH1F('h_ptDisp_dxy_0', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptOr_dxy_0 	= R.TH1F('h_ptOr_dxy_0', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])

    h_ptVtx_dxy_4 	= R.TH1F('h_ptVtx_dxy_4', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptDisp_dxy_4 	= R.TH1F('h_ptDisp_dxy_4', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptOr_dxy_4 	= R.TH1F('h_ptOr_dxy_4', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])

    h_ptVtx_dxy_7 	= R.TH1F('h_ptVtx_dxy_7', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptDisp_dxy_7 	= R.TH1F('h_ptDisp_dxy_7', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptOr_dxy_7 	= R.TH1F('h_ptOr_dxy_7', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])

    h_ptVtx_dxy_11 	= R.TH1F('h_ptVtx_dxy_11', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptDisp_dxy_11 	= R.TH1F('h_ptDisp_dxy_11', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptOr_dxy_11 	= R.TH1F('h_ptOr_dxy_11', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])

    h_ptVtx_dxy_15 	= R.TH1F('h_ptVtx_dxy_15', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptDisp_dxy_15 	= R.TH1F('h_ptDisp_dxy_15', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_ptOr_dxy_15 	= R.TH1F('h_ptOr_dxy_15', 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])

    h_ptDisp_dxy_emu_pt_dxy_thresh = []

    for i in [0,4,7,11,15]:
	h_temp_list = []
	for j in range(4):
		h_temp = R.TH1F('h_ptDisp_dxy_pt_'+str(i)+'_emudxy_'+str(j), 'Efficiency vs Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
		h_temp_list.append(h_temp)
	h_ptDisp_dxy_emu_pt_dxy_thresh.append(h_temp_list)

    h_gen_Vx		= R.TH1F('h_gen_Vx', 	'Gen Muon Vx', 	vxyz_bins[0], vxyz_bins[1], vxyz_bins[2])
    h_gen_Vy		= R.TH1F('h_gen_Vy', 	'Gen Muon Vy', 	vxyz_bins[0], vxyz_bins[1], vxyz_bins[2])
    h_gen_Vz		= R.TH1F('h_gen_Vz', 	'Gen Muon Vz', 	vxyz_bins[0], vxyz_bins[1], vxyz_bins[2])
    h_gen_eta		= R.TH1F('h_gen_eta', 	'Gen Muon Eta', eta_bins[0], eta_bins[1], eta_bins[2])
    h_gen_phi		= R.TH1F('h_gen_phi', 	'Gen Muon Phi', phi_bins[0], phi_bins[1], phi_bins[2])
#    h_gen_ptinv   	= R.TH1F('h_gen_ptinv', 'Gen 1/Pt',	120, -0.6, 0.6)

    ### Direct TEfficiency plots

    h_Teff_ptDisp_eta_10 = R.TEfficiency('h_Teff_ptDisp_eta_10', 'Efficiency vs Eta', len(EC_eta_bins)-1, array('d',EC_eta_bins))
    h_Teff_ptVtx_eta_10 = R.TEfficiency('h_Teff_ptVtx_eta_10', 'Efficiency vs Eta', len(EC_eta_bins)-1, array('d',EC_eta_bins))
    h_Teff_ptOr_eta_10 = R.TEfficiency('h_Teff_ptOr_eta_10', 'Efficiency vs Eta', len(EC_eta_bins)-1, array('d',EC_eta_bins))

    h_Teff_ptDisp_pt_10 = R.TEfficiency('h_Teff_ptDisp_pt_10', 'Efficiency vs Pt', pt_bins[0], pt_bins[1], pt_bins[2])
    h_Teff_ptVtx_pt_10 = R.TEfficiency('h_Teff_ptVtx_pt_10', 'Efficiency vs Pt', pt_bins[0], pt_bins[1], pt_bins[2])
    h_Teff_ptOr_pt_10 = R.TEfficiency('h_Teff_ptOr_pt_10', 'Efficiency vs Pt', pt_bins[0], pt_bins[1], pt_bins[2])

    h_Teff_ptDisp_dxy_10 = R.TEfficiency('h_Teff_ptDisp_dxy_10', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_Teff_ptVtx_dxy_10 = R.TEfficiency('h_Teff_ptVtx_dxy_10', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_Teff_ptOr_dxy_10 = R.TEfficiency('h_Teff_ptOr_dxy_10', 'Efficiency vs Gen Dxy', dxy_bins[0], dxy_bins[1], dxy_bins[2])

#Dimuon histograms

    ##Underlying distributions and correlations in signal MC at gen

    h_dimu_dxy1_dxy2	= R.TH2F('h_dimu_dxy1_dxy2', 	'Lead dxy vs sublead dxy; Lead Dxy [cm]; Sublead Dxy [cm]', 	dxy_bins[0], dxy_bins[1], dxy_bins[2], dxy_bins[0], dxy_bins[1], dxy_bins[2])
    h_dimu_pt1_pt2	= R.TH2F('h_dimu_pt1_pt2', 	'Lead Pt vs sublead Pt; Lead Pt [GeV]; Sublead Pt [GeV]', 	pt_bins[0], pt_bins[1], pt_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])
    h_dimu_dxy1_pt2	= R.TH2F('h_dimu_dxy1_pt2', 	'Lead dxy vs sublead Pt; Lead Dxy [cm]; Sublead Pt [GeV]', 	dxy_bins[0], dxy_bins[1], dxy_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])
    h_dimu_dxy2_pt1	= R.TH2F('h_dimu_dxy2_pt1', 	'Sublead dxy vs lead pt; Sublead Dxy [cm]; Lead Pt [GeV]', 	dxy_bins[0], dxy_bins[1], dxy_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])
    h_dimu_dxy1_pt1	= R.TH2F('h_dimu_dxy1_pt1', 	'Lead dxy vs Lead Pt; Lead Dxy [cm]; Lead Pt [GeV]', 		dxy_bins[0], dxy_bins[1], dxy_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])
    h_dimu_dxy2_pt2	= R.TH2F('h_dimu_dxy2_pt2', 	'Sublead dxy vs sublead Pt; Sublead Dxy [cm]; Sublead Pt [GeV]',dxy_bins[0], dxy_bins[1], dxy_bins[2], pt_bins[0], pt_bins[1], pt_bins[2])

    #For efficiencies

    h_den_dxy1_dxy2 = []
    h_num_dxy1_dxy2 = []

    L1ptthresh = [0, 4, 7, 11, 15]

    h_blank_dxy1_dxy2 		= R.TH2F('h_blank_dxy1_dxy2_L1pt', 'Blank Lead Dxy vs Sublead Dxy for L1pt; Lead Gen Dxy [cm]; Sublead Gen Dxy [cm]', dxy_bins[0], dxy_bins[1], dxy_bins[2], dxy_bins[0], dxy_bins[1], dxy_bins[2])

    for i in range(len(L1ptthresh)):
	h_den_templist = []
	h_num_templist = []
	for j in range(len(L1ptthresh)):
		h_den_temp 	= R.TH2F('h_den_dxy1_dxy2_L1pt_'+str(L1ptthresh[i])+'_'+str(L1ptthresh[j]), 'Lead Dxy vs Sublead Dxy for L1pt = ('+str(L1ptthresh[i])+', '+str(L1ptthresh[j])+'); Lead Gen Dxy [cm]; Sublead Gen Dxy [cm]', dxy_bins[0], dxy_bins[1], dxy_bins[2], dxy_bins[0], dxy_bins[1], dxy_bins[2])
		h_den_templist.append(h_den_temp)
		h_num_temp 	= R.TH2F('h_num_dxy1_dxy2_L1pt_'+str(L1ptthresh[i])+'_'+str(L1ptthresh[j]), 'Lead Dxy vs Sublead Dxy for L1pt = ('+str(L1ptthresh[i])+', '+str(L1ptthresh[j])+'); Lead Gen Dxy [cm]; Sublead Gen Dxy [cm]', dxy_bins[0], dxy_bins[1], dxy_bins[2], dxy_bins[0], dxy_bins[1], dxy_bins[2])
		h_num_templist.append(h_num_temp)
	h_den_dxy1_dxy2.append(h_den_templist)
	h_num_dxy1_dxy2.append(h_num_templist)

    ####debug 

#Event counters

    L1_SingleMu20_BMTF 		= 0
    mucopy_count       		= 0
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
        Emu_br = R.L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
        uGT_br = R.L1Analysis.L1AnalysisL1UpgradeDataFormat()
        Gen_br = R.L1Analysis.L1AnalysisGeneratorDataFormat()
#        Kmt_br = R.L1Analysis.L1AnalysisBMTFOutputDataFormat()

        chains['Evt'][iCh].SetBranchAddress('Event',               R.AddressOf(Evt_br))
        chains['Unp'][iCh].SetBranchAddress('L1UpgradeBmtfMuon',   R.AddressOf(Unp_br))
        chains['Emu'][iCh].SetBranchAddress('L1UpgradeEmtfMuon',  R.AddressOf(Emu_br))
        chains['uGT'][iCh].SetBranchAddress('L1Upgrade',  	   R.AddressOf(uGT_br))
        chains['Gen'][iCh].SetBranchAddress('Generator',           R.AddressOf(Gen_br))
#        chains['Emu'][iCh].SetBranchAddress('L1UpgradeBmtfOutput', R.AddressOf(Kmt_br))

	

        print '\nEntering loop over events for chain %d' % iCh
        for jEvt in range(chains['Emu'][iCh].GetEntries()):

	    L1_SingleMu20_BMTF_flag = False

            if iEvt >= MAX_EVT: break

	    iEvt +=1

            if iEvt % PRT_EVT is 0: print '\nEvent # %d (%dth in chain)' % (iEvt, jEvt+1)

            chains['Evt'][iCh].GetEntry(jEvt)
            chains['Unp'][iCh].GetEntry(jEvt)
            chains['Emu'][iCh].GetEntry(jEvt)
            chains['uGT'][iCh].GetEntry(jEvt)
            chains['Gen'][iCh].GetEntry(jEvt)

            # ## Use these lines if you don't explicitly define the DataFormat and then do SetBranchAddress above
            # Evt_br = chains['Evt'][iCh].Event
            # Unp_br = chains['Unp'][iCh].L1UpgradeBmtfMuon
            # Emu_br = chains['Emu'][iCh].L1UpgradeBmtfMuon

            if iEvt % PRT_EVT is 0: print '  * Run %d, LS %d, event %d' % (int(Evt_br.run), int(Evt_br.lumi), int(Evt_br.event))

            nUnpMu = int(Unp_br.nTfMuons)
            nEmuMu = int(Emu_br.nTfMuons)
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
	    EMTF_Mus			=[]
	    EMTF_Mus_unique		=[]
	
	    #########################################
            ###  Generator information for muons  ###
            #########################################

	    for i in range(nGenPart):
	
		if (abs(Gen_br.partId[i])!=13): continue
		if (Gen_br.partStat[i]!=1): continue

		### weighting MuGun events ###
		
		weight = 1.
	
		pt = float(Gen_br.partPt[i])

		##############################

		if abs(float(Gen_br.partVx[i])) > 120. : continue
		if abs(float(Gen_br.partVy[i])) > 120. : continue
		if abs(float(Gen_br.partVz[i])) > 100. : continue
		if abs(float(Gen_br.partEta[i])) > 3. or abs(float(Gen_br.partEta[i])) < 0.6 : continue
		if pt < 2 : continue

		h_gen_Vx.Fill(float(Gen_br.partVx[i]), weight)
		h_gen_Vy.Fill(float(Gen_br.partVy[i]), weight)
		h_gen_Vz.Fill(float(Gen_br.partVz[i]), weight)
		h_gen_eta.Fill(float(Gen_br.partEta[i]), weight)
		h_gen_phi.Fill(float(Gen_br.partPhi[i]), weight)

#		if (Gen_br.partParent[i]!=6000113): continue

		GenMus.append(i)

#		if IsinacceptBarrel(Gen_br, i) : MuType = 1 
		if IsinacceptEndcap(Gen_br, i) : MuType = 3
		else : 	continue

		GenMus_acc.append([i, MuType])
		
		GenMuPt_acc.append([Gen_br.partPt[i], i, MuType])
		
		dxy = getDxy(float(Gen_br.partVx[i]), float(Gen_br.partVy[i]), float(Gen_br.partPhi[i]))
		pt  = float(Gen_br.partPt[i])
		

		eta, phi = getGenEtaPhi(Gen_br, i, MuType)

		h_gen_pt.Fill(pt)
		h_gen_dxy_den.Fill(dxy, weight)

		if Gen_br.partPt[i] >= 5. :
			h_gen_dxy.Fill(dxy, weight)
		if Gen_br.partPt[i] >= 15. :
			h_gen_dxy_pt10.Fill(dxy, weight)
			h_gen_pt_pt10.Fill(pt, weight)
		if Gen_br.partPt[i] >= 25. :
			h_gen_dxy_pt20.Fill(dxy, weight)
			h_gen_pt_pt20.Fill(pt, weight)
		if Gen_br.partPt[i] >= 35. :
			h_gen_dxy_pt30.Fill(dxy, weight)

	    GenMuPt_acc.sort(reverse=True)

            #################################
            ###  Unpacked (legacy) muons  ###
            #################################
            for i in range(nUnpMu):
                BX      = int(Unp_br.tfMuonBx[i])
                qual    = int(Unp_br.tfMuonHwQual[i])
                ptVtx   = float(Unp_br.tfMuonHwPt[i] - 1)*0.5   ## Vertex-constrained (standard) pT is stored in 0.5 GeV steps
#                ptDispl = float(Unp_br.tfMuonHwPtUnconstrained[i] - 1)  ## Is there an offset by 1 for displaced muons? - AWB 2019.05.29, ALB: it is empty.
                eta     = float(Unp_br.tfMuonHwEta[i])*0.010875
                
                if (BX  !=  0): continue
                if (qual < 12): continue
		if (ptVtx>=20. and abs(eta)<0.83): 
			L1_SingleMu20_BMTF_flag = True


            #################################
            ###  Emulated uGT muons  ###
            #################################

	    mucopy_flag = False		
		
            for i in range(nuGTMu):
                BX      = int(uGT_br.muonBx[i])
                qual    = int(uGT_br.muonQual[i])
                ptVtx   = float(uGT_br.muonEt[i])  ## Vertex-constrained (standard) pT is stored in 0.5 GeV steps
                ptDisp  = float(uGT_br.muonEtUnconstrained[i])  ## Vertex-unconstrained pT

#		if ptDisp == 0 : continue
		
		mu1 = R.TLorentzVector()
		mu1.SetPtEtaPhiM(ptDisp, float(uGT_br.muonEtaAtVtx[i]), float(uGT_br.muonPhiAtVtx[i]), 105.7e-3)
		eta = mu1.Eta()
		phi = mu1.Phi()
                                
                if (BX  !=  0): continue
#                if (qual < 11): continue

		EmuMus.append(i)

		savemuindex = i
		ptmax = -1.
		
		for j in range(i+1, nuGTMu):
			mu2 = R.TLorentzVector()
			mu2.SetPtEtaPhiM(float(uGT_br.muonEtUnconstrained[j]), float(uGT_br.muonEtaAtVtx[j]), float(uGT_br.muonPhiAtVtx[j]), 105.7e-3)
			pt_j   = max(float(uGT_br.muonEtUnconstrained[j]), float(uGT_br.muonEt[j]))
			eta_j  = mu2.Eta()
			phi_j  = mu2.Phi()

			if (eta == eta_j and phi == phi_j):
				mucopy_flag = True
				
				if (ptmax < pt_j): 
					ptmax = pt_j
					savemuindex = j
				else :	savemuindex = i

		if i != savemuindex: continue

		## Only for EMTF

		if abs(eta) <= 1.245 or abs(eta) >= 2.45: continue
	
		EmuMus_unique.append(i)

	    #Flagging events with nearly identical Emulated muons

	    if mucopy_flag:
		mucopy_count +=1

            #################################
            ###  Emulated EMTF muons  ###
            #################################

            for i in range(nEmuMu):
                BX      = int(Emu_br.tfMuonBx[i])
                qual    = int(Emu_br.tfMuonHwQual[i])
                ptVtx   = float(Emu_br.tfMuonHwPt[i]-1.)*0.5  ## Vertex-constrained (standard) pT is stored in 0.5 GeV steps
                ptDisp  = float(Emu_br.tfMuonHwPtUnconstrained[i]-1.)  ## Vertex-unconstrained pT

#		if ptDisp == 0 : continue
		
		mu1 = R.TLorentzVector()
		mu1.SetPtEtaPhiM(ptDisp, float(Emu_br.tfMuonHwEta[i])*0.010875, getPhi(float(Emu_br.tfMuonGlobalPhi[i])), 105.7e-3)
		eta = mu1.Eta()
		phi = mu1.Phi()
                                
                if (BX  !=  0): continue
		if float(Emu_br.tfMuonHwPtUnconstrained[i] == 0.) : continue
#                if (qual < 11): continue

		EMTF_Mus.append(i)

		savemuindex = i
		ptmax = -1.
		
		for j in range(i+1, nEmuMu):
			mu2 = R.TLorentzVector()
			mu2.SetPtEtaPhiM(ptDisp, float(Emu_br.tfMuonHwEta[j])*0.010875, getPhi(float(Emu_br.tfMuonGlobalPhi[j])), 105.7e-3)
			pt_j   = max(float(Emu_br.tfMuonHwPtUnconstrained[j]-1.), float(Emu_br.tfMuonHwPt[j]-1.)*0.5)
			eta_j  = mu2.Eta()
			phi_j  = mu2.Phi()

			if (eta == eta_j and phi == phi_j):
				mucopy_flag = True
				
				if (ptmax < pt_j): 
					ptmax = pt_j
					savemuindex = j
				else :	savemuindex = i

		if i != savemuindex: continue

		## Only for EMTF

                if abs(eta) < 1.2 or abs(eta) > 2.5: continue

                EMTF_Mus_unique.append(i)

	    # Better matching subroutine

	    dR_list=[]
	    dR_list_corr=[]
	    EMTF_dR_list=[]

	    for el in GenMus_acc:
		
		i = el[0]
		mutype = el[1]
		
		GenMuCorr = R.TLorentzVector()
		temp_eta, temp_phi = getGenEtaPhi(Gen_br, i, mutype)
		GenMuCorr.SetPtEtaPhiM(float(Gen_br.partPt[i]), temp_eta, temp_phi, 105.7e-3)
		eta1corr = GenMuCorr.Eta()
		phi1corr = GenMuCorr.Phi()

		for j in EmuMus_unique:
			EmuMuvec = R.TLorentzVector()
			EmuMuvec.SetPtEtaPhiM(float(uGT_br.muonEtUnconstrained[j]), float(uGT_br.muonEtaAtVtx[j]), float(uGT_br.muonPhiAtVtx[j]), 105.7e-3)
	
			dR_corr = EmuMuvec.DeltaR(GenMuCorr)
			if dR_corr < 0.6: dR_list_corr.append([dR_corr, i, j, mutype])


	    for el in EMTF_Mus_unique:
		EmuMuvec = R.TLorentzVector()
		EmuMuvec.SetPtEtaPhiM(float(Emu_br.tfMuonHwPtUnconstrained[el]-1.), float(Emu_br.tfMuonHwEta[el])*0.010875, getPhi(float(Emu_br.tfMuonGlobalPhi[el])) , 105.7e-3)
		for j in EmuMus_unique:
                        EmuMuvec2 = R.TLorentzVector()
                        EmuMuvec2.SetPtEtaPhiM(float(uGT_br.muonEtUnconstrained[j]), float(uGT_br.muonEtaAtVtx[j]), float(uGT_br.muonPhiAtVtx[j]), 105.7e-3)

                        dR = EmuMuvec2.DeltaR(EmuMuvec)
                        if dR < 0.5: EMTF_dR_list.append([dR, el, j])

#	    mismatch_flag = False

	    if len(dR_list_corr):
#		print "Presorted dR_list: ", dR_list_corr
		dR_list_corr.sort()

	    dR_list = dR_list_corr

	    k = False
	    if len(dR_list): k = True
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

	    if len(EMTF_dR_list):
		
#		print "Presorted EMTF_dR_list: ", EMTF_dR_list

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

#	    if len(dR_list) :	    print "Cleaned sorted dR_list: ", dR_list
#	    if len(EMTF_dR_list) :	    print "Cleaned sorted EMTF_dR_list: ", EMTF_dR_list

	    #Dxy, pt for matched gen muon
	    for i in range(len(dR_list)):
		Emu_idx = dR_list[i][2]
		ptVtx 	= float(uGT_br.muonEt[Emu_idx])
		ptIVtx 	= float(uGT_br.muonIEt[Emu_idx])
		ptDisp	= float(uGT_br.muonEtUnconstrained[Emu_idx])
		ptIDisp	= float(uGT_br.muonIEtUnconstrained[Emu_idx])
		etaEmu 	= float(uGT_br.muonEtaAtVtx[Emu_idx])
#		if float(uGT_br.muonIEtUnconstrained[Emu_idx]) == 0 : ptDisp = 0.
#		else :	ptDisp	= float(uGT_br.muonIEtUnconstrained[Emu_idx]-1.)
		Emu_dxy = float(uGT_br.muonDxy[Emu_idx])

		pt =  max(ptDisp, ptVtx)
		idx	= dR_list[i][1]
		dR  	= dR_list[i][0]
		mutype 	= dR_list[i][3]
		dxy 	= getDxy(float(Gen_br.partVx[idx]), float(Gen_br.partVy[idx]), float(Gen_br.partPhi[idx]))

		matched_gen_mu_count+=1
		pt_gen = float(Gen_br.partPt[idx])		
		q_gen	= float(Gen_br.partCh[idx])
		q_emu	= float(uGT_br.muonChg[Emu_idx])

		### weighting MuGun events ###
		
		weight = 1.

		##############################

		if pt :
			if pt_gen >= 5. : 
				h_eff_gen_dxy.Fill(dxy, weight)
				h_emupt_gendxy.Fill(dxy, pt, weight)
			h_eff_gen_pt.Fill(pt_gen, weight)
		if pt >= 10. :
			if pt_gen >= 15. : 
				h_eff_gen_dxy_pt10.Fill(dxy, weight)
			h_eff_gen_pt_pt10.Fill(pt_gen, weight)
		if pt >= 20. :
			if pt_gen >= 25. : 
				h_eff_gen_dxy_pt20.Fill(dxy, weight)
			h_eff_gen_pt_pt20.Fill(pt_gen, weight)

		if pt >= 30. :
			if pt_gen >= 35. : 
				h_eff_gen_dxy_pt30.Fill(dxy, weight)

	    	emu_dxy = float(uGT_br.muonDxy[dR_list[i][2]])
		h_emu_dxy_gen_dxy.Fill(emu_dxy, dxy, weight)

		if ptVtx >= 0. : h_ptVtx_dxy_0.Fill(dxy, weight)
		if ptVtx >= 4. : h_ptVtx_dxy_4.Fill(dxy, weight)
		if ptVtx >= 7. : h_ptVtx_dxy_7.Fill(dxy, weight)
		if ptVtx >= 11. : h_ptVtx_dxy_11.Fill(dxy, weight)
		if ptVtx >= 15. : h_ptVtx_dxy_15.Fill(dxy, weight)

		if ptDisp >= 0. : h_ptDisp_dxy_0.Fill(dxy, weight)
		if ptDisp >= 4. : h_ptDisp_dxy_4.Fill(dxy, weight)
		if ptDisp >= 7. : h_ptDisp_dxy_7.Fill(dxy, weight)
		if ptDisp >= 11. : h_ptDisp_dxy_11.Fill(dxy, weight)
		if ptDisp >= 15. : h_ptDisp_dxy_15.Fill(dxy, weight)

		if ptVtx >= 0 or ptDisp >= 0 : h_ptOr_dxy_0.Fill(dxy, weight)
		if ptVtx >= 4 or ptDisp >= 4 : h_ptOr_dxy_4.Fill(dxy, weight)
		if ptVtx >= 7 or ptDisp >= 7 : h_ptOr_dxy_7.Fill(dxy, weight)
		if ptVtx >= 11 or ptDisp >= 11 : h_ptOr_dxy_11.Fill(dxy, weight)
		if ptVtx >= 15 or ptDisp >= 15 : h_ptOr_dxy_15.Fill(dxy, weight)
	
		ptthreshlist = [0,4,7,11,15]
		

		if dxy < 100. and abs(float(Gen_br.partVz[idx])) < 100. and pt_gen >= 2.:


			eta, phi = getGenEtaPhi(Gen_br, idx, mutype)
			for k in range(5):
				for j in range(4) :
					if ptDisp >= ptthreshlist[k] and Emu_dxy >= j : h_ptDisp_dxy_emu_pt_dxy_thresh[k][j].Fill(dxy, weight)
			
			#Relative L1 pt variables to gen pt
#			h_ptVtx_genpt_rel.Fill(pt_gen, (ptVtx-pt_gen)/pt_gen)
#			h_ptIVtx_genpt_rel.Fill(pt_gen, (ptIVtx-pt_gen)/pt_gen)
#			h_ptDisp_genpt_rel.Fill(pt_gen, (ptDisp-pt_gen)/pt_gen)
#			h_ptIDisp_genpt_rel.Fill(pt_gen, (ptIDisp-pt_gen)/pt_gen)

#			curvVtx 	= (q_emu/ptVtx - q_gen/pt_gen)/(q_gen/pt_gen)
#			curvIVtx 	= (q_emu/ptIVtx - q_gen/pt_gen)/(q_gen/pt_gen)
#			curvDisp 	= (q_emu/ptDisp - q_gen/pt_gen)/(q_gen/pt_gen)
#			curvIDisp 	= (q_emu/ptIDisp - q_gen/pt_gen)/(q_gen/pt_gen)
# 			curvGen		= q_gen/pt_gen
		
#			h_ptVtx_genpt_dxy.Fill(dxy, curvVtx)	
#			h_ptIVtx_genpt_dxy.Fill(dxy, curvIVtx)	
#			h_ptDisp_genpt_dxy.Fill(dxy, curvDisp)	
#			h_ptIDisp_genpt_dxy.Fill(dxy, curvIDisp)	

#			h_ptVtx_genpt_qrelpt.Fill(pt_gen, curvVtx)	
#			h_ptIVtx_genpt_qrelpt.Fill(pt_gen, curvIVtx)	
#			h_ptDisp_genpt_qrelpt.Fill(pt_gen, curvDisp)	
#			h_ptIDisp_genpt_qrelpt.Fill(pt_gen, curvIDisp)	
	
			#L1 Pt variables for debug
#			h_ptVtx_debug_denom.Fill(ptVtx)
#			h_ptIVtx_debug_denom.Fill(ptIVtx)
#			h_ptDisp_debug_denom.Fill(ptDisp)
#			h_ptIDisp_debug_denom.Fill(ptIDisp)
#			h_genpt_debug_denom.Fill(pt_gen)

#			h_uGTMu_ptVtx_genpt.Fill(pt_gen, ptVtx)
#			h_uGTMu_ptIVtx_genpt.Fill(pt_gen, ptIVtx)
#			h_uGTMu_ptDisp_genpt.Fill(pt_gen, ptDisp)
#			h_uGTMu_ptIDisp_genpt.Fill(pt_gen, ptIDisp)

			if abs(eta) >= 1.2 and abs(eta) <= 1.6 and  ptDisp == 0 	:
				for el in EMTF_dR_list:
					if Emu_idx == el[2] : ptDisp = float(Emu_br.tfMuonHwPtUnconstrained[el[1]]-1.)
					

			h_Teff_ptVtx_pt_10.FillWeighted(0, weight, pt_gen)
			h_Teff_ptDisp_pt_10.FillWeighted(0, weight, pt_gen)
			h_Teff_ptOr_pt_10.FillWeighted(0, weight, pt_gen)

			h_Teff_ptVtx_eta_10.Fill(0, eta)
			h_Teff_ptDisp_eta_10.Fill(0, eta)
			h_Teff_ptOr_eta_10.Fill(0, eta)

			h_Teff_ptVtx_dxy_10.Fill(0, dxy)
			h_Teff_ptDisp_dxy_10.Fill(0, dxy)
			h_Teff_ptOr_dxy_10.Fill(0, dxy)

			if ptVtx >= 10. : 

				h_Teff_ptVtx_pt_10.FillWeighted(1, weight, pt_gen)
				h_Teff_ptVtx_eta_10.Fill(1, eta)
				h_Teff_ptVtx_dxy_10.Fill(1, dxy)

			if ptDisp >= 10. : 

				h_Teff_ptDisp_pt_10.FillWeighted(1, weight, pt_gen)
				h_Teff_ptDisp_eta_10.Fill(1, eta)
				h_Teff_ptDisp_dxy_10.Fill(1, dxy)

			if ptDisp >= 10. or ptVtx >= 10.: 

				h_Teff_ptOr_pt_10.FillWeighted(1, weight, pt_gen)
				h_Teff_ptOr_eta_10.Fill(1, eta)
				h_Teff_ptOr_dxy_10.Fill(1, dxy)

	    #end of loop over matched muons

	    #Efficiency for Dimuons

	    Dimus_acc = []
	    
	    if len(GenMus_acc) >= 2 :
		for i in range(len(GenMus_acc)):
			idx1 = GenMus_acc[i][0]
			vx1 = float(Gen_br.partVx[idx1])
			vy1 = float(Gen_br.partVy[idx1])
			vz1 = float(Gen_br.partVz[idx1])

			pt1 = float(Gen_br.partPt[idx1])

			for j in range(i+1, len(GenMus_acc)):
				idx2 = GenMus_acc[j][0]
				vx2 = float(Gen_br.partVx[idx2])
				vy2 = float(Gen_br.partVy[idx2])
				vz2 = float(Gen_br.partVz[idx2])
				pt2 = float(Gen_br.partPt[idx2])

				if vx1==vx2 and vy1==vy2 and vz1==vz2:
					if pt1 >= pt2 	:	el = [idx1, idx2]
					else 		:	el = [idx2, idx1]
					Dimus_acc.append(el)

	    if len(Dimus_acc) :

#	Debug block
#
#		print 'Dimus in acceptance: ', Dimus_acc
#		if len(Dimus_acc) > 2 : 
#			for k in range(len(Dimus_acc)):
#				idx1 = Dimus_acc[k][0]
#				idx2 = Dimus_acc[k][1]
#				v1 = [float(Gen_br.partVx[idx1]),float(Gen_br.partVy[idx1]),float(Gen_br.partVz[idx1])]
#				v2 = [float(Gen_br.partVx[idx2]),float(Gen_br.partVy[idx2]),float(Gen_br.partVz[idx2])]
#				print 'Muon Coordinates:'
#				print 'Idx1 : ', idx1, '\t', v1
#				print 'Idx2 : ', idx2, '\t', v2


	    	matched_dimus=[]

		GenDimu_Evt +=1 
		for j in range(len(Dimus_acc)):
			temp = []
			for i in range(len(dR_list)):
				Gen_idx = dR_list[i][1]
				if Gen_idx in Dimus_acc[j] :
					el = [Gen_idx, dR_list[i][2]]
					if el not in temp:
						temp.append(el)
			if len(temp) == 2 : matched_dimus.append(temp)
			
			idx1 = Dimus_acc[j][0]
			idx2 = Dimus_acc[j][1]
			print 'Dimu pair : ', Dimus_acc[j]
			print 'Muon Index:'
			print 'Idx1 : ', idx1
			print 'Idx2 : ', idx2
	
			dxy1 = getDxy(float(Gen_br.partVx[idx1]), float(Gen_br.partVy[idx1]), float(Gen_br.partPhi[idx1]))
			dxy2 = getDxy(float(Gen_br.partVx[idx2]), float(Gen_br.partVy[idx2]), float(Gen_br.partPhi[idx2]))

			pt1 = float(Gen_br.partPt[idx1])
			pt2 = float(Gen_br.partPt[idx2])

			h_dimu_dxy1_dxy2.Fill(dxy1, dxy2)
			h_dimu_pt1_pt2.Fill(pt1, pt2)
			h_dimu_dxy1_pt2.Fill(dxy1, pt2)
			h_dimu_dxy2_pt1.Fill(dxy2, pt1)
			h_dimu_dxy1_pt1.Fill(dxy1, pt1)
			h_dimu_dxy2_pt2.Fill(dxy2, pt2)

#	    	if len(matched_dimus):
#			print 'Matched dimuon indices: ', len(matched_dimus), "\t", matched_dimus

		for dim_pair in matched_dimus:
			genidx1 = dim_pair[0][0]
			genidx2 = dim_pair[1][0]
			emuidx1 = dim_pair[0][1]
			emuidx2 = dim_pair[1][1]

			print 'Matched Muon Index:'
			print 'Idx1 : ', genidx1
			print 'Idx2 : ', genidx2

			dxy1 = getDxy(float(Gen_br.partVx[genidx1]),float(Gen_br.partVy[genidx1]),float(Gen_br.partPhi[genidx1]))
			dxy2 = getDxy(float(Gen_br.partVx[genidx2]),float(Gen_br.partVy[genidx2]),float(Gen_br.partPhi[genidx2]))
	
			genpt1 = float(Gen_br.partPt[genidx1])
			genpt2 = float(Gen_br.partPt[genidx2])

			emupt1 = float(uGT_br.muonEtUnconstrained[emuidx1])
			emupt2 = float(uGT_br.muonEtUnconstrained[emuidx2])

			ptlead = max(emupt1, emupt2)
			ptsublead = min(emupt1, emupt2)

			for l in range(1, len(L1ptthresh)):
				L1pt1 = L1ptthresh[l]
				for m in range(l):
					L1pt2 = L1ptthresh[m]
					if genpt1 >= L1pt1 + 5. and genpt2 >= L1pt2 + 5.	:	
						h_den_dxy1_dxy2[l][m].Fill(dxy1, dxy2)
						if ptlead >= L1pt1 and ptsublead >= L1pt2	:
							h_num_dxy1_dxy2[l][m].Fill(dxy1, dxy2)


	    ##End of sub loop over Gen, Emu or Unp muons

#	    print "GenMus: ", len(GenMus)
#	    print "GenMusAcc: ", len(GenMus_acc)
#	    print "EmuMus: ", len(EmuMus)
#	    print "EmuMus unique: ", len(EmuMus_unique)
	
	    ###########################
            ###  Updating counters  ###
            ###########################
 	    
	    if (L1_SingleMu20_BMTF_flag): L1_SingleMu20_BMTF +=1

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

    h_gen_dxy.SetLineWidth(2)
    h_gen_dxy.Write()
    
    h_gen_pt.SetLineWidth(2)
    h_gen_pt.Write()

    h_emu_dxy_gen_dxy.SetLineWidth(2)
    h_emu_dxy_gen_dxy.Write()

    h_eff_gen_dxy.Sumw2()
    h_Teff_gen_dxy = R.TEfficiency(h_eff_gen_dxy, h_gen_dxy)
    h_eff_gen_dxy.Divide(h_eff_gen_dxy, h_gen_dxy, 1.,1., "B")
    h_Teff_gen_dxy.SetLineWidth(2)
    h_eff_gen_dxy.SetLineWidth(2)
    h_eff_gen_dxy.Write()

    h_eff_gen_dxy_pt10.Sumw2()
    h_Teff_gen_dxy_pt10 = R.TEfficiency(h_eff_gen_dxy_pt10, h_gen_dxy_pt10)
    h_eff_gen_dxy_pt10.Divide(h_eff_gen_dxy_pt10, h_gen_dxy_pt10, 1., 1., "B")
    h_Teff_gen_dxy_pt10.SetLineWidth(2)
    h_eff_gen_dxy_pt10.SetLineWidth(2)
    h_eff_gen_dxy_pt10.Write()

    h_eff_gen_dxy_pt20.Sumw2()
    h_Teff_gen_dxy_pt20 = R.TEfficiency(h_eff_gen_dxy_pt20, h_gen_dxy_pt20)
    h_eff_gen_dxy_pt20.Divide(h_eff_gen_dxy_pt20, h_gen_dxy_pt20, 1., 1., "B")
    h_Teff_gen_dxy_pt20.SetLineWidth(2)
    h_eff_gen_dxy_pt20.SetLineWidth(2)
    h_eff_gen_dxy_pt20.Write()
    
    h_eff_gen_dxy_pt30.Sumw2()
    h_Teff_gen_dxy_pt30 = R.TEfficiency(h_eff_gen_dxy_pt30, h_gen_dxy_pt30)
    h_eff_gen_dxy_pt30.Divide(h_eff_gen_dxy_pt30, h_gen_dxy_pt30, 1., 1., "B")
    h_Teff_gen_dxy_pt30.SetLineWidth(2)
    h_eff_gen_dxy_pt30.SetLineWidth(2)
    h_eff_gen_dxy_pt30.Write()

    h_eff_gen_pt.Sumw2()
    h_eff_gen_pt.Divide(h_eff_gen_pt, h_gen_pt, 1.,1., "B")
    h_eff_gen_pt.SetLineWidth(2)
    h_eff_gen_pt.Write()

    h_eff_gen_pt_pt10.Sumw2()
    h_eff_gen_pt_pt10.Divide(h_eff_gen_pt_pt10, h_gen_pt, 1.,1., "B")
    h_eff_gen_pt_pt10.SetLineWidth(2)
    h_eff_gen_pt_pt10.Write()

    h_eff_gen_pt_pt20.Sumw2()
    h_eff_gen_pt_pt20.Divide(h_gen_pt, h_eff_gen_pt_pt20, 1., 1., "B")
    h_eff_gen_pt_pt20.SetLineWidth(2)
    h_eff_gen_pt_pt20.Write()

    h_Teff_gen_dxy.Write()
    h_Teff_gen_dxy_pt10.Write()
    h_Teff_gen_dxy_pt20.Write()
    h_Teff_gen_dxy_pt30.Write()

####################################

    h_gen_dxy_den.GetXaxis().SetTitle("Gen Dxy [cm]")
    h_gen_dxy_den.GetYaxis().SetTitle("Efficiency")
    h_gen_dxy_den.GetYaxis().SetRangeUser(0, 1.2)

    h_Teff_ptVtx_dxy_0 = R.TEfficiency(h_ptVtx_dxy_0, h_gen_dxy_den)
    h_Teff_ptDisp_dxy_0 = R.TEfficiency(h_ptDisp_dxy_0, h_gen_dxy_den)
    h_Teff_ptOr_dxy_0 = R.TEfficiency(h_ptOr_dxy_0, h_gen_dxy_den)

    h_Teff_ptVtx_dxy_4 = R.TEfficiency(h_ptVtx_dxy_4, h_gen_dxy_den)
    h_Teff_ptDisp_dxy_4 = R.TEfficiency(h_ptDisp_dxy_4, h_gen_dxy_den)
    h_Teff_ptOr_dxy_4 = R.TEfficiency(h_ptOr_dxy_4, h_gen_dxy_den)

    h_Teff_ptVtx_dxy_7 = R.TEfficiency(h_ptVtx_dxy_7, h_gen_dxy_den)
    h_Teff_ptDisp_dxy_7 = R.TEfficiency(h_ptDisp_dxy_7, h_gen_dxy_den)
    h_Teff_ptOr_dxy_7 = R.TEfficiency(h_ptOr_dxy_7, h_gen_dxy_den)

    h_Teff_ptVtx_dxy_11 = R.TEfficiency(h_ptVtx_dxy_11, h_gen_dxy_den)
    h_Teff_ptDisp_dxy_11 = R.TEfficiency(h_ptDisp_dxy_11, h_gen_dxy_den)
    h_Teff_ptOr_dxy_11 = R.TEfficiency(h_ptOr_dxy_11, h_gen_dxy_den)

    h_Teff_ptVtx_dxy_15 = R.TEfficiency(h_ptVtx_dxy_15, h_gen_dxy_den)
    h_Teff_ptDisp_dxy_15 = R.TEfficiency(h_ptDisp_dxy_15, h_gen_dxy_den)
    h_Teff_ptOr_dxy_15 = R.TEfficiency(h_ptOr_dxy_15, h_gen_dxy_den)


    h_Teff_ptVtx_dxy_0.SetName("h_Eff_ptVtx_0_gen_dxy")
    h_Teff_ptVtx_dxy_0.SetLineWidth(2)
    h_Teff_ptVtx_dxy_0.SetLineColor(R.kRed)
    h_Teff_ptDisp_dxy_0.SetName("h_Eff_ptDisp_0_gen_dxy")
    h_Teff_ptDisp_dxy_0.SetLineWidth(2)
    h_Teff_ptDisp_dxy_0.SetLineColor(R.kBlue)
    h_Teff_ptOr_dxy_0.SetName("h_Eff_ptVtxOrDisp_0_gen_dxy")
    h_Teff_ptOr_dxy_0.SetLineWidth(2)
    h_Teff_ptOr_dxy_0.SetLineColor(R.kGreen)

    h_Teff_ptVtx_dxy_4.SetName("h_Eff_ptVtx_4_gen_dxy")
    h_Teff_ptVtx_dxy_4.SetLineWidth(2)
    h_Teff_ptVtx_dxy_4.SetLineColor(R.kRed)
    h_Teff_ptDisp_dxy_4.SetName("h_Eff_ptDisp_4_gen_dxy")
    h_Teff_ptDisp_dxy_4.SetLineWidth(2)
    h_Teff_ptDisp_dxy_4.SetLineColor(R.kBlue)
    h_Teff_ptOr_dxy_4.SetName("h_Eff_ptVtxOrDisp_4_gen_dxy")
    h_Teff_ptOr_dxy_4.SetLineWidth(2)
    h_Teff_ptOr_dxy_4.SetLineColor(R.kGreen)

    h_Teff_ptVtx_dxy_7.SetName("h_Eff_ptVtx_7_gen_dxy")
    h_Teff_ptVtx_dxy_7.SetLineWidth(2)
    h_Teff_ptVtx_dxy_7.SetLineColor(R.kRed)
    h_Teff_ptDisp_dxy_7.SetName("h_Eff_ptDisp_7_gen_dxy")
    h_Teff_ptDisp_dxy_7.SetLineWidth(2)
    h_Teff_ptDisp_dxy_7.SetLineColor(R.kBlue)
    h_Teff_ptOr_dxy_7.SetName("h_Eff_ptVtxOrDisp_7_gen_dxy")
    h_Teff_ptOr_dxy_7.SetLineWidth(2)
    h_Teff_ptOr_dxy_7.SetLineColor(R.kGreen)

    h_Teff_ptVtx_dxy_11.SetName("h_Eff_ptVtx_11_gen_dxy")
    h_Teff_ptVtx_dxy_11.SetLineWidth(2)
    h_Teff_ptVtx_dxy_11.SetLineColor(R.kRed)
    h_Teff_ptDisp_dxy_11.SetName("h_Eff_ptDisp_11_gen_dxy")
    h_Teff_ptDisp_dxy_11.SetLineWidth(2)
    h_Teff_ptDisp_dxy_11.SetLineColor(R.kBlue)
    h_Teff_ptOr_dxy_11.SetName("h_Eff_ptVtxOrDisp_11_gen_dxy")
    h_Teff_ptOr_dxy_11.SetLineWidth(2)
    h_Teff_ptOr_dxy_11.SetLineColor(R.kGreen)

    h_Teff_ptVtx_dxy_15.SetName("h_Eff_ptVtx_15_gen_dxy")
    h_Teff_ptVtx_dxy_15.SetLineWidth(2)
    h_Teff_ptVtx_dxy_15.SetLineColor(R.kRed)
    h_Teff_ptDisp_dxy_15.SetName("h_Eff_ptDisp_15_gen_dxy")
    h_Teff_ptDisp_dxy_15.SetLineWidth(2)
    h_Teff_ptDisp_dxy_15.SetLineColor(R.kBlue)
    h_Teff_ptOr_dxy_15.SetName("h_Eff_ptVtxOrDisp_15_gen_dxy")
    h_Teff_ptOr_dxy_15.SetLineWidth(2)
    h_Teff_ptOr_dxy_15.SetLineColor(R.kGreen)

    h_Teff_ptDisp_dxy_emu_pt_dxy_thresh = []

    ptthreshlist =[0,4,7,11,15]
	
    for i in range(5):
	h_Teff_temp_list = []
	for j in range(4):
		h_ptDisp_dxy_emu_pt_dxy_thresh[i][j].Write()
		h_Teff_temp = R.TEfficiency(h_ptDisp_dxy_emu_pt_dxy_thresh[i][j], h_gen_dxy_den)
		h_Teff_temp.SetLineWidth(2)
		h_Teff_temp.SetName('h_Teff_ptDisp_dxy_pt_'+str(ptthreshlist[i])+'_emudxy_'+str(j))
		h_Teff_temp.SetMarkerStyle(22)
		h_Teff_temp.SetMarkerColor(R.kBlack)
		h_Teff_temp.SetLineColor(j+1)
		h_Teff_temp.Write()
		h_Teff_temp_list.append(h_Teff_temp)
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh.append(h_Teff_temp_list)
		

    h_Teff_ptVtx_dxy_0.Write()
    h_Teff_ptVtx_dxy_4.Write()
    h_Teff_ptVtx_dxy_7.Write()
    h_Teff_ptVtx_dxy_11.Write()
    h_Teff_ptVtx_dxy_15.Write()
    h_Teff_ptDisp_dxy_0.Write()
    h_Teff_ptDisp_dxy_4.Write()
    h_Teff_ptDisp_dxy_7.Write()
    h_Teff_ptDisp_dxy_11.Write()
    h_Teff_ptDisp_dxy_15.Write()
    h_Teff_ptOr_dxy_0.Write()
    h_Teff_ptOr_dxy_4.Write()
    h_Teff_ptOr_dxy_7.Write()
    h_Teff_ptOr_dxy_11.Write()
    h_Teff_ptOr_dxy_15.Write()

####################################

    h_Teff_ptVtx_pt_10.SetName("h_Eff_ptVtx_10_gen_pt")
    h_Teff_ptVtx_pt_10.SetLineWidth(2)
    h_Teff_ptVtx_pt_10.SetLineColor(R.kBlack)
    h_Teff_ptDisp_pt_10.SetName("h_Eff_ptDisp_10_gen_pt")
    h_Teff_ptDisp_pt_10.SetLineWidth(2)
    h_Teff_ptDisp_pt_10.SetLineColor(R.kRed)
    h_Teff_ptOr_pt_10.SetName("h_Eff_ptVtxOrDisp_10_gen_pt")
    h_Teff_ptOr_pt_10.SetLineWidth(2)
    h_Teff_ptOr_pt_10.SetLineColor(R.kBlue)

    h_Teff_ptVtx_pt_10.Write()
    h_Teff_ptDisp_pt_10.Write()
    h_Teff_ptOr_pt_10.Write()

####################################

    h_Teff_ptVtx_dxy_10.SetName("h_Eff_ptVtx_10_gen_dxy")
    h_Teff_ptVtx_dxy_10.SetLineWidth(2)
    h_Teff_ptVtx_dxy_10.SetLineColor(R.kBlack)
    h_Teff_ptDisp_dxy_10.SetName("h_Eff_ptDisp_10_gen_dxy")
    h_Teff_ptDisp_dxy_10.SetLineWidth(2)
    h_Teff_ptDisp_dxy_10.SetLineColor(R.kRed)
    h_Teff_ptOr_dxy_10.SetName("h_Eff_ptVtxOrDisp_10_gen_dxy")
    h_Teff_ptOr_dxy_10.SetLineWidth(2)
    h_Teff_ptOr_dxy_10.SetLineColor(R.kBlue)

    h_Teff_ptVtx_dxy_10.Write()
    h_Teff_ptDisp_dxy_10.Write()
    h_Teff_ptOr_dxy_10.Write()

####################################

    h_Teff_ptVtx_eta_10.SetName("h_Eff_ptVtx_10_gen_eta")
    h_Teff_ptVtx_eta_10.SetLineWidth(2)
    h_Teff_ptVtx_eta_10.SetLineColor(R.kBlack)
    h_Teff_ptDisp_eta_10.SetName("h_Eff_ptDisp_10_gen_eta")
    h_Teff_ptDisp_eta_10.SetLineWidth(2)
    h_Teff_ptDisp_eta_10.SetLineColor(R.kRed)
    h_Teff_ptOr_eta_10.SetName("h_Eff_ptVtxOrDisp_10_gen_eta")
    h_Teff_ptOr_eta_10.SetLineWidth(2)
    h_Teff_ptOr_eta_10.SetLineColor(R.kBlue)

    h_Teff_ptVtx_eta_10.Write()
    h_Teff_ptDisp_eta_10.Write()
    h_Teff_ptOr_eta_10.Write()

####################################

    h_emupt_gendxy.Write()

    #making plots
    c1 = R.TCanvas("c1","Eff",800,600)
    h_eff_gen_dxy.SetStats(0)

    h_Teff_gen_dxy.SetLineColor(R.kBlack)
    h_Teff_gen_dxy.SetMarkerStyle(20)
    h_Teff_gen_dxy.SetMarkerColor(R.kBlack)

    h_Teff_gen_dxy_pt10.SetLineColor(R.kBlue)
    h_Teff_gen_dxy_pt10.SetMarkerStyle(20)
    h_Teff_gen_dxy_pt10.SetMarkerColor(R.kBlue)

    h_Teff_gen_dxy_pt20.SetLineColor(R.kRed)
    h_Teff_gen_dxy_pt20.SetMarkerStyle(20)
    h_Teff_gen_dxy_pt20.SetMarkerColor(R.kRed)

    h_Teff_gen_dxy_pt30.SetLineColor(R.kGreen)
    h_Teff_gen_dxy_pt30.SetMarkerStyle(20)
    h_Teff_gen_dxy_pt30.SetMarkerColor(R.kGreen)

    h_dxy_blank.SetMaximum(1.02)
    h_dxy_blank.SetStats(0)
    h_dxy_blank.GetXaxis().SetTitle("Gen Dxy [cm]")
    h_dxy_blank.GetYaxis().SetTitle("Efficiency")
    h_dxy_blank.Draw()
    h_dxy_blank.Write()

    h_Teff_gen_dxy.Draw("same")
    h_Teff_gen_dxy_pt10.Draw("same")
    h_Teff_gen_dxy_pt20.Draw("same")
    h_Teff_gen_dxy_pt30.Draw("same")

    leg = R.TLegend(0.5,0.6,0.7,0.8)
    leg.AddEntry(h_Teff_gen_dxy, "L1 P_{T} > 0 GeV", "p")
    leg.AddEntry(h_Teff_gen_dxy_pt10, "L1 P_{T} > 10 GeV","p")
    leg.AddEntry(h_Teff_gen_dxy_pt20, "L1 P_{T} > 20 GeV","p")
    leg.AddEntry(h_Teff_gen_dxy_pt30, "L1 P_{T} > 30 GeV","p")
    leg.SetLineWidth(0)
    leg.Draw("same")

    c1.Write()
    c1.Close()

    c2 = R.TCanvas("c2","Eff",800,600)
    h_eff_gen_dxy.SetStats(0)
    h_eff_gen_dxy.GetXaxis().SetTitle("Gen Dxy [cm]")
    h_eff_gen_dxy.GetYaxis().SetTitle("Efficiency")
    h_eff_gen_dxy.GetYaxis().SetRangeUser(0,1.02)

    h_eff_gen_dxy.SetLineColor(R.kBlack)
    h_eff_gen_dxy_pt10.SetLineColor(R.kBlue)
    h_eff_gen_dxy_pt20.SetLineColor(R.kRed)
    h_eff_gen_dxy_pt30.SetLineColor(R.kGreen)

    h_eff_gen_dxy.Draw()
    h_eff_gen_dxy_pt10.Draw("same")
    h_eff_gen_dxy_pt20.Draw("same")
    h_eff_gen_dxy_pt30.Draw("same")

    leg = R.TLegend(0.5,0.6,0.7,0.8)
    leg.AddEntry(h_eff_gen_dxy, "L1 P_{T} > 0 GeV")
    leg.AddEntry(h_eff_gen_dxy_pt10, "L1 P_{T} > 10 GeV")
    leg.AddEntry(h_eff_gen_dxy_pt20, "L1 P_{T} > 20 GeV")
    leg.AddEntry(h_eff_gen_dxy_pt30, "L1 P_{T} > 30 GeV")
    leg.SetLineWidth(0)
    leg.Draw("same")

    c2.Write()
    c2.Close()

#############################

    c_0 = R.TCanvas("c_0","Dxy Eff L1_pt 0",800,600)
    h_dxy_blank.Draw()
    h_Teff_ptVtx_dxy_0.Draw("same")
    h_Teff_ptDisp_dxy_0.Draw("same")
    h_Teff_ptOr_dxy_0.Draw("same")

    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 P_{T} > 0 GeV","C")
    leg.AddEntry(h_Teff_ptVtx_dxy_0, "Constrained")
    leg.AddEntry(h_Teff_ptDisp_dxy_0, "Unconstrained")
    leg.AddEntry(h_Teff_ptOr_dxy_0, "OR of constrained and unconstrained")
    for i in range(4):
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[0][i].Draw("same")
	leg.AddEntry(h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[0][i], "Unconstrained, L1 Dxy >= "+str(i))
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_0.Write()
    c_0.Close()

#######

    c_4 = R.TCanvas("c_4","Dxy Eff L1_pt 4",800,600)

    h_dxy_blank.Draw()
    h_Teff_ptVtx_dxy_4.Draw("same")
    h_Teff_ptDisp_dxy_4.Draw("same")
    h_Teff_ptOr_dxy_4.Draw("same")

    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 P_{T} > 4 GeV","C")
    leg.AddEntry(h_Teff_ptVtx_dxy_4, "Constrained")
    leg.AddEntry(h_Teff_ptDisp_dxy_4, "Unconstrained")
    leg.AddEntry(h_Teff_ptOr_dxy_4, "OR of constrained and unconstrained")
    for i in range(4):
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[1][i].Draw("same")
	leg.AddEntry(h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[1][i], "Unconstrained, L1 Dxy >= "+str(i))
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_4.Write()
    c_4.Close()

#######

    c_7 = R.TCanvas("c_7","Dxy Eff L1_pt 7",800,600)

    h_dxy_blank.Draw()
    h_Teff_ptVtx_dxy_7.Draw("same")
    h_Teff_ptDisp_dxy_7.Draw("same")
    h_Teff_ptOr_dxy_7.Draw("same")

    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 P_{T} > 7 GeV","C")
    leg.AddEntry(h_Teff_ptVtx_dxy_7, "Constrained")
    leg.AddEntry(h_Teff_ptDisp_dxy_7, "Unconstrained")
    leg.AddEntry(h_Teff_ptOr_dxy_7, "OR of constrained and unconstrained")
    for i in range(4):
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[2][i].Draw("same")
	leg.AddEntry(h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[2][i], "Unconstrained, L1 Dxy >= "+str(i))
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_7.Write()
    c_7.Close()

#######

    c_11 = R.TCanvas("c_11","Dxy Eff L1_pt 11",800,600)

    h_dxy_blank.Draw()
    h_Teff_ptVtx_dxy_11.Draw("same")
    h_Teff_ptDisp_dxy_11.Draw("same")
    h_Teff_ptOr_dxy_11.Draw("same")

    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 P_{T} > 11 GeV","C")
    leg.AddEntry(h_Teff_ptVtx_dxy_11, "Constrained")
    leg.AddEntry(h_Teff_ptDisp_dxy_11, "Unconstrained")
    leg.AddEntry(h_Teff_ptOr_dxy_11, "OR of constrained and unconstrained")
    for i in range(4):
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[3][i].Draw("same")
	leg.AddEntry(h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[3][i], "Unconstrained, L1 Dxy >= "+str(i))
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_11.Write()
    c_11.Close()

#######

    c_15 = R.TCanvas("c_15","Dxy Eff L1_pt 15",800,600)

    h_dxy_blank.Draw()
    h_Teff_ptVtx_dxy_15.Draw("same")
    h_Teff_ptDisp_dxy_15.Draw("same")
    h_Teff_ptOr_dxy_15.Draw("same")

    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 P_{T} > 15 GeV","C")
    leg.AddEntry(h_Teff_ptVtx_dxy_15, "Constrained")
    leg.AddEntry(h_Teff_ptDisp_dxy_15, "Unconstrained")
    leg.AddEntry(h_Teff_ptOr_dxy_15, "OR of constrained and unconstrained")
    for i in range(4):
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[4][i].Draw("same")
	leg.AddEntry(h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[4][i], "Unconstrained, L1 Dxy >= "+str(i))
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_15.Write()
    c_15.Close()

###############################################

    ptthreshlist = [0,4,7,11,15]

    c_dxy_0 = R.TCanvas("c_dxy_0","Dxy Eff L1 dxy 0",800,600)

    h_dxy_blank.Draw()
    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 Dxy >= 0","C")
    for i in range(5):
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][0].SetLineColor(i+1)
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][0].Draw("same")
	leg.AddEntry(h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][0], "Unconstrained L1 P_{T} > "+str(ptthreshlist[i]))
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_dxy_0.Write()
    c_dxy_0.Close()

#########

    c_dxy_1 = R.TCanvas("c_dxy_1","Dxy Eff L1 dxy 1",800,600)

    h_dxy_blank.Draw()
    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 Dxy >= 1","C")
    for i in range(5):
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][1].SetLineColor(i+1)
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][1].Draw("same")
	leg.AddEntry(h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][1], "Unconstrained L1 P_{T} > "+str(ptthreshlist[i]))
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_dxy_1.Write()
    c_dxy_1.Close()

#########

    c_dxy_2 = R.TCanvas("c_dxy_2","Dxy Eff L1 dxy 2",800,600)

    h_dxy_blank.Draw()
    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 Dxy >= 2","C")
    for i in range(5):
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][2].SetLineColor(i+1)
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][2].Draw("same")
	leg.AddEntry(h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][2], "Unconstrained L1 P_{T} > "+str(ptthreshlist[i]))
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_dxy_2.Write()
    c_dxy_2.Close()

#########

    c_dxy_3 = R.TCanvas("c_dxy_3","Dxy Eff L1 dxy 3",800,600)

    h_dxy_blank.Draw()
    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 Dxy >= 3","C")
    for i in range(5):
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][3].SetLineColor(i+1)
	h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][3].Draw("same")
	leg.AddEntry(h_Teff_ptDisp_dxy_emu_pt_dxy_thresh[i][3], "Unconstrained L1 P_{T} > "+str(ptthreshlist[i]))
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_dxy_3.Write()
    c_dxy_3.Close()

################################################

    c_pt_10 = R.TCanvas("c_pt_10","Pt Eff L1_pt 10",800,600)

    h_pt_blank.SetMaximum(1.02)
    h_pt_blank.SetStats(0)
    h_pt_blank.GetXaxis().SetTitle("Gen Pt [GeV]")
    h_pt_blank.GetXaxis().SetRangeUser(0, 65.)
    h_pt_blank.GetYaxis().SetTitle("Efficiency")
    h_pt_blank.Write()
    h_pt_blank.Draw()

    h_Teff_ptVtx_pt_10.Draw("same")
    h_Teff_ptDisp_pt_10.Draw("same")
    h_Teff_ptOr_pt_10.Draw("same")

    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("L1 P_{T} > 10 GeV","C")
    leg.AddEntry(h_Teff_ptVtx_pt_10, "Constrained")
    leg.AddEntry(h_Teff_ptDisp_pt_10, "Unconstrained")
    leg.AddEntry(h_Teff_ptOr_pt_10, "OR of constrained and unconstrained")
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_pt_10.Write()
    c_pt_10.Close()

################################################

    c_dxy_10 = R.TCanvas("c_dxy_10","Dxy Eff L1_pt 10",800,600)

    h_dxy_blank.SetMaximum(1.02)
    h_dxy_blank.SetStats(0)
    h_dxy_blank.GetXaxis().SetTitle("Gen Dxy [cm]")
    h_dxy_blank.GetXaxis().SetRangeUser(0, 300.)
    h_dxy_blank.GetYaxis().SetTitle("Efficiency")
    h_dxy_blank.Write()
    h_dxy_blank.Draw()

    h_Teff_ptVtx_dxy_10.Draw("same")
    h_Teff_ptDisp_dxy_10.Draw("same")
    h_Teff_ptOr_dxy_10.Draw("same")

    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("Gen Dxy for L1 P_{T} > 10 GeV","C")
    leg.AddEntry(h_Teff_ptVtx_dxy_10, "Constrained")
    leg.AddEntry(h_Teff_ptDisp_dxy_10, "Unconstrained")
    leg.AddEntry(h_Teff_ptOr_dxy_10, "OR of constrained and unconstrained")
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_dxy_10.Write()
    c_dxy_10.Close()

################################################

    c_eta_10 = R.TCanvas("c_eta_10","Eta Eff L1_pt 10",800,600)

    h_eta_blank.SetMaximum(1.02)
    h_eta_blank.SetStats(0)
    h_eta_blank.GetXaxis().SetTitle("Gen Eta")
    h_eta_blank.GetXaxis().SetRangeUser(-3., 3.)
    h_eta_blank.GetYaxis().SetTitle("Efficiency")
    h_eta_blank.Write()
    h_eta_blank.Draw()

    h_Teff_ptVtx_eta_10.Draw("same")
    h_Teff_ptDisp_eta_10.Draw("same")
    h_Teff_ptOr_eta_10.Draw("same")

    leg = R.TLegend(0.2,0.4,0.8,0.9)
    leg.SetHeader("Gen Eta for L1 P_{T} > 10 GeV","C")
    leg.AddEntry(h_Teff_ptVtx_eta_10, "Constrained")
    leg.AddEntry(h_Teff_ptDisp_eta_10, "Unconstrained")
    leg.AddEntry(h_Teff_ptOr_eta_10, "OR of constrained and unconstrained")
    leg.SetLineWidth(0)
    leg.Draw("same")

    c_eta_10.Write()
    c_eta_10.Close()

################################################

    h_gen_Vx.Write()
    h_gen_Vy.Write()
    h_gen_Vz.Write()
    h_gen_eta.Write()
    h_gen_phi.Write()
#    h_gen_ptinv.Write()

############### Dimuon plots ##################

#    h_4mu_lxy.Write()

    h_dimu_dxy1_dxy2.Write()
    h_dimu_pt1_pt2.Write()
    h_dimu_dxy1_pt2.Write()
    h_dimu_dxy2_pt1.Write()
    h_dimu_dxy1_pt1.Write()
    h_dimu_dxy2_pt2.Write()

# Efficiencies

    h_Teff_dimu_dxy_dxy = []

    for i in range(len(L1ptthresh)):
	h_Teff_templist = []
	for j in range(len(L1ptthresh)):
		if (R.TEfficiency.CheckConsistency(h_num_dxy1_dxy2[i][j], h_den_dxy1_dxy2[i][j])):
			h_Teff_temp = R.TEfficiency(h_num_dxy1_dxy2[i][j], h_den_dxy1_dxy2[i][j])
		else :	h_Teff_temp = R.TEfficiency(h_blank_dxy1_dxy2, h_den_dxy1_dxy2[i][j])
		h_Teff_templist.append(h_Teff_temp)
    	h_Teff_dimu_dxy_dxy.append(h_Teff_templist)

    for i in range(1, len(L1ptthresh)):
	for j in range(i):
		h_num_dxy1_dxy2[i][j].Write()
		h_den_dxy1_dxy2[i][j].Write()
		h_Teff_dimu_dxy_dxy[i][j].Write()

################################################

    out_file.Close()
    del chains

    print '\nWrote out file: plots/'+out_file_str+'.root'
    print '\n Events run over: ', iEvt
    print 'Number of L1_SingleMu20_BMTF: ', L1_SingleMu20_BMTF, '\t rate: ', float(L1_SingleMu20_BMTF)*scale[evtclassid]/iEvt

    print 'Event counter'
    print 'Total GenMus: ', GenMuCount
    print 'GenMus in Acceptance: ', GenMuinAcc
    print 'EmuMus passing qual: ', EmuMus_count
    print 'EmuMus passing qual and unique: ', EmuMus_unique_count
    print 'Events with non-unique EmuMus: ', mucopy_count
    
    print 'Mismatch between corr and uncorr Gen Muons matched to qual emus: ', mismatch_evt

    print 'Matched gen mus : ', matched_gen_mu_count
    print 'Events with matched gen Mus: ', matched_gen_mu_evt

    print 'Events with Gen Dimuons in acceptance: ', GenDimu_Evt

if __name__ == '__main__':
    main()



#                h_pt_displ_emu.Fill( min( max( ptDispl+0.01, pt_bins[1]+0.01), pt_bins[2]-0.01) )
#    h_pt_vtx_unp   = R.TH1F('h_pt_vtx_unp',   'Legacy BMTF vertex-constrained pT spectrum',         pt_bins[0], pt_bins[1], pt_bins[2])
#    h_pt_vtx_emu   = R.TH1F('h_pt_vtx_emu',   'Kalman BMTF vertex-constrained pT spectrum',         pt_bins[0], pt_bins[1], pt_bins[2])
#    h_phi_mu_emu   = R.TH1F('h_phi_mu_emu',   'Kalman BMTF muon Phi',         		phi_bins[0], phi_bins[1], phi_bins[2])
#    h_pt_displ_emu = R.TH1F('h_pt_displ_emu', 'Kalman BMTF non-vertex-constrained pT spectrum',     pt_bins[0], pt_bins[1], pt_bins[2])
#    h_pt_dxy_emu   = R.TH2F('h_pt_dxy_emu',   'Kalman BMTF non-vertex-constrained pT vs dxy',       pt_bins[0], pt_bins[1], pt_bins[2], dxy_bins[0], dxy_bins[1], dxy_bins[2])
#    h_pt_displ_kmt = R.TH1F('h_pt_displ_kmt', 'Internal Kalman non-vertex-constrained pT spectrum', pt_bins[0], pt_bins[1], pt_bins[2])
#    h_chi2_kmt = R.TH1F('h_chi2_kmt', 'Kalman BMTF track #chi^{2} distribution', chi_bins[0], chi_bins[1], chi_bins[2])

 #   h_Vx_gen = R.TH1F('h_Vx_gen', 'Generator level Vertex x coordinate',     Vx_bins[0], Vx_bins[1], Vx_bins[2])

#    h_pt_vtx_unp.SetLineWidth(2)
#    h_pt_vtx_unp.SetLineColor(R.kBlack)
#    h_pt_vtx_unp.Write()
    
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

#    c2 = R.TCanvas("c2","Eff",800,600)
#    h_eff_gen_pt.SetStats(0)
#    h_eff_gen_pt.GetXaxis().SetTitle("Gen P_{T} [GeV]")
#    h_eff_gen_pt.GetYaxis().SetTitle("Efficiency")
#    h_eff_gen_pt.GetYaxis().SetRangeUser(0,1.02)

#    h_eff_gen_pt_pt10.SetLineColor(R.kBlue)
#    h_eff_gen_pt_pt20.SetLineColor(R.kRed)

#    h_eff_gen_pt.Draw()
#    h_eff_gen_pt_pt10.Draw("same")
#    h_eff_gen_pt_pt20.Draw("same")

#    leg = R.TLegend(0.5,0.6,0.7,0.8)
#    leg.AddEntry("h_eff_gen_pt", "L1 P_{T} > 0 GeV")
#    leg.AddEntry("h_eff_gen_pt_pt10", "L1 P_{T} > 10 GeV")
#    leg.AddEntry("h_eff_gen_pt_pt20", "L1 P_{T} > 20 GeV")
#    leg.SetLineWidth(0)
#    leg.Draw("same")


#	    with open('EtaPhiDebug.txt', 'a') as f:
#		for i in range(len(dR_list)):
#			Gen_idx = dR_list[i][1]
#			Emu_idx = dR_list[i][2]
#
#			Lxy = sqrt(float(Gen_br.partVx[Gen_idx])**2 + float(Gen_br.partVy[Gen_idx])**2)
#			if Lxy < 30. : continue
#
#			Gen_Eta = float(Gen_br.partEta[Gen_idx])
#			Gen_Eta_prop = getGenEtaBarrel(Gen_br, Gen_idx)
##			Gen_Phi = float(Gen_br.partPhi[Gen_idx])
#			Gen_Phi_prop = getGenPhiBarrel(Gen_br, Gen_idx)	
#			

#			Emu_Eta = float(uGT_br.muonEta[Emu_idx])
#			Emu_Eta_Vtx = float(uGT_br.muonEtaAtVtx[Emu_idx])
#			Emu_Phi = float(uGT_br.muonPhi[Emu_idx])
#			Emu_Phi_Vtx = float(uGT_br.muonPhiAtVtx[Emu_idx])

#			f.write('\n Dxy : %.4f' %(Lxy))
#			f.write('\nGen: EtaVtx, PhiVtx : \t %.4f %.4f' %(Gen_Eta, Gen_Phi))
#			f.write('\nGen: Eta, Phi at MS2 : \t %.4f %.4f' %(Gen_Eta_prop, Gen_Phi_prop))
#			f.write('\n')
#			f.write('\nEmu: EtaVtx, PhiVtx : \t %.4f %.4f' %(Emu_Eta_Vtx, Emu_Phi_Vtx))
#			f.write('\nEmu: Eta, Phi at MS2 : \t %.4f %.4f' %(Emu_Eta, Emu_Phi))
#			f.write('\n--------------------------------------\n')

	    #Plotting dR between Gen Muon and corrected Gen Muon.

#	    for i in GenMus_acc:
#		phi1 = float(Gen_br.partPhi[i])
#		eta1 = float(Gen_br.partEta[i])
#		phi2 = 0 #getGenPhiBarrel(Gen_br, i)
#		eta2 = 0 #getGenEtaBarrel(Gen_br, i)
#		pt = float(Gen_br.partPt[i])

#		dR = deltaR(phi1, phi2, eta1, eta2)
#		dphi = dPhi(phi1, phi2)

#		y = float(Gen_br.partVy[i])
#		x = float(Gen_br.partVx[i])
#
#		alpha = atan(y/x)
#		if y > 0. and x < 0. : alpha = alpha + pi
#		elif y < 0. and x < 0. : alpha = alpha - pi
#		dAlphi = dPhi(alpha, phi2)
	
#		h_gen_dphi_pt.Fill(dphi, pt)
#		h_gen_dalphi_pt.Fill(dAlphi, pt)
#		h_gen_phi_phit.Fill(phi1, phi2)
#		h_gen_alpha_phi.Fill(phi1, alpha)

#		h_gen_dR.Fill(dR)
#		Lxy = sqrt(Gen_br.partVx[i]**2+Gen_br.partVy[i]**2)
		
#		h_gen_dphi_lxy.Fill(dphi, Lxy)
#		h_gen_dalphi_lxy.Fill(dAlphi, Lxy)

#		deta = abs(eta1 - eta2)
#		h_gen_deta_lxy.Fill(deta, Lxy)
#
#		h_gen_dR_Lxy.Fill(Lxy, dR)

#		muHw = R.TLorentzVector()
#		muHw.SetPtEtaPhiM(float(uGT_br.muonIEtUnconstrained[i]-1)*0.5, float(uGT_br.muonIEta[i])*0.010875, (float(uGT_br.muonIPhi[i]))*0.010908, 105.7e-3)
#		pt_phy = mu1.Pt()		
#		pt_Hw = muHw.Pt()
#		ptvtx_phy = float(uGT_br.muonEt[i])
#		ptvtx_Hw = float(uGT_br.muonIEt[i]-1)*0.5
#		eta_phy = mu1.Eta()
#		eta_Hw =  muHw.Eta()
#		phi_phy = mu1.Phi()
#		phi_Hw = muHw.Phi()
		
#		h_uGTMu_pt_pt.Fill(pt_Hw, pt_phy)
#		h_uGTMu_ptvtx_ptvtx.Fill(ptvtx_Hw, ptvtx_phy)
#		h_uGTMu_eta_eta.Fill(eta_Hw, eta_phy)
#		h_uGTMu_phi_phi.Fill(phi_Hw, phi_phy)

#		h_uGTMu_pt_diff.Fill(pt_phy, pt_phy-pt_Hw)
#		h_uGTMu_ptvtx_diff.Fill(ptvtx_phy, ptvtx_phy-ptvtx_Hw)
#		h_uGTMu_phi_diff.Fill(phi_phy, phi_phy-phi_Hw)
#		h_uGTMu_eta_diff.Fill(eta_phy, eta_phy-eta_Hw)


#    h_ptVtx_debug_denom		= R.TH1F('h_ptVtx_debug_denom',		'Gen Pt',		pt_bins[0], pt_bins[1], pt_bins[2])
#    h_ptIVtx_debug_denom	= R.TH1F('h_ptIVtx_debug_denom',	'Gen Pt',		pt_bins[0], pt_bins[1], pt_bins[2])
#    h_ptDisp_debug_denom	= R.TH1F('h_ptDisp_debug_denom',	'Gen Pt',		pt_bins[0], pt_bins[1], pt_bins[2])
#    h_ptIDisp_debug_denom	= R.TH1F('h_ptIDisp_debug_denom',	'Gen Pt',		pt_bins[0], pt_bins[1], pt_bins[2])
#    h_genpt_debug_denom		= R.TH1F('h_genpt_debug_denom',		'Gen Pt',		pt_bins[0], pt_bins[1], pt_bins[2])

#    h_ptVtx_genpt_rel		= R.TH2F('h_ptVtx_genpt_rel', 		'GenPt',		pt_bins[0], pt_bins[1], pt_bins[2], 50, -5., 5.)
#    h_ptIVtx_genpt_rel		= R.TH2F('h_ptIVtx_genpt_rel', 		'GenPt',		pt_bins[0], pt_bins[1], pt_bins[2], 50, -5., 5.)
#    h_ptDisp_genpt_rel		= R.TH2F('h_ptDisp_genpt_rel', 		'GenPt',		pt_bins[0], pt_bins[1], pt_bins[2], 50, -5., 5.)
#    h_ptIDisp_genpt_rel		= R.TH2F('h_ptIDisp_genpt_rel', 	'GenPt',		pt_bins[0], pt_bins[1], pt_bins[2], 50, -5., 5.)

#    h_ptVtx_genpt_dxy		= R.TH2F('h_ptVtx_genpt_dxy', 			'Gen Dxy [GeV]',		dxy_bins[0], dxy_bins[1], dxy_bins[2], 100, -2., 2.)
#    h_ptIVtx_genpt_dxy		= R.TH2F('h_ptIVtx_genpt_dxy',	 		'Gen Dxy [GeV]',		dxy_bins[0], dxy_bins[1], dxy_bins[2], 100, -2., 2.)
#    h_ptDisp_genpt_dxy		= R.TH2F('h_ptDisp_genpt_dxy', 			'Gen Dxy [GeV]',		dxy_bins[0], dxy_bins[1], dxy_bins[2], 100, -2., 2.)
#    h_ptIDisp_genpt_dxy		= R.TH2F('h_ptIDisp_genpt_dxy', 		'Gen Dxy [GeV]',		dxy_bins[0], dxy_bins[1], dxy_bins[2], 100, -2., 2.)

#    h_ptVtx_genpt_qrelpt	= R.TH2F('h_ptVtx_genpt_qrelpt', 		'Gen Pt',		pt_bins[0], pt_bins[1], pt_bins[2], 100, -2., 2.)
#    h_ptIVtx_genpt_qrelpt	= R.TH2F('h_ptIVtx_genpt_qrelpt', 		'Gen Pt',		pt_bins[0], pt_bins[1], pt_bins[2], 100, -2., 2.)
#    h_ptDisp_genpt_qrelpt	= R.TH2F('h_ptDisp_genpt_qrelpt', 		'Gen Pt',		pt_bins[0], pt_bins[1], pt_bins[2], 100, -2., 2.)
#    h_ptIDisp_genpt_qrelpt	= R.TH2F('h_ptIDisp_genpt_qrelpt', 		'Gen Pt',		pt_bins[0], pt_bins[1], pt_bins[2], 100, -2., 2.)


#    h_uGTMu_pt_pt.GetXaxis().SetTitle("Hardware pT [GeV]")
#    h_uGTMu_pt_pt.GetYaxis().SetTitle("Physical pT [GeV]")
#    h_uGTMu_pt_pt.Write()

#    h_uGTMu_phi_phi.GetXaxis().SetTitle("Hardware phi")
#    h_uGTMu_phi_phi.GetYaxis().SetTitle("Physical phi")
#    h_uGTMu_phi_phi.Write()

#    h_uGTMu_eta_eta.GetXaxis().SetTitle("Hardware eta")
#    h_uGTMu_eta_eta.GetYaxis().SetTitle("Physical eta")
#    h_uGTMu_eta_eta.Write()

#    h_uGTMu_ptvtx_ptvtx.GetXaxis().SetTitle("Hardware pT [GeV]")
#    h_uGTMu_ptvtx_ptvtx.GetYaxis().SetTitle("Physical pT [GeV]")
#    h_uGTMu_ptvtx_ptvtx.Write()

######

#    h_uGTMu_pt_diff.GetYaxis().SetTitle("Physical pT [GeV]")
#    h_uGTMu_pt_diff.Write()

#    h_uGTMu_ptvtx_diff.GetYaxis().SetTitle("Physical pT [GeV]")
#    h_uGTMu_ptvtx_diff.Write()

#    h_uGTMu_phi_diff.GetYaxis().SetTitle("Physical phi")
#    h_uGTMu_phi_diff.Write()

#    h_uGTMu_eta_diff.GetYaxis().SetTitle("Physical eta")
#    h_uGTMu_eta_diff.Write()
    
###############################################

#Debug hists

#    h_ptVtx_debug_denom.SetLineWidth(2)
#    h_ptVtx_debug_denom.SetLineColor(R.kBlue)
#    h_ptVtx_debug_denom.Write()

#    h_ptIVtx_debug_denom.SetLineWidth(2)
#    h_ptIVtx_debug_denom.SetLineColor(R.kGreen)
#    h_ptIVtx_debug_denom.Write()

#    h_ptDisp_debug_denom.SetLineWidth(2)
#    h_ptDisp_debug_denom.SetLineColor(R.kRed)
#    h_ptDisp_debug_denom.Write()

#    h_ptIDisp_debug_denom.SetLineWidth(2)
#    h_ptIDisp_debug_denom.SetLineColor(6)
#    h_ptIDisp_debug_denom.Write()

#    h_genpt_debug_denom.SetLineWidth(2)
#    h_genpt_debug_denom.SetLineColor(R.kBlack)
#    h_genpt_debug_denom.Write()

#    h_uGTMu_ptVtx_genpt.Write()
#    h_uGTMu_ptIVtx_genpt.Write()
#    h_uGTMu_ptDisp_genpt.Write()
#    h_uGTMu_ptIDisp_genpt.Write()

#    h_ptVtx_genpt_rel.Write()
#    h_ptIVtx_genpt_rel.Write()
#    h_ptDisp_genpt_rel.Write()
#    h_ptIDisp_genpt_rel.Write()

##################

#    h_ptVtx_genpt_dxy.GetXaxis().SetTitle("Gen Dxy [cm]")
#    h_ptVtx_genpt_dxy.GetYaxis().SetTitle("#frac{q/p_{T}(L1) - q/p_{T}(Gen)}{q/p_{T}(Gen)} (Constrained, Physical)")
#    h_ptVtx_genpt_dxy.Write()

#    h_ptIVtx_genpt_dxy.GetXaxis().SetTitle("Gen Dxy [cm]")
#    h_ptIVtx_genpt_dxy.GetYaxis().SetTitle("#frac{q/p_{T}(L1) - q/p_{T}(Gen)}{q/p_{T}(Gen)} (Constrained, Hardware)")
#    h_ptIVtx_genpt_dxy.Write()

#    h_ptDisp_genpt_dxy.GetXaxis().SetTitle("Gen Dxy [cm]")
#    h_ptDisp_genpt_dxy.GetYaxis().SetTitle("#frac{q/p_{T}(L1) - q/p_{T}(Gen)}{q/p_{T}(Gen)} (Unconstrained, Physical)")
#    h_ptDisp_genpt_dxy.Write()
		
#    h_ptIDisp_genpt_dxy.GetXaxis().SetTitle("Gen Dxy [cm]")
#    h_ptIDisp_genpt_dxy.GetYaxis().SetTitle("#frac{q/p_{T}(L1) - q/p_{T}(Gen)}{q/p_{T}(Gen)} (Unconstrained, Hardware)")
#    h_ptIDisp_genpt_dxy.Write()

###################

#    h_ptVtx_genpt_qrelpt.GetXaxis().SetTitle("Gen p_{T} [GeV]")
#    h_ptVtx_genpt_qrelpt.GetYaxis().SetTitle("#frac{q/p_{T}(L1) - q/p_{T}(Gen)}{q/p_{T}(Gen)} (Constrained, Physical)")
#    h_ptVtx_genpt_qrelpt.Write()

#    h_ptIVtx_genpt_qrelpt.GetXaxis().SetTitle("Gen p_{T} [GeV]")
#    h_ptIVtx_genpt_qrelpt.GetYaxis().SetTitle("#frac{q/p_{T}(L1) - q/p_{T}(Gen)}{q/p_{T}(Gen)} (Constrained, Hardware)")
#    h_ptIVtx_genpt_qrelpt.Write()

#    h_ptDisp_genpt_qrelpt.GetXaxis().SetTitle("Gen p_{T} [GeV]")
#    h_ptDisp_genpt_qrelpt.GetYaxis().SetTitle("#frac{q/p_{T}(L1) - q/p_{T}(Gen)}{q/p_{T}(Gen)} (Unconstrained, Physical)")
#    h_ptDisp_genpt_qrelpt.Write()

#    h_ptIDisp_genpt_qrelpt.GetXaxis().SetTitle("Gen p_{T} [GeV]")
#    h_ptIDisp_genpt_qrelpt.GetYaxis().SetTitle("#frac{q/p_{T}(L1) - q/p_{T}(Gen)}{q/p_{T}(Gen)} (Constrained, Hardware)")
#    h_ptIDisp_genpt_qrelpt.Write()

##############################################

