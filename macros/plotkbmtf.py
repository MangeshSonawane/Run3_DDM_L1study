from ROOT import *
from math import *
import os, sys, time, getopt, glob

DRCUT = 0.5
L1PTCUT = 10.

flist = glob.glob("data/*.root")

t_gen = TChain("l1GeneratorTree/L1GenTree")
t_unp = TChain("l1UpgradeTfMuonTree/L1UpgradeTfMuonTree")
t_emu = TChain("l1UpgradeTfMuonEmuTree/L1UpgradeTfMuonTree")

for fname in flist:
    print fname
    t_gen.Add(fname)
    t_unp.Add(fname)
    t_emu.Add(fname)
#    break

bins_Dxy = [30.,60.,90.,300.]
string_Dxy = ["Dxy < 30 cm", "30 < Dxy < 60 cm", "60 < Dxy < 90 cm","90 < Dxy < 300 cm"]
colors_unp_Dxy = [1,4,2,3]
colors_emu_Dxy = [1,4,2,3]

h_tot = []
h_unp_eff = []
h_emu_eff = []
h_pt_gen = []
h_pt_unp = []
h_pt_emu = []
h_pt_emud = []
h_dpt_unp = []
h_dpt_emud = []
h_dptx_unp = []
h_dptx_emud = []
h_dxy = []
h_dxy_gen20 = []
h_dxy_eff10 = []
for ib, b in enumerate(bins_Dxy):
    h_tot.append(TH1F("h_tot"+str(ib),"",50,0.,100))
    h_unp_eff.append(TH1F("h_unp_eff"+str(ib),"",50,0.,100))
    h_emu_eff.append(TH1F("h_emu_eff"+str(ib),"",50,0.,100))
    h_pt_gen.append(TH1F("h_pt_gen"+str(ib),"",50,0,200))
    h_pt_unp.append(TH1F("h_pt_unp"+str(ib),"",50,0,200))
    h_pt_emu.append(TH1F("h_pt_emu"+str(ib),"",50,0,200))
    h_pt_emud.append(TH1F("h_pt_emud"+str(ib),"",50,0,200))
    h_dpt_unp.append(TH1F("h_dpt_unp"+str(ib),"",100,-100,100))
    h_dpt_emud.append(TH1F("h_dpt_emud"+str(ib),"",100,-100,100))
    h_dptx_unp.append(TH1F("h_dptx_unp"+str(ib),"",100,-100,100))
    h_dptx_emud.append(TH1F("h_dptx_emud"+str(ib),"",100,-100,100))
    h_dxy.append(TH1F("h_dxy"+str(ib),"",100,-100,100))
    h_dxy_gen20.append(TH1F("h_dxy_gen20"+str(ib),"",100,-100,100))
    h_dxy_eff10.append(TH1F("h_dxy_eff10"+str(ib),"",100,-100,100))

h2_tot = TH1F("h2_tot","",30,0.,150)
h2_unp_eff = TH1F("h2_unp_eff","",30,0.,150)
h2_emu_eff = TH1F("h2_emu_eff","",30,0.,150)
h2_unpx_eff = TH1F("h2_unpx_eff","",30,0.,150)
h2_emux_eff = TH1F("h2_emux_eff","",30,0.,150)
h2_unpx2_eff = TH1F("h2_unpx2_eff","",30,0.,150)
h2_emux2_eff = TH1F("h2_emux2_eff","",30,0.,150)
h2_emu_l1dxy15_eff = TH1F("h2_emux_l1dxy15_eff","",30,0.,150)
h2_emu_l1dxy45_eff = TH1F("h2_emux_l1dxy45_eff","",30,0.,150)
h2_emu_l1dxy75_eff = TH1F("h2_emux_l1dxy75_eff","",30,0.,150)

hh1 = TH1F("hh1","",50,0.,200)
hh2 = TH1F("hh2","",50,0.,200)
hh3 = TH1F("hh3","",50,0.,200)
hh4 = TH1F("hh4","",50,0.,200)

#Evt_br = L1Analysis.L1AnalysisEventDataFormat()
Gen_br = L1Analysis.L1AnalysisGeneratorDataFormat()
Unp_br = L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
Emu_br = L1Analysis.L1AnalysisL1UpgradeTfMuonDataFormat()
Kmt_br = L1Analysis.L1AnalysisBMTFOutputDataFormat()

t_gen.SetBranchAddress("Generator",AddressOf(Gen_br))
t_unp.SetBranchAddress("L1UpgradeBmtfMuon",AddressOf(Unp_br))
t_emu.SetBranchAddress("L1UpgradeBmtfMuon",AddressOf(Emu_br))
t_emu.SetBranchAddress("L1UpgradeBmtfOutput",AddressOf(Kmt_br))

def sign(x): return 1 if x >= 0 else -1

def impactParameter(vx,vy,pt,phi,ch):
    r = 88. * pt
    cx = vx + ch * r * sin(phi)
    cy = vy - ch * r * cos(phi)
    ip = abs(r - sqrt(pow(cx,2)+pow(cy,2)))
    return ip

def impactParameterSimple(vx,vy,pt,phi,ch):
    ip = abs(vx*sin(phi)-vy*cos(phi))
    return ip

def inacceptance(vx,vy,vz,eta):
    Lxy = sqrt(pow(vx,2)+pow(vy,2))
    Lz = vz
#    if Lxy >= 350.: return False
#    if abs(Lz) >= 300.: return False
    if Lxy >= 700.: return False
    if abs(Lz) >= 650.: return False
    maxeta = -log(tan(0.5*atan((700.-Lxy)/(650.-Lz))))
    mineta = -log(tan(0.5*(pi-atan((700.-Lxy)/(650.+Lz)))))
    if eta < mineta: return False
    if eta > maxeta: return False
    return True

def phiL1(globalPhiHw):
    phi = globalPhiHw/287.5*pi if globalPhiHw<287.5 else \
                 (globalPhiHw-575.)/287.5*pi
    return phi
    
def dPhi(phi1,phi2):
    return acos(cos(phi1-phi2))

def getdrmin(brg,brm):
    drmin = 10.
    drmin_im = -1
    for im in range(brm.nTfMuons):
        if brm.tfMuonBx[im] != 0: continue
        phi_mu = phiL1(float(brm.tfMuonGlobalPhi[im]))
        eta_mu = float(brm.tfMuonHwEta[im])*0.010875
        dphi = dPhi(phi_mu,float(brg.partPhi[ip]))
        deta = eta_mu - float(brg.partEta[ip])
        dr = sqrt(pow(dphi,2)+pow(deta,2))
        if dr < drmin:
            drmin = dr
            drmin_im = im
    return drmin_im, drmin

def addoverflow(H):
    nb = H.GetNbinsX()
    x1 = H.GetBinContent(nb)
    x2 = H.GetBinContent(nb+1)
    x3 = x1+x2
    H.SetBinContent(nb,x3)
    H.SetBinContent(nb+1,0.)

def normalise(H):
    tot = H.Integral()
    if tot == 0: tot = 1e-6
    H.Scale(1./tot)

for i in range(t_gen.GetEntries()):
#    if i>20: break
    t_gen.GetEntry(i)
    t_unp.GetEntry(i)
    t_emu.GetEntry(i)
    
    for ip in range(Gen_br.nPart):
        if abs(Gen_br.partId[ip]) != 13 or Gen_br.partParent[ip] != 6000113: continue
#        if abs(Gen_br.partEta[ip]) > 0.9: continue
        if not inacceptance(Gen_br.partVx[ip],Gen_br.partVy[ip],Gen_br.partVz[ip],Gen_br.partEta[ip]): continue
        
        ptGen = Gen_br.partPt[ip]
        chGen = -sign(Gen_br.partId[ip])
        Dxy = impactParameterSimple(Gen_br.partVx[ip],Gen_br.partVy[ip],ptGen,Gen_br.partPhi[ip],chGen)
        
        drmin_unp_im, drmin_unp = getdrmin(Gen_br, Unp_br)
        drmin_emu_im, drmin_emu = getdrmin(Gen_br, Emu_br)
                
        ptL1u = -999.
        if drmin_unp_im > -1:
            ptL1u = (float(Unp_br.tfMuonHwPt[drmin_unp_im])-1.)*0.5
        ptL1 = -999.
        ptL1d = -999.
        if drmin_emu_im > -1:
            ptL1 = (float(Emu_br.tfMuonHwPt[drmin_emu_im])-1.)*0.5
            ptL1d = float(Emu_br.tfMuonHwPtDispl[drmin_emu_im])
            ptL1d = max(ptL1d,ptL1)
            dxy_emu = abs(float(Kmt_br.dxy[drmin_emu_im])/10.)*1.5
            
        hh1.Fill(ptGen)
        if drmin_unp < DRCUT:  hh2.Fill(ptL1u) 
        if drmin_emu < DRCUT:
            hh3.Fill(ptL1)
            hh4.Fill(ptL1d)
            
        if ptGen > 20.:
            h2_tot.Fill(Dxy)
            if drmin_unp < DRCUT and ptL1u > 15.: h2_unpx2_eff.Fill(Dxy)
            if drmin_emu < DRCUT and ptL1d > 15.: h2_emux2_eff.Fill(Dxy)
            if drmin_unp < DRCUT and ptL1u > L1PTCUT: h2_unp_eff.Fill(Dxy)
            if drmin_emu < DRCUT and ptL1d > L1PTCUT: h2_emu_eff.Fill(Dxy)
            if drmin_unp < DRCUT: h2_unpx_eff.Fill(Dxy)
            if drmin_emu < DRCUT: h2_emux_eff.Fill(Dxy)
            if drmin_emu < DRCUT and ptL1d > L1PTCUT:
                if dxy_emu > 15.: h2_emu_l1dxy15_eff.Fill(Dxy)
                if dxy_emu > 45.: h2_emu_l1dxy45_eff.Fill(Dxy)
                if dxy_emu > 75.: h2_emu_l1dxy75_eff.Fill(Dxy)
        for ib, b in enumerate(bins_Dxy):
            b1 = 0.
            if ib > 0: b1 = bins_Dxy[ib-1]
            if Dxy < b and Dxy >= b1:
                h_tot[ib].Fill(ptGen)
                h_pt_gen[ib].Fill(ptGen)
                if drmin_unp < DRCUT:
                    h_pt_unp[ib].Fill(ptL1u)
                    if ptGen > 20. and ptGen < 30.: 
#                    if ptGen > 20.: 
                        h_dpt_unp[ib].Fill(ptL1u-ptGen)
                        if ptL1u > L1PTCUT: h_dptx_unp[ib].Fill(ptL1u-ptGen)
                    if ptL1u > L1PTCUT: h_unp_eff[ib].Fill(ptGen)
                if drmin_emu < DRCUT:
                    h_pt_emu[ib].Fill(ptL1)
                    h_pt_emud[ib].Fill(ptL1d)
                    if ptGen > 20. and ptGen < 30.:
#                    if ptGen > 20.:
                        h_dpt_emud[ib].Fill(ptL1d-ptGen)
                        if ptL1d > L1PTCUT: h_dptx_emud[ib].Fill(ptL1d-ptGen)
                    if ptL1d > L1PTCUT: h_emu_eff[ib].Fill(ptGen)
                    h_dxy[ib].Fill(dxy_emu-Dxy)
                    if ptGen > 20.:
                        h_dxy_gen20[ib].Fill(dxy_emu-Dxy)
                        if ptL1d > L1PTCUT:
                            h_dxy_eff10[ib].Fill(dxy_emu-Dxy)

              

gStyle.SetOptStat(kFALSE)

'''

c1 = TCanvas("c1")
c1.SetGrid(1,1)
legc1 = TLegend(0.1,0.9,0.9,1.)
legc1.SetNColumns(5)
dummy = TH1F("dummy","",10,0,100)
dummy.SetTitle(";pTgen;Efficiency")
dummy.SetMaximum(1.)
dummy.Draw()
teffu = []
for ib, b in enumerate(bins_Dxy):
    teffu.append(TEfficiency(h_unp_eff[ib],h_tot[ib]))
    teffu[ib].SetLineColor(colors_unp_Dxy[ib])
    teffu[ib].SetMarkerStyle(24)
    teffu[ib].SetMarkerColor(colors_unp_Dxy[ib])
    opt = "PZ same"
    teffu[ib].Draw(opt)

teff = []
for ib, b in enumerate(bins_Dxy):
    teff.append(TEfficiency(h_emu_eff[ib],h_tot[ib]))
    teff[ib].SetLineColor(colors_emu_Dxy[ib])
    teff[ib].SetMarkerStyle(20)
    teff[ib].SetMarkerColor(colors_emu_Dxy[ib])
    opt = "PZsame" if ib>0 else "PZsame"
    teff[ib].Draw(opt)
    legc1.AddEntry(teff[ib],string_Dxy[ib],"p")
legc1.AddEntry(teffu[0],"legacy BMTF","p")
legc1.Draw()

c2 = TCanvas("c2")
addoverflow(hh1)
addoverflow(hh2)
addoverflow(hh3)
addoverflow(hh4)
normalise(hh1)
normalise(hh2)
normalise(hh3)
normalise(hh4)
hh1.SetLineColor(1)
hh2.SetLineColor(3)
hh3.SetLineColor(2)
hh4.SetLineColor(4)
hh1.SetLineWidth(2)
hh2.SetLineWidth(2)
hh3.SetLineWidth(2)
hh4.SetLineWidth(2)
#hh1.SetMaximum(1.1*max(hh1.GetMaximum(),hh2.GetMaximum(),hh3.GetMaximum(),hh4.GetMaximum()))
hh1.SetMaximum(0.3)
leg = TLegend(0.5,0.7,0.9,0.9)
leg.AddEntry(hh1,"generator","l")
leg.AddEntry(hh2,"legacy BMTF","l")
leg.AddEntry(hh3,"kBMTF constrained","l")
leg.AddEntry(hh4,"kBMTF max(const.,unconst.)","l")

hh1.SetTitle("H to 2X to 4mu (MH/MX/ctau = 200 GeV / 50 GeV / 2000 mm);gen or L1 pT;normalised to area")
hh1.Draw("hist")
hh2.Draw("hist same")
hh3.Draw("hist same")
hh4.Draw("hist same")
leg.Draw()

c3 = TCanvas("c3")
c3.Divide(2,2)
for ib, b in enumerate(bins_Dxy):
    c3.cd(ib+1)
    addoverflow(h_pt_gen[ib])
    addoverflow(h_pt_unp[ib])
    addoverflow(h_pt_emu[ib])
    addoverflow(h_pt_emud[ib])
    normalise(h_pt_gen[ib])
    normalise(h_pt_unp[ib])
    normalise(h_pt_emu[ib])
    normalise(h_pt_emud[ib])
    h_pt_gen[ib].SetMaximum(1.1*max(h_pt_gen[ib].GetMaximum(),h_pt_unp[ib].GetMaximum(),h_pt_emu[ib].GetMaximum(),h_pt_emud[ib].GetMaximum()))
#    h_pt_gen[ib].SetMaximum(0.3)
    h_pt_gen[ib].SetLineColor(1)
    h_pt_gen[ib].SetTitle(string_Dxy[ib]+";gen or L1 pT;normalised to area")
    h_pt_gen[ib].Draw("hist")
    h_pt_unp[ib].SetLineColor(3)
    h_pt_unp[ib].Draw("hist same")
    h_pt_emu[ib].SetLineColor(2)
    h_pt_emu[ib].Draw("hist same")
    h_pt_emud[ib].SetLineColor(4)
    h_pt_emud[ib].Draw("hist same")
    leg.Draw()

c4 = TCanvas("c4")
c4.Divide(2,2)
for ib, b in enumerate(bins_Dxy):
    c4.cd(ib+1)
    h_dpt_unp[ib].SetMaximum(1.1*max(h_dpt_unp[ib].GetMaximum(),h_dpt_emud[ib].GetMaximum()))
    h_dpt_unp[ib].SetLineColor(2)
    h_dpt_unp[ib].SetTitle(string_Dxy[ib]+";pt(L1)-pt(Gen);")
    h_dpt_unp[ib].Draw()
    h_dpt_emud[ib].SetLineColor(4)
    h_dpt_emud[ib].Draw("same")
    h_dptx_unp[ib].SetLineColor(2)
    h_dptx_unp[ib].SetFillColor(2)
    h_dptx_unp[ib].SetFillStyle(3001)
    h_dptx_unp[ib].Draw("same")
    h_dptx_emud[ib].SetLineColor(4)
    h_dptx_emud[ib].SetFillColor(4)
    h_dptx_emud[ib].SetFillStyle(3002)
    h_dptx_emud[ib].Draw("same")


c5 = TCanvas("c5")
c5.SetGrid(1,1)
dummy5 = TH1F("dummy5","",10,0,150)
dummy5.SetTitle(";Dxy;Efficiency")
dummy5.SetMaximum(1.)
dummy5.Draw()
legc5 = TLegend(0.1,0.9,0.9,1.)
legc5.SetNColumns(4)
t2effu = TEfficiency(h2_unp_eff,h2_tot)
t2effd = TEfficiency(h2_emu_eff,h2_tot)
t2effu.SetMarkerStyle(24)
t2effd.SetMarkerStyle(20)
t2effu.Draw("PZ same")
t2effd.Draw("PZ same")

t2effux = TEfficiency(h2_unpx_eff,h2_tot)
t2effdx = TEfficiency(h2_emux_eff,h2_tot)
t2effux.SetMarkerStyle(24)
t2effdx.SetMarkerStyle(20)
t2effux.SetMarkerColor(2)
t2effdx.SetMarkerColor(2)
t2effux.SetLineColor(2)
t2effdx.SetLineColor(2)
t2effux.Draw("PZ same")
t2effdx.Draw("PZ same")

t2effux2 = TEfficiency(h2_unpx2_eff,h2_tot)
t2effdx2 = TEfficiency(h2_emux2_eff,h2_tot)
t2effux2.SetMarkerStyle(24)
t2effdx2.SetMarkerStyle(20)
t2effux2.SetMarkerColor(4)
t2effdx2.SetMarkerColor(4)
t2effux2.SetLineColor(4)
t2effdx2.SetLineColor(4)
t2effux2.Draw("PZ same")
t2effdx2.Draw("PZ same")

legc5.AddEntry(t2effdx,"no cut on pT(L1)","p")
legc5.AddEntry(t2effd,"pT(L1) > 10. GeV","p")
legc5.AddEntry(t2effdx2,"pT(L1) > 15. GeV","p")
#legc5.AddEntry(t2effd,"kBMTF max","p")
legc5.AddEntry(t2effu,"legacy BMTF","p")
legc5.Draw()

c1.SaveAs("c1.png")
c1.SaveAs("c1.root")
c2.SaveAs("c2.png")
c2.SaveAs("c2.root")
c3.SaveAs("c3.png")
c3.SaveAs("c3.root")
c4.SaveAs("c4.png")
c4.SaveAs("c4.root")
c5.SaveAs("c5.png")
c5.SaveAs("c5.root")

'''

c6 = TCanvas("c6")
c6.Divide(2,2)
leg6 = TLegend(0.6,0.5,0.9,0.9)
leg6.AddEntry(h_dxy[0],"L1 muon matched to gen","f")
leg6.AddEntry(h_dxy_gen20[0],"+ pt(gen) > 20 GeV","f")
leg6.AddEntry(h_dxy_eff10[0],"+ pt(L1) > 10 GeV","f")
for ib, b in enumerate(bins_Dxy):
    c6.cd(ib+1)
    h_dxy[ib].SetLineColor(1)
    h_dxy[ib].SetTitle(string_Dxy[ib]+";dxy(L1)-dxy(Gen);")
    h_dxy[ib].Draw()
    h_dxy_gen20[ib].SetLineColor(4)
    h_dxy_gen20[ib].SetFillColor(4)
    h_dxy_gen20[ib].SetFillStyle(3001)
    h_dxy_gen20[ib].Draw("same")
    h_dxy_eff10[ib].SetLineColor(2)
    h_dxy_eff10[ib].SetFillColor(2)
    h_dxy_eff10[ib].SetFillStyle(3002)
    h_dxy_eff10[ib].Draw("same")
    leg6.Draw()
c6.SaveAs("c6.png")
c6.SaveAs("c6.root")
    
c7 = TCanvas("c7")
c7.SetGrid(1,1)
dummy7 = TH1F("dummy7","",10,0,150)
dummy7.SetTitle(";Dxy;Efficiency")
dummy7.SetMaximum(1.)
dummy7.Draw()
legc7 = TLegend(0.1,0.9,0.9,1.)
legc7.SetNColumns(5)
t2effu = TEfficiency(h2_unp_eff,h2_tot)
t2effd = TEfficiency(h2_emu_eff,h2_tot)
t2effu.SetMarkerStyle(24)
t2effd.SetMarkerStyle(20)
t2effu.Draw("PZ same")
t2effd.Draw("PZ same")

t2effd15 = TEfficiency(h2_emu_l1dxy15_eff,h2_tot)
t2effd15.SetMarkerStyle(20)
t2effd15.SetMarkerColor(2)
t2effd15.SetLineColor(2)
t2effd15.Draw("PZ same")

t2effd45 = TEfficiency(h2_emu_l1dxy45_eff,h2_tot)
t2effd45.SetMarkerStyle(20)
t2effd45.SetMarkerColor(4)
t2effd45.SetLineColor(4)
t2effd45.Draw("PZ same")

t2effd75 = TEfficiency(h2_emu_l1dxy75_eff,h2_tot)
t2effd75.SetMarkerStyle(20)
t2effd75.SetMarkerColor(6)
t2effd75.SetLineColor(6)
t2effd75.Draw("PZ same")


legc7.AddEntry(t2effd,"no Dxy(L1) cut","p")
legc7.AddEntry(t2effd15,"Dxy(L1) > 15. GeV","p")
legc7.AddEntry(t2effd45,"Dxy(L1) > 45. GeV","p")
legc7.AddEntry(t2effd75,"Dxy(L1) > 75. GeV","p")
legc7.AddEntry(t2effu,"legacy BMTF","p")
legc7.Draw()

legc7bis = TLegend(0.7,0.7,0.85,0.85)
legc7bis.AddEntry(0,"pT(gen) > 20 GeV","")
legc7bis.AddEntry(0,"pT(L1) > 10 GeV","")
legc7bis.Draw()

c7.SaveAs("c7.png")
c7.SaveAs("c7.root")

