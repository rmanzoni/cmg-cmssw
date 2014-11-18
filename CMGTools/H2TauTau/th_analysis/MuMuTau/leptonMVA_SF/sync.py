import math, sys, array, optparse
import os, ROOT, random
import numpy as num
from ROOT import TFile, TH1F, gDirectory, TMVA, TTree, Double, TLorentzVector, Double
import CMGTools.H2TauTau.config_sf as tool
from CMGTools.RootTools.utils.DeltaR import deltaR,deltaPhi

parser = optparse.OptionParser()
parser.add_option('--phys', action="store", dest="phys", default='data')
options, args = parser.parse_args()

if "/smearer_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/H2TauTau/python/proto/plotter/smearer.cc+" % os.environ['CMSSW_BASE']);
if "/mcCorrections_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/H2TauTau/python/proto/plotter/mcCorrections.cc+" % os.environ['CMSSW_BASE']);

print '[INFO] Physics Proecss = ', options.phys

mva_muon_th = 0.5

mva_muonreader = TMVA.Reader("!Color:Silent=T:Verbose=F")
mva_mvar_map   = {}

for var in ['bdt_muon_dxy','bdt_muon_dz','bdt_muon_dB3D', 'bdt_muon_mva_ch_iso','bdt_muon_mva_neu_iso','bdt_muon_mva_jet_dr','bdt_muon_mva_ptratio','bdt_muon_mva_csv']:
    mva_mvar_map[var] = array.array('f',[0])
    mva_muonreader.AddVariable(var, mva_mvar_map[var])

mva_muonreader.BookMVA('mva_muon_data', 'training/weights/TMVAClassification_BDT_muon.weights.xml')


process = [options.phys]

db = tool.ReadFile(process, 'mmt')
filedict = db.returnFile()

    
if __name__ == '__main__':

    outputfile = 'root_process/Zmumu_' + options.phys + '.root'
    file = TFile(outputfile,'recreate')
    t = TTree('Tree','Tree')
        
    muon_pt = num.zeros(1, dtype=float)
    muon_eta = num.zeros(1, dtype=float)
    muon_phi = num.zeros(1, dtype=float)
    muon_mva = num.zeros(1, dtype=float)
    muon_flag = num.zeros(1, dtype=int)
    muon_dxy = num.zeros(1, dtype=float)
    muon_dz = num.zeros(1, dtype=float)
    muon_dB3D = num.zeros(1, dtype=float)
    muon_ch_iso = num.zeros(1, dtype=float)
    muon_neu_iso = num.zeros(1, dtype=float)
    muon_jet_dr = num.zeros(1, dtype=float)
    muon_ptratio = num.zeros(1, dtype=float)
    muon_csv = num.zeros(1, dtype=float)

    smuon_pt = num.zeros(1, dtype=float)
    smuon_eta = num.zeros(1, dtype=float)
    smuon_phi = num.zeros(1, dtype=float)
    smuon_mva = num.zeros(1, dtype=float)
    smuon_flag = num.zeros(1, dtype=int)    
    smuon_dxy = num.zeros(1, dtype=float)
    smuon_dz = num.zeros(1, dtype=float)
    smuon_dB3D = num.zeros(1, dtype=float)
    smuon_ch_iso = num.zeros(1, dtype=float)
    smuon_neu_iso = num.zeros(1, dtype=float)
    smuon_jet_dr = num.zeros(1, dtype=float)
    smuon_ptratio = num.zeros(1, dtype=float)
    smuon_csv = num.zeros(1, dtype=float)

    evt_weight = num.zeros(1, dtype=float)
    evt_id = num.zeros(1, dtype=int)
    evt_top_weight = num.zeros(1, dtype=float)
    evt_Mmm = num.zeros(1, dtype=float)
    evt_met = num.zeros(1, dtype=float)
    evt_isMC = num.zeros(1, dtype=int)
    evt_run = num.zeros(1, dtype=int)
    evt_evt = num.zeros(1, dtype=int)
    evt_lum = num.zeros(1, dtype=int)
    
    t.Branch('muon_pt',muon_pt,'muon_pt/D')
    t.Branch('muon_eta',muon_eta,'muon_eta/D')
    t.Branch('muon_phi',muon_phi,'muon_phi/D')
    t.Branch('muon_mva', muon_mva, 'muon_mva/D')
    t.Branch('muon_flag', muon_flag, 'muon_flag/I')
    t.Branch('muon_dxy', muon_dxy, 'muon_dxy/D')
    t.Branch('muon_dz', muon_dz, 'muon_dz/D')
    t.Branch('muon_dB3D', muon_dB3D, 'muon_dB3D/D')
    t.Branch('muon_ch_iso', muon_ch_iso, 'muon_ch_iso/D')
    t.Branch('muon_neu_iso', muon_neu_iso, 'muon_neu_iso/D')
    t.Branch('muon_jet_dr', muon_jet_dr, 'muon_jet_dr/D')
    t.Branch('muon_ptratio', muon_ptratio, 'muon_ptratio/D')
    t.Branch('muon_csv', muon_csv, 'muon_csv/D')
                            

    t.Branch('smuon_pt',smuon_pt,'smuon_pt/D')
    t.Branch('smuon_eta',smuon_eta,'smuon_eta/D')
    t.Branch('smuon_phi',smuon_phi,'smuon_phi/D')
    t.Branch('smuon_mva', smuon_mva, 'smuon_mva/D')
    t.Branch('smuon_flag', smuon_flag, 'smuon_flag/I')
    t.Branch('smuon_dxy', smuon_dxy, 'smuon_dxy/D')
    t.Branch('smuon_dz', smuon_dz, 'smuon_dz/D')
    t.Branch('smuon_dB3D', smuon_dB3D, 'smuon_dB3D/D')
    t.Branch('smuon_ch_iso', smuon_ch_iso, 'smuon_ch_iso/D')
    t.Branch('smuon_neu_iso', smuon_neu_iso, 'smuon_neu_iso/D')
    t.Branch('smuon_jet_dr', smuon_jet_dr, 'smuon_jet_dr/D')
    t.Branch('smuon_ptratio', smuon_ptratio, 'smuon_ptratio/D')
    t.Branch('smuon_csv', smuon_csv, 'smuon_csv/D')


    t.Branch('evt_weight', evt_weight, 'evt_weight/D')
    t.Branch('evt_id', evt_id, 'evt_id/I')
    t.Branch('evt_top_weight', evt_top_weight, 'evt_top_weight/D')
    t.Branch('evt_Mmm', evt_Mmm, 'evt_Mmm/D')
    t.Branch('evt_met', evt_met, 'evt_met/D')
    t.Branch('evt_isMC', evt_isMC, 'evt_isMC/I')
    t.Branch('evt_run', evt_run, 'evt_run/I')
    t.Branch('evt_evt', evt_evt, 'evt_evt/I')
    t.Branch('evt_lum', evt_lum, 'evt_lum/I')


   
    for index, ifile in enumerate(filedict):

        pname = ifile[0]
        filename = ifile[1]
        lum_weight = ifile[2]
        ptype = ifile[3]

        print '[INFO] ', index, filename, 'is processing'

        myfile = TFile(filename)

        main = gDirectory.Get('H2TauTauTreeProducerMMT')
        mchain = gDirectory.Get('H2TauTauTreeProducerMMT_muon')
        gchain = gDirectory.Get('H2TauTauTreeProducerMMT_gen')
        
        ptr_m = 0        
        ptr_ng = 0

        Total = main.GetEntries()

        top_inclusive = 1.
        
        if pname == 'tt0l' or pname=='tt1l' or pname=='tt2l':
            total_entry = 0
            
            for jentry in xrange(main.GetEntries()):

                ientry = main.LoadTree(jentry)
                nb = main.GetEntry(jentry)

                total_entry += tool.returnTopWeight(pname, main.top_pt, main.atop_pt)

            top_inclusive = main.GetEntries()/total_entry
            
        for jentry in xrange(main.GetEntries()):

            ientry = main.LoadTree(jentry)
            nb = main.GetEntry(jentry)

            evt_flag = False

            if jentry%20000==0:
                print '[INFO]', jentry, '/', main.GetEntries() #nmuon, nelectron, ntau, nvmuon, nvelectron, nvtau


            nmuon      = int(main.nmuon)

            if pname != 'data':
                ngen = int(main.nGen)


            gp = []
            if pname != 'data':
                for igen in xrange(ptr_ng, ptr_ng+ngen):
                    
                    gchain.LoadTree(igen)
                    gchain.GetEntry(igen)
                    
                    gj = tool.easyobj_gen(gchain.gen_pt,
                                          gchain.gen_eta,
                                          gchain.gen_phi,
                                          gchain.gen_pdgid)
                    gp.append(gj)



            
            signal_muon = []
            
            for im in xrange(ptr_m, ptr_m + nmuon):
                mchain.LoadTree(im)
                mchain.GetEntry(im)
                
                muon_ipdg = -99
                muon_min_dr = 100

                for gen in gp:
                    _dr_ = deltaR(gen.eta, gen.phi, mchain.muon_eta, mchain.muon_phi)
                    if _dr_ < 0.5 and muon_min_dr > _dr_:
                        muon_min_dr = _dr_
                        muon_ipdg = gen.pdgid

                matchid = 0
                matchany = 0
                if abs(muon_ipdg)==5:
                    matchany = 2


                mva_mvar_map['bdt_muon_mva_ch_iso'][0] = mchain.muon_mva_ch_iso
                mva_mvar_map['bdt_muon_mva_neu_iso'][0] = mchain.muon_mva_neu_iso
                mva_mvar_map['bdt_muon_mva_csv'][0] = mchain.muon_mva_csv

                cor_dxy = mchain.muon_mva_dxy
                cor_dz = mchain.muon_mva_dz
                cor_jet_dr = mchain.muon_mva_jet_dr
                cor_ptratio = mchain.muon_mva_ptratio
                cor_dB3D = abs(mchain.muon_dB3D)

                if pname != 'data':
                    cor_dxy = ROOT.scaleDxyMC(mchain.muon_mva_dxy, int(muon_ipdg), mchain.muon_pt, mchain.muon_eta, matchid, matchany)
                    cor_dz = ROOT.scaleDzMC(mchain.muon_mva_dz, int(muon_ipdg), mchain.muon_pt, mchain.muon_eta, matchid, matchany)
                    cor_jet_dr = ROOT.correctJetDRMC(mchain.muon_mva_jet_dr, int(muon_ipdg), mchain.muon_pt, mchain.muon_eta, matchid, matchany)
                    cor_ptratio = ROOT.correctJetPtRatioMC(mchain.muon_mva_ptratio, int(muon_ipdg), mchain.muon_pt, mchain.muon_eta, matchid, matchany)
                    cor_dB3D = ROOT.scaleSip3dMC(abs(mchain.muon_dB3D), int(muon_ipdg), mchain.muon_pt, mchain.muon_eta, matchid, matchany)
                    
                mva_mvar_map['bdt_muon_dxy'][0] = cor_dxy
                mva_mvar_map['bdt_muon_dz'][0] = cor_dz
                mva_mvar_map['bdt_muon_dB3D'][0] = math.log(cor_dB3D)
                mva_mvar_map['bdt_muon_mva_jet_dr'][0] = cor_jet_dr
                mva_mvar_map['bdt_muon_mva_ptratio'][0] = cor_ptratio
                
                mva_iso_muon = mva_muonreader.EvaluateMVA('mva_muon_data')

                if not mchain.muon_id:
                    continue


                _flag_ = (mva_iso_muon > mva_muon_th)

                muon = tool.mobj(mchain.muon_pt,
                                 mchain.muon_eta,
                                 mchain.muon_phi,
                                 mchain.muon_mass,                                 
                                 mchain.muon_charge,
                                 mchain.muon_trigmatch,
                                 mchain.muon_trig_weight,
                                 mchain.muon_id_weight,
                                 mva_iso_muon,
                                 cor_dxy,
                                 cor_dz,
                                 cor_dB3D,
                                 mchain.muon_mva_ch_iso,
                                 mchain.muon_mva_neu_iso,
                                 cor_jet_dr,
                                 cor_ptratio,
                                 mchain.muon_mva_csv,
                                 _flag_
                                 )
                

                signal_muon.append(muon)

            
            if not (len(signal_muon)==2):
                    
                ptr_m += nmuon
                if pname != 'data': ptr_ng += ngen
                continue


#            import pdb; pdb.set_trace()

            ptr_m += nmuon
            if pname != 'data': ptr_ng += ngen

            tagid = (random.random() > 0.5)

            tag_muon = signal_muon[int(tagid)]
            prob_muon = signal_muon[int(not tagid)]

            if not tag_muon.flag:
                continue
            
            if not tag_muon.charge*prob_muon.charge==-1:
                continue

            if not (tag_muon.pt > 40. and tag_muon.trigmatch and prob_muon.trigmatch):
                continue

            if tool.diobj(tag_muon, prob_muon).returnmass() < 20:
                continue

            weight = 1.
            isMC = False
                        
            if pname == 'data':
                pass
            else:
                weight = main.weight*tag_muon.trig*tag_muon.id*prob_muon.trig*prob_muon.id*lum_weight
                isMC = True                    

            muon_pt [0] = tag_muon.pt
            muon_eta [0] = tag_muon.eta
            muon_phi [0] = tag_muon.phi

            muon_dxy[0] = tag_muon.dxy
            muon_dz[0] = tag_muon.dz
            muon_dB3D[0] = math.log(abs(tag_muon.dB3D))
            muon_ch_iso[0] = tag_muon.chiso
            muon_neu_iso[0] = tag_muon.niso
            muon_jet_dr[0] = tag_muon.jet_dr
            muon_ptratio[0] = tag_muon.ptratio
            muon_csv[0] = tag_muon.csv



            smuon_pt [0] = prob_muon.pt
            smuon_eta [0] = prob_muon.eta
            smuon_phi [0] = prob_muon.phi

            smuon_dxy[0] = prob_muon.dxy
            smuon_dz[0] = prob_muon.dz
            smuon_dB3D[0] = math.log(abs(prob_muon.dB3D))
            smuon_ch_iso[0] = prob_muon.chiso
            smuon_neu_iso[0] = prob_muon.niso
            smuon_jet_dr[0] = prob_muon.jet_dr
            smuon_ptratio[0] = prob_muon.ptratio
            smuon_csv[0] = prob_muon.csv


            muon_mva[0] = tag_muon.mva
            muon_flag[0] = tag_muon.flag

            smuon_mva[0] = prob_muon.mva
            smuon_flag[0] = prob_muon.flag

            evt_weight [0] = weight
            evt_Mmm [0] = tool.diobj(tag_muon, prob_muon).returnmass()            
            evt_isMC [0] = isMC
            evt_id [0] = ptype
            evt_run[0] = main.run
            evt_evt[0] = main.evt
            evt_lum[0] = main.lumi
            evt_met[0] = main.pfmet
            
            t.Fill()

        
    file.Write()
    file.Close()




    
