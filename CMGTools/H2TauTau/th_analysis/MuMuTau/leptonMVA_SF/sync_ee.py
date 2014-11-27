import math, sys, array, optparse
import os, ROOT, random
import numpy as num
from ROOT import TFile, TH1F, gDirectory, TMVA, TTree, Double, TLorentzVector, Double
import CMGTools.H2TauTau.config as tool
from CMGTools.RootTools.utils.DeltaR import deltaR,deltaPhi

parser = optparse.OptionParser()
parser.add_option('--phys', action="store", dest="phys", default='data')
options, args = parser.parse_args()

if "/smearer_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/H2TauTau/python/proto/plotter/smearer.cc+" % os.environ['CMSSW_BASE']);
if "/mcCorrections_cc.so" not in ROOT.gSystem.GetLibraries(): 
    ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/H2TauTau/python/proto/plotter/mcCorrections.cc+" % os.environ['CMSSW_BASE']);

print '[INFO] Physics Proecss = ', options.phys

mva_electron_barrel = 0.2
mva_electron_endcap = 0.2

mva_electronreader = TMVA.Reader("!Color:Silent=T:Verbose=F")
mva_mvar_map   = {}

for var in ['bdt_electron_mva_score', 'bdt_electron_mva_ch_iso','bdt_electron_mva_neu_iso','bdt_electron_mva_jet_dr','bdt_electron_mva_ptratio','bdt_electron_mva_csv', 'bdt_electron_dB3D']:
#for var in ['bdt_electron_mva_score', 'bdt_electron_mva_ch_iso','bdt_electron_mva_neu_iso','bdt_electron_mva_jet_dr','bdt_electron_mva_ptratio','bdt_electron_mva_csv', 'bdt_electron_sip3D']:
    mva_evar_map[var] = array.array('f',[0])
    mva_electronreader.AddVariable(var, mva_evar_map[var])

mva_electronreader.BookMVA('mva_electron_data', 'training/weights/TMVAClassification_BDT_electron.weights.xml')

process = [options.phys]

db = tool.ReadFile(process, 'mmt')
filedict = db.returnFile()

    
if __name__ == '__main__':

    outputfile = 'root_process/Zmumu_' + options.phys + '.root'
    file = TFile(outputfile,'recreate')
    t = TTree('Tree','Tree')
        
    electron_pt = num.zeros(1, dtype=float)
    electron_eta = num.zeros(1, dtype=float)
    electron_phi = num.zeros(1, dtype=float)
    electron_mass = num.zeros(1, dtype=float)
    electron_jetpt = num.zeros(1, dtype=float)
    electron_dxy = num.zeros(1, dtype=float)
    electron_dz = num.zeros(1, dtype=float)
    electron_dB3D = num.zeros(1, dtype=float)
    electron_id = num.zeros(1, dtype=int)
    electron_iso = num.zeros(1, dtype=int)
    electron_reliso = num.zeros(1, dtype=float)
    electron_MT = num.zeros(1, dtype=float)
    electron_charge = num.zeros(1, dtype=int)
    electron_pdg = num.zeros(1, dtype=int)
    electron_ptratio = num.zeros(1, dtype=float)
    electron_mva = num.zeros(1, dtype=float)
    electron_mva_ch_iso = num.zeros(1, dtype=float)
    electron_mva_neu_iso = num.zeros(1, dtype=float)
    electron_mva_jet_dr = num.zeros(1, dtype=float)
    electron_mva_ptratio = num.zeros(1, dtype=float)
    electron_mva_csv = num.zeros(1, dtype=float)
    electron_new_mva = num.zeros(1, dtype=float)
    electron_flag = num.zeros(1, dtype=int)
    
    selectron_pt = num.zeros(1, dtype=float)
    selectron_eta = num.zeros(1, dtype=float)
    selectron_phi = num.zeros(1, dtype=float)
    selectron_mass = num.zeros(1, dtype=float)
    selectron_jetpt = num.zeros(1, dtype=float)
    selectron_dxy = num.zeros(1, dtype=float)
    selectron_dz = num.zeros(1, dtype=float)
    selectron_dB3D = num.zeros(1, dtype=float)
    selectron_id = num.zeros(1, dtype=int)
    selectron_iso = num.zeros(1, dtype=int)
    selectron_reliso = num.zeros(1, dtype=float)
    selectron_MT = num.zeros(1, dtype=float)
    selectron_charge = num.zeros(1, dtype=int)
    selectron_pdg = num.zeros(1, dtype=int)
    selectron_ptratio = num.zeros(1, dtype=float)
    selectron_mva = num.zeros(1, dtype=float)
    selectron_mva_ch_iso = num.zeros(1, dtype=float)
    selectron_mva_neu_iso = num.zeros(1, dtype=float)
    selectron_mva_jet_dr = num.zeros(1, dtype=float)
    selectron_mva_ptratio = num.zeros(1, dtype=float)
    selectron_mva_csv = num.zeros(1, dtype=float)
    selectron_new_mva = num.zeros(1, dtype=float)
    selectron_flag = num.zeros(1, dtype=int)    

    evt_weight = num.zeros(1, dtype=float)
    evt_top_weight = num.zeros(1, dtype=float)
    evt_Mmm = num.zeros(1, dtype=float)
    evt_met = num.zeros(1, dtype=float)
    evt_isMC = num.zeros(1, dtype=int)
    evt_id = num.zeros(1, dtype=int)
    evt_run = num.zeros(1, dtype=int)
    evt_evt = num.zeros(1, dtype=int)
    evt_lum = num.zeros(1, dtype=int)
    
    t.Branch('electron_pt',electron_pt,'electron_pt/D')
    t.Branch('electron_eta',electron_eta,'electron_eta/D')
    t.Branch('electron_phi',electron_phi,'electron_phi/D')
    t.Branch('electron_mass', electron_mass, 'electron_mass/D')
    t.Branch('electron_jetpt',electron_jetpt, 'electron_jetpt/D')
    t.Branch('electron_dxy',electron_dxy, 'electron_dxy/D')
    t.Branch('electron_dz',electron_dz, 'electron_dz/D')
    t.Branch('electron_dB3D',electron_dB3D, 'electron_dB3D/D')
    t.Branch('electron_id', electron_id, 'electron_id/I')
    t.Branch('electron_iso', electron_iso, 'electron_iso/I')
    t.Branch('electron_reliso', electron_reliso, 'electron_reliso/D')
    t.Branch('electron_MT', electron_MT, 'electron_MT/D')
    t.Branch('electron_charge', electron_charge, 'electron_charge/I')
    t.Branch('electron_pdg',electron_pdg,'electron_pdg/I')
    t.Branch('electron_ptratio', electron_ptratio, 'electron_ptratio/D')
    t.Branch('electron_mva', electron_mva, 'electron_mva/D')
    t.Branch('electron_mva_ch_iso', electron_mva_ch_iso, 'electron_mva_ch_iso/D')
    t.Branch('electron_mva_neu_iso', electron_mva_neu_iso, 'electron_mva_neu_iso/D')
    t.Branch('electron_mva_jet_dr', electron_mva_jet_dr, 'electron_mva_jet_dr/D')
    t.Branch('electron_mva_ptratio', electron_mva_ptratio, 'electron_mva_ptratio/D')
    t.Branch('electron_mva_csv', electron_mva_csv, 'electron_mva_csv/D')
    t.Branch('electron_new_mva', electron_new_mva, 'electron_new_mva/D')
    t.Branch('electron_flag', electron_flag, 'electron_flag/I')

    t.Branch('selectron_pt',selectron_pt,'selectron_pt/D')
    t.Branch('selectron_eta',selectron_eta,'selectron_eta/D')
    t.Branch('selectron_phi',selectron_phi,'selectron_phi/D')
    t.Branch('selectron_mass', selectron_mass, 'selectron_mass/D')
    t.Branch('selectron_jetpt',selectron_jetpt, 'selectron_jetpt/D')
    t.Branch('selectron_dxy',selectron_dxy, 'selectron_dxy/D')
    t.Branch('selectron_dz',selectron_dz, 'selectron_dz/D')
    t.Branch('selectron_dB3D',selectron_dB3D, 'selectron_dB3D/D')
    t.Branch('selectron_id', selectron_id, 'selectron_id/I')
    t.Branch('selectron_iso', selectron_iso, 'selectron_iso/I')
    t.Branch('selectron_reliso', selectron_reliso, 'selectron_reliso/D')
    t.Branch('selectron_MT', selectron_MT, 'selectron_MT/D')
    t.Branch('selectron_charge', selectron_charge, 'selectron_charge/I')
    t.Branch('selectron_pdg',selectron_pdg,'selectron_pdg/I')
    t.Branch('selectron_ptratio', selectron_ptratio, 'selectron_ptratio/D')
    t.Branch('selectron_mva', selectron_mva, 'selectron_mva/D')
    t.Branch('selectron_mva_ch_iso', selectron_mva_ch_iso, 'selectron_mva_ch_iso/D')
    t.Branch('selectron_mva_neu_iso', selectron_mva_neu_iso, 'selectron_mva_neu_iso/D')
    t.Branch('selectron_mva_jet_dr', selectron_mva_jet_dr, 'selectron_mva_jet_dr/D')
    t.Branch('selectron_mva_ptratio', selectron_mva_ptratio, 'selectron_mva_ptratio/D')
    t.Branch('selectron_mva_csv', selectron_mva_csv, 'selectron_mva_csv/D')
    t.Branch('selectron_new_mva', selectron_new_mva, 'selectron_new_mva/D')
    t.Branch('selectron_flag', selectron_flag, 'selectron_flag/I')

    t.Branch('evt_weight', evt_weight, 'evt_weight/D')
    t.Branch('evt_top_weight', evt_top_weight, 'evt_top_weight/D')
    t.Branch('evt_Mmm', evt_Mmm, 'evt_Mmm/D')
    t.Branch('evt_met', evt_met, 'evt_met/D')
    t.Branch('evt_isMC', evt_isMC, 'evt_isMC/I')
    t.Branch('evt_id', evt_id, 'evt_id/I')
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
        echain = gDirectory.Get('H2TauTauTreeProducerMMT_electron')
        jchain = gDirectory.Get('H2TauTauTreeProducerMMT_jet')
        gchain = gDirectory.Get('H2TauTauTreeProducerMMT_gen')
        
        ptr_m = 0        
        ptr_nj = 0
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
                print '[INFO]', jentry, '/', main.GetEntries() #nelectron, nelectron, ntau, nvelectron, nvelectron, nvtau



            nelectron      = int(main.nelectron)
            nelectron  = int(main.nelectron)
            njets      = int(main.nJets)

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



            
            signal_electron = []
            
            for im in xrange(ptr_m, ptr_m + nelectron):
                echain.LoadTree(im)
                echain.GetEntry(im)
                
                electron_ipdg = -99
                electron_min_dr = 100

                for gen in gp:
                    _dr_ = deltaR(gen.eta, gen.phi, echain.electron_eta, echain.electron_phi)
                    if _dr_ < 0.5 and electron_min_dr > _dr_:
                        electron_min_dr = _dr_
                        electron_ipdg = gen.pdgid

                matchid = 0
                matchany = 0
                if abs(electron_ipdg)==5:
                    matchany = 2


                mva_mvar_map['bdt_electron_mva_ch_iso'][0] = echain.electron_mva_ch_iso
                mva_mvar_map['bdt_electron_mva_neu_iso'][0] = echain.electron_mva_neu_iso
                mva_mvar_map['bdt_electron_mva_csv'][0] = echain.electron_mva_csv

                cor_dxy = echain.electron_mva_dxy
                cor_dz = echain.electron_mva_dz
                cor_jet_dr = echain.electron_mva_jet_dr
                cor_ptratio = echain.electron_mva_ptratio
                cor_sip3D = echain.electron_sip3D
                cor_dB3D = abs(echain.electron_dB3D)

                if pname != 'data':
                    cor_dxy = ROOT.scaleDxyMC(echain.electron_mva_dxy, int(electron_ipdg), echain.electron_pt, echain.electron_eta, matchid, matchany)
                    cor_dz = ROOT.scaleDzMC(echain.electron_mva_dz, int(electron_ipdg), echain.electron_pt, echain.electron_eta, matchid, matchany)
                    cor_jet_dr = ROOT.correctJetDRMC(echain.electron_mva_jet_dr, int(electron_ipdg), echain.electron_pt, echain.electron_eta, matchid, matchany)
                    cor_ptratio = ROOT.correctJetPtRatioMC(echain.electron_mva_ptratio, int(electron_ipdg), echain.electron_pt, echain.electron_eta, matchid, matchany)
                    cor_sip3D = ROOT.scaleSip3dMC(echain.electron_sip3D, int(electron_ipdg), echain.electron_pt, echain.electron_eta, matchid, matchany)
                    cor_dB3D = ROOT.scaleSip3dMC(abs(echain.electron_dB3D), int(electron_ipdg), echain.electron_pt, echain.electron_eta, matchid, matchany)
                    
                mva_mvar_map['bdt_electron_dxy'][0] = cor_dxy
                mva_mvar_map['bdt_electron_dz'][0] = cor_dz
                mva_mvar_map['bdt_electron_dB3D'][0] = math.log(cor_dB3D)
                mva_mvar_map['bdt_electron_mva_jet_dr'][0] = cor_jet_dr
                mva_mvar_map['bdt_electron_mva_ptratio'][0] = cor_ptratio
                
                mva_iso_electron = mva_electronreader.EvaluateMVA('mva_electron_data')

                if not echain.electron_id:
                    continue


                _flag_ = ((abs(echain.electron_eta) < 1.479 and mva_iso_electron > mva_electron_barrel) or (abs(echain.electron_eta) > 1.479 and mva_iso_electron > mva_electron_endcap))

                electron = tool.mobj(echain.electron_pt,
                                 echain.electron_eta,
                                 echain.electron_phi,
                                 echain.electron_mass,
                                 echain.electron_jetpt,
                                 echain.electron_njet,
                                 echain.electron_charge,
                                 echain.electron_trigmatch,
                                 echain.electron_trig_weight,
                                 echain.electron_id_weight,
                                 echain.electron_id,
                                 echain.electron_iso,
                                 echain.electron_reliso,
                                 echain.electron_MT,
                                 cor_dxy,
                                 echain.electron_mva_dxy,
                                 cor_dz,
                                 echain.electron_mva_dz,
                                 cor_dB3D,
                                 abs(echain.electron_dB3D),
                                 cor_sip3D,
                                 echain.electron_sip3D,
                                 echain.electron_jetcsv,
                                 echain.electron_jetcsv_10,
                                 echain.electron_mva,
                                 echain.electron_mva_ch_iso,
                                 echain.electron_mva_neu_iso,
                                 cor_jet_dr,
                                 echain.electron_mva_jet_dr,
                                 cor_ptratio,
                                 echain.electron_mva_ptratio,
                                 echain.electron_mva_csv,
                                 mva_iso_electron,
                                 _flag_
                                 )
                
                
                signal_electron.append(electron)



            
            if not (len(signal_electron)==2):
                    
                ptr_m += nelectron
                if pname != 'data': ptr_ng += ngen
                ptr_nj += njets
                continue


#            import pdb; pdb.set_trace()

            ptr_m += nelectron
            ptr_nj += njets
            if pname != 'data': ptr_ng += ngen

            tagid = (random.random() > 0.5)

            tag_electron = signal_electron[int(tagid)]
            prob_electron = signal_electron[int(not tagid)]

            
            if not tag_electron.charge*prob_electron.charge==-1:
#                print 'OS failes'
                continue

            if not (tag_electron.pt > 20. and tag_electron.trigmatch and prob_electron.trigmatch):
#                print 'trigger matching failes'
                continue

                    

#            print 'Mll', tool.diobj(tag_electron, prob_electron).returnmass()

            weight = 1.
            isMC = False
                        
            if pname == 'data':
                pass
            else:
                weight = main.weight*tag_electron.trig*tag_electron.id*prob_electron.trig*prob_electron.id*lum_weight
                isMC = True                    

            electron_pt [0] = tag_electron.pt
            electron_eta [0] = tag_electron.eta
            electron_phi [0] = tag_electron.phi
            electron_mass [0] = tag_electron.mass
            electron_jetpt [0] = tag_electron.jetpt
            electron_id [0] = tag_electron.isid
            electron_iso [0] = tag_electron.isiso
            electron_reliso [0] = tag_electron.reliso
            electron_MT [0] = tag_electron.MT
            electron_charge [0] = tag_electron.charge
            electron_dxy [0] = tag_electron.dxy
            electron_dz [0] = tag_electron.dz
            electron_dB3D [0] = math.log(tag_electron.dB3D)


            electron_ipdg = 0
            electron_min_dr = 100
            
            if pname != 'data':
                
                for gen in gp:
                    if tag_electron.returndR(gen) < 0.5:
                        electron_min_dr = tag_electron.returndR(gen)
                        electron_ipdg = gen.pdgid
                        
            electron_pdg[0] = electron_ipdg


            selectron_pt [0] = prob_electron.pt
            selectron_eta [0] = prob_electron.eta
            selectron_phi [0] = prob_electron.phi
            selectron_mass [0] = prob_electron.mass
            selectron_jetpt [0] = prob_electron.jetpt
            selectron_id [0] = prob_electron.isid
            selectron_iso [0] = prob_electron.isiso
            selectron_reliso [0] = prob_electron.reliso
            selectron_MT [0] = prob_electron.MT
            selectron_charge [0] = prob_electron.charge
            selectron_dxy [0] = prob_electron.dxy
            selectron_dz [0] = prob_electron.dz
            selectron_dB3D [0] = math.log(prob_electron.dB3D)


            selectron_ipdg = 0
            selectron_min_dr = 100
            
            if pname != 'data':
                
                for gen in gp:
                    if prob_electron.returndR(gen) < 0.5:
                        selectron_min_dr = prob_electron.returndR(gen)
                        selectron_ipdg = gen.pdgid
                        
            selectron_pdg[0] = selectron_ipdg

            if pname == 'tt0l' or pname=='tt1l' or pname=='tt2l':
                evt_weight [0] = weight*top_inclusive*tool.returnTopWeight(pname, main.top_pt, main.atop_pt)
                evt_top_weight [0] = top_inclusive*tool.returnTopWeight(pname, main.top_pt, main.atop_pt)
            else:
                evt_weight [0] = weight
                evt_top_weight [0] = 1.
                    
                

            electron_mva[0] = tag_electron.mva
            electron_mva_ch_iso[0] = tag_electron.mva_ch_iso
            electron_mva_neu_iso[0] = tag_electron.mva_neu_iso
            electron_mva_jet_dr[0] = tag_electron.mva_jet_dr
            electron_mva_ptratio[0] = tag_electron.mva_ptratio
            electron_mva_csv[0] = tag_electron.mva_csv
            electron_new_mva[0] = tag_electron.new_mva
            electron_flag[0] = tag_electron.flag
            
            
            selectron_mva[0] = prob_electron.mva
            selectron_mva_ch_iso[0] = prob_electron.mva_ch_iso
            selectron_mva_neu_iso[0] = prob_electron.mva_neu_iso
            selectron_mva_jet_dr[0] = prob_electron.mva_jet_dr
            selectron_mva_ptratio[0] = prob_electron.mva_ptratio
            selectron_mva_csv[0] = prob_electron.mva_csv
            selectron_new_mva[0] = prob_electron.new_mva
            selectron_flag[0] = prob_electron.flag

            evt_Mmm [0] = tool.diobj(tag_electron, prob_electron).returnmass()            
            evt_isMC [0] = isMC
            evt_id [0] = ptype
            evt_run[0] = main.run
            evt_evt[0] = main.evt
            evt_lum[0] = main.lumi
            evt_met[0] = main.pfmet
            
            t.Fill()

        
    file.Write()
    file.Close()




    
