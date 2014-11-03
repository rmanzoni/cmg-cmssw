import math, sys, array
import numpy as num
from ROOT import TFile, gDirectory, TH1F, gStyle, gROOT, TTree, TMVA, Double
import math, copy, sys, array, optparse

process = ['WW', 'WZ', 'ZZ',
           'tt1l', 'tt2l',
           'DY', 'DY1', 'DY2', 'DY3', 'DY4',
           'Wjet', 'W1jet', 'W2jet', 'W3jet', 'W4jet',
           'tH_YtMinus1', 'TTW', 'TTZ', 'TTH',
           'data']

#process = ['tH_YtMinus1']


#nstep = 20
nstep = 50
scan_width = 0.2


region = ['signal','antiE','antiMu','antiEMu']

directory = 'root_process'


### For options
parser = optparse.OptionParser()
parser.add_option('--kNN', action="store", dest="kNN", default='50')
options, args = parser.parse_args()

#print '[INFO] kNN = ', options.kNN


gROOT.SetBatch(True)

muonreader = [0 for ii in range(len(process))]
electronreader = [0 for ii in range(len(process))]

#mva_muon_barrel = 0.0089
#mva_electron_barrel = 0.0649

#mva_muon_endcap = 0.0621
#mva_electron_endcap = 0.0891

def evaluate(eta, mva, id, th_barrel, th_endcap, lepton, region):

#    import pdb; pdb.set_trace()

    if region in ['signal', 'antiMu', 'antiE', 'antiEMu']:
        pass
    else:
        print 'Non expected region', region

    flag = False
    
    if lepton=='electron':
        if abs(eta) < 1.479:
            if region in ['signal', 'antiMu']:
                if id > 0.5 and mva > th_barrel: flag = True
            else:
                if not (id > 0.5 and mva > th_barrel): flag = True
        elif abs(eta) > 1.479:
            if region in ['signal', 'antiMu']:
                if id > 0.5 and mva > th_endcap: flag = True
            else:
                if not (id > 0.5 and mva > th_endcap): flag = True
                
    elif lepton=='muon':
        if abs(eta) < 1.479:
            if region in ['signal', 'antiE']:
                if id > 0.5 and mva > th_barrel: flag = True
            else:
                if not (id > 0.5 and mva > th_barrel): flag = True
        elif abs(eta) > 1.479:
            if region in ['signal', 'antiE']:
                if id > 0.5 and mva > th_endcap: flag = True
            else:
                if not (id > 0.5 and mva > th_endcap): flag = True

    return flag



def returnkNN(iregion, iprocess, weight_electron, weight_muon):

    kNN_weight = 1.
    if iregion=='antiE':
        if weight_electron==1:
            kNN_weight = 0
            print '[WARNING] 0 weight for e', iprocess, iregion
        else:
            kNN_weight = weight_electron/(1-weight_electron)
    elif iregion=='antiMu':
        if weight_muon==1:
            kNN_weight = 0
            print '[WARNING] warning, 0 weight for mu', iprocess, iregion
        else:
            kNN_weight = weight_muon/(1-weight_muon)
    elif iregion=='antiEMu':
        if weight_electron==1 or weight_muon==1:
            kNN_weight = 0
            print '[WARNING] warning, 0 weight for mu*e', iprocess, iregion
        else:
            kNN_weight = weight_muon*weight_electron/((1-weight_muon)*(1-weight_electron))
    elif iregion=='signal':
        kNN_weight = 1.

    return kNN_weight


for index, pn in enumerate(process):

    e_xml = 'kNN_training/weights/KNN_' + pn + '_electron_' + options.kNN + '.xml'
    m_xml = 'kNN_training/weights/KNN_' + pn + '_muon_' + options.kNN + '.xml'

#    if pn in ['WZ','ZZ','tt1l','tt2l','data']:
    if pn in []:
        pass
    else:
#        print '[INFO] The process', pn, 'uses the kNN weight for data ...'
        e_xml = 'kNN_training/weights/KNN_data_electron_' + options.kNN + '.xml'
        m_xml = 'kNN_training/weights/KNN_data_muon_' + options.kNN + '.xml'

    print '[INFO] xml file = ', e_xml, m_xml

    muonreader[index] = TMVA.Reader("!Color:Silent=T:Verbose=F")
    electronreader[index] = TMVA.Reader("!Color:Silent=T:Verbose=F")        
    mvar_map   = {}
    evar_map   = {}

    
#    for var in ['lepton_pt', 'lepton_kNN_jetpt', 'evt_njet']:
    for var in ['lepton_pt', 'evt_njet']:
        mvar_map[var] = array.array('f',[0])
        muonreader[index].AddVariable(var, mvar_map[var])
        
        evar_map[var] = array.array('f',[0])
        electronreader[index].AddVariable(var, evar_map[var])

    mvaname = 'muon_' + pn
    muonreader[index].BookMVA(mvaname, m_xml)

    mvaname = 'electron_' + pn
    electronreader[index].BookMVA(mvaname, e_xml)




#_muon_barrel_ = 0.0089
#_electron_barrel_ = 0.0649

#_muon_endcap_ = 0.0621
#_electron_endcap_ = 0.0891

_muon_barrel_ = -0.0383
_electron_barrel_ = 0.0557

_muon_endcap_ = -0.0891
_electron_endcap_ = 0.0523

outputfile = 'Scan.root'
    
file = TFile(outputfile,'recreate')
t = TTree('Tree','Tree')

muon_b_mva = num.zeros(1, dtype=num.float32)
muon_e_mva = num.zeros(1, dtype=num.float32)
electron_b_mva = num.zeros(1, dtype=num.float32)
electron_e_mva = num.zeros(1, dtype=num.float32)
nsig = num.zeros(1, dtype=num.float32)
nbkg = num.zeros(1, dtype=num.float32)
nbkg_ttbar = num.zeros(1, dtype=num.float32)
nbkg_diboson = num.zeros(1, dtype=num.float32)
nbkg_ttv = num.zeros(1, dtype=num.float32)
nbkg_tth = num.zeros(1, dtype=num.float32)
nbkg_red = num.zeros(1, dtype=num.float32)
significance = num.zeros(1, dtype=num.float32)
n_ie = num.zeros(1, dtype=num.float32)
n_im = num.zeros(1, dtype=num.float32)
n_max = num.zeros(1, dtype=num.float32)

t.Branch('muon_b_mva',muon_b_mva,'muon_b_mva/F')
t.Branch('muon_e_mva',muon_e_mva,'muon_e_mva/F')
t.Branch('electron_b_mva',electron_b_mva,'electron_b_mva/F')
t.Branch('electron_e_mva',electron_e_mva,'electron_e_mva/F')
t.Branch('nsig',nsig,'nsig/F')
t.Branch('nbkg',nbkg,'nbkg/F')
t.Branch('nbkg_ttbar',nbkg_ttbar,'nbkg_ttbar/F')
t.Branch('nbkg_diboson',nbkg_diboson,'nbkg_diboson/F')
t.Branch('nbkg_ttv',nbkg_ttv,'nbkg_ttv/F')
t.Branch('nbkg_tth',nbkg_tth,'nbkg_tth/F')
t.Branch('nbkg_red',nbkg_red,'nbkg_red/F')
t.Branch('significance',significance,'significance/F')
t.Branch('n_ie',n_ie,'n_ie/F')
t.Branch('n_im',n_im,'n_im/F')
t.Branch('n_max',n_max,'n_max/F')



for ie in range(0, nstep+1):
    for im in range(0, nstep+1):

        ie_barrel = _electron_barrel_ - scan_width + 2*ie*scan_width/nstep
        ie_endcap = _electron_endcap_ - scan_width + 2*ie*scan_width/nstep

        im_barrel = _muon_barrel_ - scan_width + 2*im*scan_width/nstep
        im_endcap = _muon_endcap_ - scan_width + 2*im*scan_width/nstep

        print 'scanning (ie, im) = (', ie, ',', im, ') : th for (e_barrel, e_endcap, mu_barrel, mu_endcap) = ', ie_barrel, ie_endcap, im_barrel, im_endcap

        
        dict_counter = {
            'Diboson':{'process':['WW','WZ','ZZ'], 'counter':0.},
#            'ttbar':{'process':['tt1l','tt2l'], 'counter':0.},
            'ttV':{'process':['TTW','TTZ'], 'counter':0.},
            'ttH':{'process':['TTH'], 'counter':0.},
            'signal':{'process':['tH_YtMinus1'], 'counter':0.}
            }

        for rindex, iregion in enumerate(region):

            if iregion is not 'signal':
                continue

            for index, iprocess in enumerate(process):

                if iprocess=='data':
                    continue

                fname = directory + '/f12_' + iregion + '_' + iprocess + '.root'
#                print fname
                myfile = TFile(fname)
                main = gDirectory.Get('Tree')
                
                evt_dict = {}

                for jentry in xrange(main.GetEntries()):
                    
                    ientry = main.LoadTree(jentry)
                    nb = main.GetEntry(jentry)
                    
                    if not (evaluate(main.muon_eta, main.muon_new_mva, main.muon_id, im_barrel, im_endcap, 'muon', iregion)==True and \
                            evaluate(main.electron_eta, main.electron_new_mva, main.electron_id, ie_barrel, ie_endcap, 'electron', iregion)==True):
                        continue


                    sumpt = main.muon_pt + main.electron_pt + main.tau_pt

                    dict = {'sumpt':sumpt,
                            'weight':main.evt_weight}
            
                    if evt_dict.has_key(main.evt_evt):
                        if evt_dict[main.evt_evt]['sumpt'] < sumpt:
                            evt_dict[main.evt_evt] = dict
                    else:
                        evt_dict[main.evt_evt] = dict


                for key, value in dict_counter.iteritems():
                    if iprocess in value['process']:
                        for ekey, evalue in evt_dict.iteritems():
                            dict_counter[key]['counter'] += float(evalue['weight'])





                
        nevent_mc = [0,0,0,0]
        nevent_data = [0,0,0,0]
        nsf = [0,0,0,0]

        for rindex, iregion in enumerate(region):

            if iregion is 'signal':
                continue

            for index, iprocess in enumerate(process):
        
#                if iprocess in ['WW', 'WZ', 'ZZ', 'tt0l', 'tt1l', 'tt2l', 'TTW', 'TTZ', 'TTH','data']:
                if iprocess in ['WW', 'WZ', 'ZZ', 'TTW', 'TTZ', 'TTH','data']:
                    pass
                else:
                    continue
        
                fname = directory + '/f12_' + iregion + '_' + iprocess + '.root'
                myfile = TFile(fname)
                main = gDirectory.Get('Tree')


                evt_dict = {}
                for jentry in xrange(main.GetEntries()):

                    ientry = main.LoadTree(jentry)
                    nb = main.GetEntry(jentry)

                    
                    if not (evaluate(main.muon_eta, main.muon_new_mva, main.muon_id, im_barrel, im_endcap, 'muon', iregion)==True and \
                            evaluate(main.electron_eta, main.electron_new_mva, main.electron_id, ie_barrel, ie_endcap, 'electron', iregion)==True):
                        continue


                    sumpt = main.muon_pt + main.electron_pt + main.tau_pt
                    
                    dict = {'sumpt':sumpt,
                            'weight':main.evt_weight,
                            'event':main.evt_evt}
                    
                    if evt_dict.has_key(main.evt_evt):
                        #                print 'old = ', evt_dict[main.evt_evt]['sumpt'], 'new = ', sumpt
                        if evt_dict[main.evt_evt]['sumpt'] < sumpt:
                            evt_dict[main.evt_evt] = dict
                    else:
                        evt_dict[main.evt_evt] = dict


                total = 0.        

                for ekey, evalue in evt_dict.iteritems():
                    total += float(evalue['weight'])


                if iprocess=='data':
                    nevent_data[rindex] += total
                else:
                    nevent_mc[rindex] += total

            sf = 1.
            if nevent_data[rindex]==0:
                print '[WARNING] 0 division = ', rindex, iregion, '(Data, MC) = ', nevent_data[rindex], nevent_mc[rindex]
            else:
                sf = Double((nevent_data[rindex] - nevent_mc[rindex])/(nevent_data[rindex]))

            nsf[rindex] = sf




        Nred = 0.

        for rindex, iregion in enumerate(region):
            if iregion is 'signal':
                continue

            for index, iprocess in enumerate(process):

                if iprocess is not 'data':
                    continue

                fname = directory + '/f12_' + iregion + '_' + iprocess + '.root'
                myfile = TFile(fname)
                main = gDirectory.Get('Tree')      

                evt_dict = {}

                for jentry in xrange(main.GetEntries()):
                    
                    ientry = main.LoadTree(jentry)
                    nb = main.GetEntry(jentry)
                    
                    if not (evaluate(main.muon_eta, main.muon_new_mva, main.muon_id, im_barrel, im_endcap, 'muon', iregion)==True and \
                            evaluate(main.electron_eta, main.electron_new_mva, main.electron_id, ie_barrel, ie_endcap, 'electron', iregion)==True):
                        continue

                    sumpt = main.muon_pt + main.electron_pt + main.tau_pt

                    weight_muon = 0.5
                    weight_electron = 0.5

                    if iregion=='antiMu' or iregion=='antiEMu':
                        
                        mvar_map['lepton_pt'][0] = main.muon_pt
#                        mvar_map['lepton_kNN_jetpt'][0] = main.muon_kNN_jetpt
                        mvar_map['evt_njet'][0] = main.evt_njet + 1
                        
                        mvaname = 'muon_' + iprocess
                        
                        weight_muon = muonreader[index].EvaluateMVA(mvaname)
                
                    if iregion=='antiE' or iregion=='antiEMu':
                    
                        evar_map['lepton_pt'][0] = main.electron_pt
#                        evar_map['lepton_kNN_jetpt'][0] = main.electron_kNN_jetpt
                        evar_map['evt_njet'][0] = main.evt_njet + 1
                        
                        mvaname = 'electron_' + iprocess
                        weight_electron = electronreader[index].EvaluateMVA(mvaname)

               
                    kNN_weight = returnkNN(iregion, iprocess, weight_electron, weight_muon)

                    weight_total = main.evt_weight*kNN_weight*nsf[rindex]
                    if iregion=='antiEMu':
                        weight_total *= -1.


                    dict = {'sumpt':sumpt,
                            'weight':weight_total}
            
                    if evt_dict.has_key(main.evt_evt):
                        if evt_dict[main.evt_evt]['sumpt'] < sumpt:
                            evt_dict[main.evt_evt] = dict
                    else:
                        evt_dict[main.evt_evt] = dict

                total_weight = 0
                for ekey, evalue in evt_dict.iteritems():
                    Nred += float(evalue['weight'])
                    total_weight += float(evalue['weight'])
            



        Nsig = 0.
        Nbkg = 0.

        print
        
        for key, value in dict_counter.iteritems():
#            print '\t', key, ('%.2f' % (float(dict_counter[key]['counter'])))

            if key=='signal':
                Nsig += float(dict_counter[key]['counter'])
            else:
                Nbkg += float(dict_counter[key]['counter'])

        s = float(Nsig)*0.919
        b = float(Nbkg+Nred)
        _sig_ = s/math.sqrt(s+b+0.1*b*b)
        print '\t Nsig, Nbkg, Nred, Nbkg_tot, sig. = ', ('%.2f' % float(s)), ('%.2f' % float(Nbkg)), ('%.2f' % float(Nred)), ('%.2f' % float(b)), ('%.2f' % float(_sig_))

#        dict_counter = {
#            'Diboson':{'process':['WW','WZ','ZZ'], 'counter':0.},
#            'ttbar':{'process':['tt1l','tt2l'], 'counter':0.},
#            'ttV':{'process':['TTW','TTZ'], 'counter':0.},
#            'ttH':{'process':['TTH'], 'counter':0.},
#            'signal':{'process':['tH_YtMinus1'], 'counter':0.}
#            }


        muon_b_mva [0] = im_barrel
        muon_e_mva [0] = im_endcap
        electron_b_mva [0] = ie_barrel
        electron_e_mva [0] = ie_endcap
        nsig[0] = s
        nbkg[0] = b
#        nbkg_ttbar[0] = float(dict_counter['ttbar']['counter'])
        nbkg_diboson[0] = float(dict_counter['Diboson']['counter'])
        nbkg_ttv[0] = float(dict_counter['ttV']['counter'])
        nbkg_tth[0] = float(dict_counter['ttH']['counter'])
        nbkg_red[0] = Nred
        significance[0] = _sig_
        n_ie[0] = ie
        n_im[0] = im
        n_max[0] = nstep

        t.Fill()

file.Write()
file.Close()
