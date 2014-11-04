import ROOT
import optparse, os
import array

### For options
parser = optparse.OptionParser()
parser.add_option('--channel', action="store", dest="channel", default='electron')
parser.add_option('--region', action="store", dest="region", default='barrel')
parser.add_option('--nbjets', action="store", dest="nbjets", default=-1)
options, args = parser.parse_args()

print 'channel = ', options.channel
print 'region = ', options.region

min_n_bjets = options.nbjets if options.nbjets >= 0 else 0
if options.channel == 'muon' and options.nbjets < 0:
    min_n_bjets = 1


min_n_bjets = int(min_n_bjets)

print 'Requiring at least', min_n_bjets, 'b jets'

TMVA_tools = ROOT.TMVA.Tools.Instance()

lname = options.channel
fname = 'Wjet_' + lname + '_training.root'

nNeighbours = 50

dir="/afs/cern.ch/user/y/ytakahas/work/th_analysis/CMSSW_5_3_14_patch2/src/CMGTools/H2TauTau/th_analysis/EMuTau/analysis/root_aux/"

file_data = ROOT.TFile(dir + fname)
tree_data = file_data.Get('kNNTrainingTree')

training_vars = ['lepton_pt', 'evt_njet']

baseline_selection = 'evt_nbjet>={min_n_bjets}&&(!evt_isMC || evt_id==0 || evt_id==1 || evt_id==24 || evt_id==25)'.format(min_n_bjets=min_n_bjets)
if options.region=='barrel':
    baseline_selection += '&& abs(lepton_eta) < 1.479'
elif options.region == 'endcap':
    baseline_selection += '&& abs(lepton_eta) > 1.479'
else:
    print 'Invalid regional selection !'

print 'baseline selection, ', baseline_selection
    
signal_selection = '(lepton_id > 0.5 && lepton_mva > lepton_mva_threshold)'
background_selection = '!' + signal_selection #(!lepton_iso || !lepton_id)'
# signal_selection = '(lepton_mva > lepton_mva_threshold)'
# background_selection = '!' + signal_selection #(!lepton_iso || !lepton_id)'

signal_selection += '&&' + baseline_selection
background_selection += '&&' + baseline_selection

num_pass = tree_data.GetEntries(signal_selection)
num_fail = tree_data.GetEntries(background_selection)

outFile = ROOT.TFile('TMVA.root', 'recreate')

factory    = ROOT.TMVA.Factory(
    "TMVAClassification", 
    outFile, 
    "!V:!Silent:Color:DrawProgressBar:Transformations=I" ) 

for var in training_vars:
    factory.AddVariable(var, 'F') 

factory.SetWeightExpression('evt_weight*evt_isMCw')


factory.SetInputTrees(tree_data, ROOT.TCut(signal_selection), ROOT.TCut(background_selection))
factory.PrepareTrainingAndTestTree( ROOT.TCut(''), ROOT.TCut(''),
                                    "nTrain_Signal={num_pass}:nTrain_Background={num_fail}:SplitMode=Random:NormMode=None:!V" )

factory.BookMethod( 
    ROOT.TMVA.Types.kKNN, "KNN50", 
    "H:nkNN={nNeighbours}:ScaleFrac=0.:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T".format(nNeighbours=nNeighbours))


factory.TrainAllMethods()

os.system('cp ' + dir + fname +' data/Wjet_' + lname + '_training_data_' + options.region + '_knn.root')

reader = ROOT.TMVA.Reader("!Color:Silent=T:Verbose=F")

var_list = []
for var in training_vars:
    var_list.append(array.array('f',[0]))
    reader.AddVariable(var, var_list[-1])

reader.BookMVA('KNN50', 'weights/TMVAClassification_KNN50.weights.xml')




var_dict = {
    'lepton_pt':{'nbins':30, 'xmin':0., 'xmax':100.},
    'lepton_eta':{'nbins':30, 'xmin':-2.5, 'xmax':2.5},
    'lepton_kNN_jetpt':{'nbins':30, 'xmin':0., 'xmax':100.},
    'evt_Mem':{'nbins':30, 'xmin':0., 'xmax':400.},
    'slepton_pt':{'nbins':30, 'xmin':0., 'xmax':100.},
    'slepton_eta':{'nbins':30, 'xmin':-2.5, 'xmax':2.5},
    'evt_njet':{'nbins':10, 'xmin':-0.5, 'xmax':9.5},
}

for var in var_dict:
    vd = var_dict[var]
    vd['hist_w'] = ROOT.TH1F(var+'_w', '', vd['nbins'], vd['xmin'], vd['xmax'])
    vd['hist_p'] = ROOT.TH1F(var+'_p', '', vd['nbins'], vd['xmin'], vd['xmax'])


n_signal = 0
n_background = 0

for evt in tree_data:
    for var, arr in zip(training_vars, var_list):
        arr[0] = getattr(evt, var)

    mva_val = reader.EvaluateMVA('KNN50')
    # print mva_val
    if evt.evt_nbjet>=min_n_bjets:# and evt.lepton_id > 0.5:
        if not evt.evt_isMC and ((options.region=='barrel' and abs(evt.lepton_eta) < 1.479) or (options.region=='endcap' and abs(evt.lepton_eta) > 1.479)):
            if evt.lepton_id > 0.5 and evt.lepton_mva > evt.lepton_mva_threshold :
                n_signal += 1
                for var in var_dict:
                    vd = var_dict[var]
                    vd['hist_p'].Fill(getattr(evt, var))

            else:
                n_background += 1
                for var in var_dict:
                    vd = var_dict[var]
                    vd['hist_w'].Fill(getattr(evt, var), mva_val/(1.-mva_val))

print 'nS', n_signal
print 'nB', n_background
print 'eff =', float(n_signal)/n_background

cv = ROOT.TCanvas()
for var in var_dict:
    vd = var_dict[var]
    if vd['hist_w'].GetMaximum() > vd['hist_p'].GetMaximum():
        vd['hist_w'].Draw('hist e')
        vd['hist_p'].SetLineColor(2)
        vd['hist_p'].Draw('same hist e')
    else:
        vd['hist_p'].SetLineColor(2)
        vd['hist_p'].Draw('hist e')
        vd['hist_w'].Draw('same hist e')
    print 'Integral ori', vd['hist_p'].Integral()
    print 'Integral pre', vd['hist_w'].Integral()
    cv.Print('plots/' + lname + '_' + var+ '_' + options.region + '.gif')




os.system('mv weights/TMVAClassification_KNN50.weights.xml weights/KNN_data_' + lname + '_' + options.region + '_50.xml')

