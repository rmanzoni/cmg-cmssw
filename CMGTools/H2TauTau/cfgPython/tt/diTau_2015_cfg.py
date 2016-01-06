import os
import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.HeppyCore.framework.config import printComps
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption

# Tau-tau analyzers
from CMGTools.H2TauTau.proto.analyzers.TauTauAnalyzer             import TauTauAnalyzer
from CMGTools.H2TauTau.proto.analyzers.H2TauTauTreeProducerTauTau import H2TauTauTreeProducerTauTau
from CMGTools.H2TauTau.proto.analyzers.TauDecayModeWeighter       import TauDecayModeWeighter
from CMGTools.H2TauTau.proto.analyzers.LeptonWeighter             import LeptonWeighter
from CMGTools.H2TauTau.proto.analyzers.TauP4Scaler                import TauP4Scaler
from CMGTools.H2TauTau.proto.analyzers.SVfitProducer              import SVfitProducer

# common configuration and sequence
from CMGTools.H2TauTau.htt_ntuple_base_cff import commonSequence, genAna, dyJetsFakeAna, puFileData, puFileMC, eventSelector

# Get all heppy options; set via '-o production' or '-o production=True'

# production = True run on batch, production = False (or unset) run locally
production = getHeppyOption('production')
# production = True

# local switches
syncntuple   = False
computeSVfit = True
pick_events  = False

dyJetsFakeAna.channel = 'tt'

### Define tau-tau specific modules

tauTauAna = cfg.Analyzer(
  class_object        = TauTauAnalyzer                            ,
  name                = 'TauTauAnalyzer'                          ,
  pt1                 = 45                                        ,
  eta1                = 2.1                                       ,
  iso1                = 1.                                        ,
  looseiso1           = 999999999.                                ,
  pt2                 = 45                                        ,
  eta2                = 2.1                                       ,
  iso2                = 1.                                        ,
  looseiso2           = 999999999.                                ,
  isolation           = 'byCombinedIsolationDeltaBetaCorrRaw3Hits',
  m_min               = 10                                        ,
  m_max               = 99999                                     ,
  dR_min              = 0.5                                       ,
  jetPt               = 30.                                       ,
  jetEta              = 4.7                                       ,
  relaxJetId          = False                                     ,
  verbose             = False                                     ,
  from_single_objects = False                                     ,
  )

tau1Calibration = cfg.Analyzer(
  TauP4Scaler       ,
  'TauP4Scaler_tau1',
  leg      = 'leg1' ,
  scaleMET = False  ,
  verbose  = False  ,
  )

tau2Calibration = cfg.Analyzer(
  TauP4Scaler       ,
  'TauP4Scaler_tau2',
  leg      = 'leg2' ,
  scaleMET = False  ,
  verbose  = False  ,
  )

tauDecayModeWeighter = cfg.Analyzer(
  TauDecayModeWeighter   ,
  'TauDecayModeWeighter' ,
  legs = ['leg1', 'leg2'],
  )

tau1Weighter = cfg.Analyzer(
  LeptonWeighter                    ,
  name        ='LeptonWeighter_tau1',
  effWeight   = None                ,
  effWeightMC = None                ,
  lepton      = 'leg1'              ,
  verbose     = False               ,
  disable     = True                ,
  )

tau2Weighter = cfg.Analyzer(
  LeptonWeighter                    ,
  name        ='LeptonWeighter_tau2',
  effWeight   = None                ,
  effWeightMC = None                ,
  lepton      = 'leg2'              ,
  verbose     = False               ,
  disable     = True                ,
  )

treeProducer = cfg.Analyzer(
  H2TauTauTreeProducerTauTau         ,
  name = 'H2TauTauTreeProducerTauTau'
  )

syncTreeProducer = cfg.Analyzer(
  H2TauTauTreeProducerTauTau                     ,
  name         = 'H2TauTauSyncTreeProducerTauTau',
  varStyle     = 'sync'                          ,
  #skimFunction = 'event.isSignal' #don't cut out any events from the sync tuple
  )

svfitProducer = cfg.Analyzer(
  SVfitProducer,
  name                       = 'SVfitProducer',
  integration                = 'MarkovChain'  , # 'VEGAS'
  integrateOverVisPtResponse = True           ,
#   visPtResponseFile          = os.environ['CMSSW_BASE']+'/src/TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root',
  visPtResponseFile          = os.environ['CMSSW_BASE']+'/src/TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF_new.root',
  verbose                    = False          ,
  l1type                     = 'tau'          ,
  l2type                     = 'tau'
  )

###################################################
### CONNECT SAMPLES TO THEIR ALIASES AND FILES  ###
# ###################################################

from CMGTools.H2TauTau.proto.samples.spring15.components_25ns import *

creator = ComponentCreator()

ggh125 = creator.makeMCComponent(
    'GGH125', 
    '/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM', 
    'CMS', 
    '.*root', 
    1.0
)

ggh125.files = ['file:/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_7_4_3/src/CMGTools/H2TauTau/prod/diTau_fullsel_tree_CMG.root']

MC_list = [ggh125]

split_factor = 1e5

for sample in MC_list:
    sample.triggers = mc_triggers_tt
    sample.splitFactor = splitFactor(sample, split_factor)

###################################################
###              ASSIGN PU to MC                ###
###################################################
for mc in MC_list:
    mc.puFileData = puFileData
    mc.puFileMC = puFileMC

###################################################
###             SET COMPONENTS BY HAND          ###
###################################################
selectedComponents = MC_list
# selectedComponents = mc_dict['HiggsGGH125']
# for c in selectedComponents : c.splitFactor *= 5

###################################################
###                  SEQUENCE                   ###
###################################################
sequence = commonSequence
sequence.insert(sequence.index(genAna), tauTauAna)
sequence.append(tau1Calibration)
sequence.append(tau2Calibration)
sequence.append(tauDecayModeWeighter)
sequence.append(tau1Weighter)
sequence.append(tau2Weighter)
if computeSVfit:
    sequence.append(svfitProducer)
sequence.append(treeProducer)
if syncntuple:
    sequence.append(syncTreeProducer)

###################################################
###             CHERRY PICK EVENTS              ###
###################################################
if pick_events:

    import csv
    fileName = '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_7_4_3/src/CMGTools/H2TauTau/cfgPython/2015-sync/Imperial.csv'
#     fileName = '/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_7_4_3/src/CMGTools/H2TauTau/cfgPython/2015-sync/CERN.csv'
    f = open(fileName, 'rb')
    reader = csv.reader(f)
    evtsToPick = []

    for i, row in enumerate(reader):
        evtsToPick += [int(j) for j in row]

    eventSelector.toSelect = evtsToPick
    sequence.insert(0, eventSelector)

###################################################
###            SET BATCH OR LOCAL               ###
###################################################
if not production:
  cache                = True
#   comp                 = my_connect.mc_dict['HiggsGGH125']
  comp                 = ggh125
  selectedComponents   = [comp]
  comp.splitFactor     = 1
  comp.fineSplitFactor = 1
#   comp.files           = comp.files[:1]
#   for comp in selectedComponents:
#     comp.splitFactor     = 1
#     comp.fineSplitFactor = 4
    


# the following is declared in case this cfg is used in input to the
# heppy.py script
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
config = cfg.Config( components   = selectedComponents,
                     sequence     = sequence          ,
                     services     = []                ,
                     events_class = Events
                     )

printComps(config.components, True)

def modCfgForPlot(config):
  config.components = []
