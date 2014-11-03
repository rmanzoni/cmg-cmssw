import copy, os
import CMGTools.RootTools.fwlite.Config as cfg
from CMGTools.RootTools.fwlite.Config import printComps

from CMGTools.H2TauTau.triggerMap import pathsAndFilters
from CMGTools.RootTools.RootTools import * 

jobmode = True

# Andrew Summer 13 (MC is identical to the previous one)
puFileMC   = '/afs/cern.ch/user/a/agilbert/public/HTT_Pileup/13-09-13/MC_Summer12_PU_S10-600bins.root'
puFileData = '/afs/cern.ch/user/a/agilbert/public/HTT_Pileup/13-09-13/Data_Pileup_2012_ReRecoPixel-600bins.root'

mc_vertexWeight = None

mc_tauEffWeight_mc = 'effTau_muTau_MC_2012ABCDSummer13'
mc_muEffWeight_mc  = 'effMu_muTau_MC_2012ABCD'
mc_tauEffWeight    = 'effTau_muTau_Data_2012ABCDSummer13'
mc_muEffWeight     = 'effMu_muTau_Data_2012ABCDSummer13'

eventSelector = cfg.Analyzer(
    'EventSelector',
    toSelect = []
    )

jsonAna = cfg.Analyzer(
    'JSONAnalyzer',
    )

triggerAna = cfg.Analyzer(
    'DiTriggerAnalyzer'    
    )

vertexAna = cfg.Analyzer(
    'VertexAnalyzer',
    goodVertices = 'goodPVFilter'  ,
    vertexWeight = mc_vertexWeight ,
    fixedWeight  = 1               ,
    verbose      = False           ,
    )

pileUpAna = cfg.Analyzer(
    'PileUpAnalyzer',
    true = True
    )


EMuTauAna = cfg.Analyzer(
    'WHEMTAnalyzer',
    triggerMap = pathsAndFilters,
    )

genTopAna = cfg.Analyzer(
    'GenTopAnalyzer',
    src     = 'genParticlesPruned',
    verbose = False
    )

tauDecayModeWeighter = cfg.Analyzer(
    'TauDecayModeWeighter',
    )

tauFakeRateWeighter = cfg.Analyzer(
    'TauFakeRateWeighter'
    )

# defined for vbfAna and eventSorter
vbfKwargs = dict( Mjj = 500, deltaEta = 3.5 )

jetAna = cfg.Analyzer(
    'JetAnalyzerEMT',
    jetCol     = 'cmgPFJetSel' ,
    jetPt      = 20.           ,
    jetEta     = 4.7           ,
    btagSFseed = 123456        ,
    relaxJetId = False         , 
    jerCorr    = False         , # jet energy resolution
    jesCorr    = 1.         , # jet energy scale in units of sigma
    )

vbfSimpleAna = cfg.Analyzer(
    'VBFSimpleAnalyzer',
    vbfMvaWeights = ''  ,
    cjvPtCut      = 30. ,
    **vbfKwargs    
    )

treeProducer = cfg.Analyzer(
    'H2TauTauTreeProducerEMT2'
    )

#########################################################################################
# sample definition
from CMGTools.H2TauTau.proto.samples.run2012.emuTau_YutaFeb12 import * 
#########################################################################################

## RICCARDO
MC_list = copy.copy( mc_diboson )
MC_list.extend( mc_ttv )
# MC_list.extend( mc_tH  )
MC_list.extend( mc_ttbarh )

# MC_list = copy.copy( mc_tH  )

# import pdb ; pdb.set_trace()
for mc in MC_list:
    mc.puFileMC = puFileMC
    mc.puFileData = puFileData

selectedComponents = []

sequence = cfg.Sequence([
  #eventSelector,
  jsonAna, 
  triggerAna, # remove for the signal
  vertexAna,
  EMuTauAna,
  jetAna,
  pileUpAna,
  genTopAna,
  treeProducer,
] )

selectedComponents = [comp for comp in selectedComponents if comp.dataset_entries > 0]

if jobmode :
  #selectedComponents = allsamples
  selectedComponents = MC_list
else       :
  comp = tH_YtMinus1
  selectedComponents = [comp]
  for comp in selectedComponents :
    comp.splitFactor = 1
    
config = cfg.Config( components = selectedComponents,
                     sequence = sequence )

printComps(config.components, True)
