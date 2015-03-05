import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.HeppyCore.framework.config import printComps

# Tau-tau analyzers
from CMGTools.H2TauTau.proto.analyzers.TauMuAnalyzer import TauMuAnalyzer
from CMGTools.H2TauTau.proto.analyzers.H2TauTauTreeProducerTauMu import H2TauTauTreeProducerTauMu
from CMGTools.H2TauTau.proto.analyzers.TauDecayModeWeighter import TauDecayModeWeighter
from CMGTools.H2TauTau.proto.analyzers.TauFakeRateWeighter import TauFakeRateWeighter
from CMGTools.H2TauTau.proto.analyzers.LeptonWeighter import LeptonWeighter
from CMGTools.H2TauTau.proto.analyzers.SVfitProducer import SVfitProducer

# common configuration and sequence
from CMGTools.H2TauTau.htt_ntuple_base_cff import commonSequence, genAna, dyJetsFakeAna, puFileData, puFileMC, eventSelector


### mu-tau specific configuration settings

# 'Nom', 'Up', 'Down', or None
shift = None
syncntuple = True
computeSVfit = True

# When ready, include weights from CMGTools.H2TauTau.proto.weights.weighttable

# mc_tauEffWeight_mc = 'effTau_muTau_MC_2012ABCDSummer13'
# mc_muEffWeight_mc = 'effMu_muTau_MC_2012ABCD'
# mc_tauEffWeight = 'effTau_muTau_Data_2012ABCDSummer13'
# mc_muEffWeight = 'effMu_muTau_Data_2012ABCDSummer13'

mc_tauEffWeight_mc = None
mc_muEffWeight_mc = None
mc_tauEffWeight = None
mc_muEffWeight = None

dyJetsFakeAna.channel = 'mt'

### Define mu-tau specific modules

tauMuAna = cfg.Analyzer(
    TauMuAnalyzer,
    name='TauMuAnalyzer',
    pt1 = 20,
    eta1 = 2.3,
    iso1 = None,
    pt2 = 20,
    eta2 = 2.1,
    iso2 = 0.1,
    m_min = 10,
    m_max = 99999,
    dR_min = 0.5,
    # triggerMap = pathsAndFilters,
    verbose = False
)

tauDecayModeWeighter = cfg.Analyzer(
    TauDecayModeWeighter,
    name='TauDecayModeWeighter',
)

tauFakeRateWeighter = cfg.Analyzer(
    TauFakeRateWeighter,
    name='TauFakeRateWeighter'
)

tauWeighter = cfg.Analyzer(
    LeptonWeighter,
    name='LeptonWeighter_tau',
    effWeight = None,
    effWeightMC = None,
    lepton = 'leg1',
    verbose = False,
    disable = True,
    )

muonWeighter = cfg.Analyzer(
    LeptonWeighter,
    name='LeptonWeighter_mu',
    effWeight = None,
    effWeightMC = None,
    lepton = 'leg2',
    verbose = False,
    disable = True,
    idWeight = None,
    isoWeight = None    
    )

treeProducer = cfg.Analyzer(
    H2TauTauTreeProducerTauMu,
    name='H2TauTauTreeProducerTauMu'
    )

syncTreeProducer = cfg.Analyzer(
    H2TauTauTreeProducerTauMu,
    name='H2TauTauSyncTreeProducerTauMu',
    varStyle='sync',
    skimFunction='event.isSignal'
    )

svfitProducer = cfg.Analyzer(
    SVfitProducer,
    name='SVfitProducer',
    integration='VEGAS',
    #integration='MarkovChain',
    verbose=True,
    #order='21', # muon first, tau second
    l1type='tau',
    l2type='muon'
    )
    
###################################################
### CONNECT SAMPLES TO THEIR ALIASES AND FILES  ###
###################################################
from CMGTools.H2TauTau.proto.samples.phys14.tauMu_Jan_Feb13 import MC_list, mc_dict

###################################################
###              ASSIGN PU to MC                ###
###################################################
for mc in MC_list:
    mc.puFileData = puFileData
    mc.puFileMC = puFileMC

###################################################
###             SET COMPONENTS BY HAND          ###
###################################################
selectedComponents = [mc_dict['HiggsGGH125']]
# for c in selectedComponents : c.splitFactor *= 5

###################################################
###                  SEQUENCE                   ###
###################################################
sequence = commonSequence
sequence.insert(sequence.index(genAna), tauMuAna)
sequence.append(tauDecayModeWeighter)
sequence.append(tauFakeRateWeighter)
sequence.append(tauWeighter)
sequence.append(muonWeighter)
if computeSVfit: 
    sequence.append(svfitProducer)
sequence.append(treeProducer)
if syncntuple:
    sequence.append(syncTreeProducer)

###################################################
###             CHERRY PICK EVENTS              ###
###################################################
eventSelector.toSelect = [
# 9900,
2295,
# 2387,
# 10708,
# 12903,
# 14346,
# 20892,
# 18330,
# 30690,
# 25028,
# 101963,
# 32014,
# 32083,
# 37201,
# 37231,
# 37331,
# 40237,
# 88452,
# 108624,
# 108916,
# 157183,
# 110384,
# 149769,
# 118933,
# 124979,
# 120862,
# 225472,
# 156341,
# 154924,
# 233043,
# 172839,
# 199620,
# 211839,
# 198722,
# 196990,
# 226702,
# 232388,
# 449583,
# 441540,
# 443310,
# 443725,
# 468223,
# 488832,
# 447815,
# 476046,
# 481577,
# 460393,
# 468923,
# 482373,
# 492342,
# 491937,
# 496893,
# 73332,
# 151521,
# 151569,
# 169921,
# 169929,
# 170071,
# 176994,
# 206107,
# 234107,
# 247021,
# 243800,
# 258131,
# 248209,
# 248575,
# 294133,
# 290588,
# 262269,
# 335559,
# 335571,
# 289438,
# 289454,
# 291817,
# 342602,
# 413549,
# 453009,
# 304657,
# 297604,
# 342384,
# 330246,
# 338285,
# 332162,
# 341067,
# 350791,
# 343316,
# 378930,
# 406526,
# 353566,
# 393965,
# 409108,
# 373476,
# 375633,
# 392602,
# 392661,
# 386694,
# 382770,
# 397288,
# 389451,
# 399110,
# 413093,
# 401117,
# 422475,
# 419874,
# 439431,
# 423943,
# 423989,
# 441062,
# 441118,
# 455156,
# 438738,
# 439950,
# 230907,
# 235401,
# 237754,
# 275445,
# 256264,
# 265678,
# 265710,
# 276245,
# 278982,
# 284689,
# 318537,
# 285859,
# 286228,
# 295617,
# 301691,
# 306541,
# 305553,
# 313502,
# 313125,
# 315728,
# 317683,
# 7981,
# 10814,
# 10365,
# 25267,
# 49826,
# 61187,
# 49713,
# 57702,
# 117234,
# 113370,
# 118238,
# 118417,
# 120001,
# 120272,
# 120400,
# 123349,
# 124555,
# 124803,
# 200322,
# 216045,
# 199212,
# 237559,
# 253619,
# 264307,
# 274801,
# 266395,
# 274665,
# 296107,
# 278518,
# 314458,
]
sequence.insert(0, eventSelector)

###################################################
###            SET BATCH OR LOCAL               ###
###################################################
# JAN - can we finally get this via command line options?
test = 1  # test = 0 run on batch, test = 1 run locally
if test == 1:
    comp = mc_dict['HiggsGGH125']
    selectedComponents = [comp]
    comp.splitFactor = 1
    # comp.files = comp.files[:1]


# the following is declared in case this cfg is used in input to the
# heppy.py script
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
config = cfg.Config(components=selectedComponents,
                    sequence=sequence,
                    services=[],
                    events_class=Events
                    )

printComps(config.components, True)

def modCfgForPlot(config):
    config.components = []
