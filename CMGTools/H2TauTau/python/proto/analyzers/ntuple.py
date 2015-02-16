#!/bin/env python

from PhysicsTools.Heppy.analyzers.core.autovars     import NTupleVariable, NTupleCollection, NTupleObjectType
from PhysicsTools.Heppy.analyzers.objects.autophobj import jetType
from PhysicsTools.HeppyCore.utils.deltar            import deltaR, deltaPhi
from CMGTools.H2TauTau.proto.analyzers.tauIDs       import tauIDs

def multi_getattr(obj, attr, ifNonePass = True, **kw):
  '''
  Get a named attribute from an object; multi_getattr(x, 'a.b.c.d') is
  equivalent to x.a.b.c.d. When a default argument is given, it is
  returned when any attribute in the chain doesn't exist; without
  it, an exception is raised when a missing attribute is encountered.
  '''
  attributes = attr.split('.')
  import pdb ; pdb.set_trace()
  for i in attributes:
    i = i.replace('()','')
    try:
      #import pdb ; pdb.set_trace()
      obj = getattr(obj, i)
      if callable(obj):
        obj = obj()
        objToReturn = obj
    except AttributeError:
      if kw.has_key('default'):
        return kw['default']
      elif not ifNonePass:
        raise
      else :
        objToReturn = None
  return obj
  
def pthiggs(event) :
  return (event.diLepton.leg1().p4()+event.diLepton.leg2().p4()+event.diLepton.met().p4()).pt()
def l1eta  (event) :
  return event.diLepton.leg1().eta()
def l2eta  (event) :
  return event.diLepton.leg2().eta()
def l1phi  (event) :
  return event.diLepton.leg1().phi()
def l2phi  (event) :
  return event.diLepton.leg2().phi()
def metphi (event) :
  return event.diLepton.met().phi()

# di-tau
diLeptonVars = [
  NTupleVariable('visMass'      , lambda event : event.diLepton.mass()                                         , float, mcOnly=False, help='di-tau visible mass'                                                                         ),
  NTupleVariable('svfitMass'    , lambda event : event.diLepton.svfitMass()                                    , float, mcOnly=False, help='di-tau SVfit mass'                                                                           ),
  NTupleVariable('pZetaMET'     , lambda event : event.diLepton.pZetaMET()                                     , float, mcOnly=False, help='p_Z invisible'                                                                               ),
  NTupleVariable('pZetaVis'     , lambda event : event.diLepton.pZetaVis()                                     , float, mcOnly=False, help='p_Z visible'                                                                                 ),
  NTupleVariable('pZetaDisc'    , lambda event : event.diLepton.pZetaDisc()                                    , float, mcOnly=False, help='p_Z discriminator'                                                                           ),
  NTupleVariable('mt'           , lambda event : event.diLepton.mTLeg2()                                       , float, mcOnly=False, help='m_T(met,leg2). Leg2 usually means the lepton. Met is the one associated to the di-tau pair'  ),
  NTupleVariable('mtleg2'       , lambda event : event.diLepton.mTLeg2()                                       , float, mcOnly=False, help='m_T(met,leg2). Met is the one associated to the di-tau pair'                                 ),
  NTupleVariable('mtleg1'       , lambda event : event.diLepton.mTLeg1()                                       , float, mcOnly=False, help='m_T(met,leg1). Met is the one associated to the di-tau pair'                                 ),
  NTupleVariable('met'          , lambda event : event.diLepton.met().pt()                                     , float, mcOnly=False, help='met pt'                                                                                      ),
  NTupleVariable('mex'          , lambda event : event.diLepton.met().px()                                     , float, mcOnly=False, help='mex'                                                                                         ),
  NTupleVariable('mey'          , lambda event : event.diLepton.met().py()                                     , float, mcOnly=False, help='mey'                                                                                         ),
  NTupleVariable('metphi'       , lambda event : event.diLepton.met().phi()                                    , float, mcOnly=False, help='met phi'                                                                                     ),
  # RIC: these don't appear to be there FIXME
  # NTupleVariable('metcov00'     , lambda event : event.diLepton.metSig().significance()(0,0))                  , float, mcOnly=False, help='met covariance matrix (0,0)'                                                                 ),
  # NTupleVariable('metcov01'     , lambda event : event.diLepton.metSig().significance()(0,1))                  , float, mcOnly=False, help='met covariance matrix (0,1)'                                                                 ),
  # NTupleVariable('metcov10'     , lambda event : event.diLepton.metSig().significance()(1,0))                  , float, mcOnly=False, help='met covariance matrix (1,0)'                                                                 ),
  # NTupleVariable('metcov11'     , lambda event : event.diLepton.metSig().significance()(1,1))                  , float, mcOnly=False, help='met covariance matrix (1,1)'                                                                 ),
  NTupleVariable('pthiggs'      , lambda event : pthiggs(event)                                                , float, mcOnly=False, help='pt Higgs: (leg1 + leg2 + MET) pt'                                                            ),
  NTupleVariable('deltaPhiL1L2' , lambda event : deltaPhi(l1phi(event), l2phi(event))                          , float, mcOnly=False, help='delta phi between leg 1 and leg2'                                                            ),
  NTupleVariable('deltaEtaL1L2' , lambda event : abs(l1eta(event)-l2eta(event))                                , float, mcOnly=False, help='delta eta between leg 1 and leg2'                                                            ),
  NTupleVariable('deltaRL1L2'   , lambda event : deltaR(l1eta(event), l1phi(event), l2eta(event), l2phi(event)), float, mcOnly=False, help='delta R between leg 1 and leg2'                                                              ),
  NTupleVariable('deltaPhiL1MET', lambda event : deltaPhi(l1phi(event), metphi(event))                         , float, mcOnly=False, help='delta phi between leg 1 and met. Met is the one associated to the di-tau pair'               ),
  NTupleVariable('deltaPhiL2MET', lambda event : deltaPhi(l2phi(event), metphi(event))                         , float, mcOnly=False, help='delta phi between leg 2 and met. Met is the one associated to the di-tau pair'               ),
]  

def VBFmjj           (event) :
  if hasattr(event,'vbf') : return event.vbf.mjj 
  else                    : return None
def VBFdeta          (event) :
  if hasattr(event,'vbf') : return event.vbf.deta 
  else                    : return None
def VBFncj           (event) :
  if hasattr(event,'vbf') : return len(event.vbf.centralJets) 
  else                    : return None
def VBFmva           (event) :
  if hasattr(event,'vbf') : return event.vbf.mva 
  else                    : return None
def VBFdphi          (event) :
  if hasattr(event,'vbf') : return event.vbf.dphi 
  else                    : return None
def VBFdijetpt       (event) :
  if hasattr(event,'vbf') : return event.vbf.dijetpt 
  else                    : return None
def VBFdijetphi      (event) :
  if hasattr(event,'vbf') : return event.vbf.dijetphi 
  else                    : return None
def VBFdphidijethiggs(event) :
  if hasattr(event,'vbf') : return event.vbf.dphidijethiggs 
  else                    : return None
def VBFvisjeteta     (event) :
  if hasattr(event,'vbf') : return event.vbf.visjeteta 
  else                    : return None
def VBFptvis         (event) :
  if hasattr(event,'vbf') : return event.vbf.ptvis 
  else                    : return None

# vbf
VBFVars =  [
  NTupleVariable('mjj'      , lambda event : VBFmjj           (event), float, mcOnly=False, help='VBF invariant mass of the two highest pt jets'                        ),
  NTupleVariable('deta'     , lambda event : VBFdeta          (event), float, mcOnly=False, help='VBF eta separation between the two highest pt jets'                   ),
  NTupleVariable('nCentral' , lambda event : VBFncj           (event), float, mcOnly=False, help='VBF number of jets with pt>30 in the eta gap between the tagging jets'),
  NTupleVariable('vbf_mva'  , lambda event : VBFmva           (event), float, mcOnly=False, help='VBF MVA score'                                                        ),
  NTupleVariable('jdphi'    , lambda event : VBFdphi          (event), float, mcOnly=False, help='VBF phi separation between the two highest pt jets'                   ),
  NTupleVariable('dijetpt'  , lambda event : VBFdijetpt       (event), float, mcOnly=False, help='VBF pt of the di-jet system'                                          ),
  NTupleVariable('dijetphi' , lambda event : VBFdijetphi      (event), float, mcOnly=False, help='VBF phi of the di-jet system'                                         ),
  NTupleVariable('hdijetphi', lambda event : VBFdphidijethiggs(event), float, mcOnly=False, help='VBF phi separation between the Higgs system and the di-jet system'    ),
  NTupleVariable('visjeteta', lambda event : VBFvisjeteta     (event), float, mcOnly=False, help='VBF visible jet eta'                                                  ),
  NTupleVariable('ptvis'    , lambda event : VBFptvis         (event), float, mcOnly=False, help='VBF visible pt'                                                       ),
]

# DY fakes
ZttGenVars =  [
  NTupleVariable('isZtt'  , lambda event : event.isZtt  , float, mcOnly=True, help='Z/H -> tautau, fully hadronic'               ),
  NTupleVariable('isZmt'  , lambda event : event.isZmt  , float, mcOnly=True, help='Z/H -> tautau, mutau'                        ),
  NTupleVariable('isZet'  , lambda event : event.isZet  , float, mcOnly=True, help='Z/H -> tautau, etau'                         ),
  NTupleVariable('isZee'  , lambda event : event.isZee  , float, mcOnly=True, help='Z/H -> tautau, ee'                           ),
  NTupleVariable('isZmm'  , lambda event : event.isZmm  , float, mcOnly=True, help='Z/H -> tautau, mm'                           ),
  NTupleVariable('isZem'  , lambda event : event.isZem  , float, mcOnly=True, help='Z/H -> tautau, em'                           ),
  NTupleVariable('isZEE'  , lambda event : event.isZEE  , float, mcOnly=True, help='Z/H -> ee'                                   ),
  NTupleVariable('isZMM'  , lambda event : event.isZMM  , float, mcOnly=True, help='Z/H -> mm'                                   ),
  NTupleVariable('isZLL'  , lambda event : event.isZLL  , float, mcOnly=True, help='Z/H -> ll (ee or mm)'                        ),
  NTupleVariable('isFake' , lambda event : event.isFake , float, mcOnly=True, help='check DYJetsFakeAnalyzer'                    ),
  NTupleVariable('genMass', lambda event : event.genMass, float, mcOnly=True, help='Z/H gen mass'                                ),
  NTupleVariable('genMet' , lambda event : event.genMet , float, mcOnly=True, help='gen MET from sum of gen neutrino momenta'    ),
  NTupleVariable('genMex' , lambda event : event.genMex , float, mcOnly=True, help='gen MEx from sum of gen neutrino momenta'    ),
  NTupleVariable('genMey' , lambda event : event.genMey , float, mcOnly=True, help='gen MEy from sum of gen neutrino momenta'    ),
  NTupleVariable('genMey' , lambda event : event.genMey , float, mcOnly=True, help='gen MET phi from sum of gen neutrino momenta'),
]

# simple particle
def fillParticle( pName, particle ):
    return [
             NTupleVariable('{pName}_pt'    .format(pName=pName), lambda event : multi_getattr(event, particle + '.pt()'    ), float, mcOnly=False, help='{PART} pt'    .format(PART=particle)),
             NTupleVariable('{pName}_eta'   .format(pName=pName), lambda event : multi_getattr(event, particle + '.eta()'   ), float, mcOnly=False, help='{PART} eta'   .format(PART=particle)),
             NTupleVariable('{pName}_phi'   .format(pName=pName), lambda event : multi_getattr(event, particle + '.phi()'   ), float, mcOnly=False, help='{PART} phi'   .format(PART=particle)),
             NTupleVariable('{pName}_charge'.format(pName=pName), lambda event : multi_getattr(event, particle + '.charge()'), int  , mcOnly=False, help='{PART} charge'.format(PART=particle)),
           ]  
    
def fillGenParticle( pName, particle, isTau = False, genVisTau = None ):
    
    if isTau :
        if genVisTau == None : 
            print 'Please specify the event\'s method that indicates \
                   the gen visible tau. It can (should) differ from %s.\
                   Usually is heppy Tau.physObj.genJet()' %(particle)
            raise
        allvars = [
                    NTupleVariable('{pName}_mass'   .format(pName=pName), lambda event : multi_getattr(event, particle  + '.mass()' ), float, mcOnly=True, help='{PART} gen mass'       .format(PART=particle)),
                    NTupleVariable('{pName}_pdgId'  .format(pName=pName), lambda event : multi_getattr(event, particle  + '.pdgId()'), float, mcOnly=True, help='{PART} pdg ID'         .format(PART=particle)),
                    NTupleVariable('{pName}_vis_pt' .format(pName=pName), lambda event : multi_getattr(event, genVisTau + '.pt()'   ), float, mcOnly=True, help='{PART} gen tau vis pt' .format(PART=particle)),
                    NTupleVariable('{pName}_vis_eta'.format(pName=pName), lambda event : multi_getattr(event, genVisTau + '.eta()'  ), float, mcOnly=True, help='{PART} gen tau vis eta'.format(PART=particle)),
                    NTupleVariable('{pName}_vis_phi'.format(pName=pName), lambda event : multi_getattr(event, genVisTau + '.phi()'  ), float, mcOnly=True, help='{PART} gen tau vis phi'.format(PART=particle)),
                  ] 
    else :             
        allvars = [
                    NTupleVariable('{pName}_mass'   .format(pName=pName), lambda event : multi_getattr(event, particle + '.mass()' ), float, mcOnly=True, help='{PART} gen mass'.format(PART=particle)),
                    NTupleVariable('{pName}_pdgId'  .format(pName=pName), lambda event : multi_getattr(event, particle + '.pdgId()'), float, mcOnly=True, help='{PART} pdg ID'  .format(PART=particle)),
                  ]
    return fillParticle(pName, particle) + allvars
   
# lepton
def fillLepton( pName, lepton ):
    return fillParticle(pName, lepton) + \
           [
             NTupleVariable('{pName}_relIso05'      .format(pName=pName), lambda event : multi_getattr(event, lepton + '.relIsoAllChargedDB05()'), float, mcOnly=False, help='{PART} relative isolation, all charged particles, delta beta corrected, radius 0.5'.format(PART=lepton)),
             NTupleVariable('{pName}_dxy'           .format(pName=pName), lambda event : multi_getattr(event, lepton + '.dxy()'                 ), float, mcOnly=False, help='{PART} dxy wrt its own vertex (vertex has been assigned beforehand)'               .format(PART=lepton)),
             NTupleVariable('{pName}_dz'            .format(pName=pName), lambda event : multi_getattr(event, lepton + '.dz()'                  ), float, mcOnly=False, help='{PART} dz wrt its own vertex (vertex has been assigned beforehand)'                .format(PART=lepton)),
             NTupleVariable('{pName}_weight'        .format(pName=pName), lambda event : multi_getattr(event, lepton + '.weight'                ), float, mcOnly=False, help='{PART} weight. Product of all weights'                                             .format(PART=lepton)),
             NTupleVariable('{pName}_triggerWeight' .format(pName=pName), lambda event : multi_getattr(event, lepton + '.triggerWeight'         ), float, mcOnly=False, help='{PART} trigger weight computed as data efficiency divided by MC efficiency'        .format(PART=lepton)),
             NTupleVariable('{pName}_triggerEffData'.format(pName=pName), lambda event : multi_getattr(event, lepton + '.triggerEffData'        ), float, mcOnly=False, help='{PART} trigger efficiency as measured in data'                                     .format(PART=lepton)),
             NTupleVariable('{pName}_triggerEffMC'  .format(pName=pName), lambda event : multi_getattr(event, lepton + '.triggerEffMC'          ), float, mcOnly=False, help='{PART} trigger efficiency as measured in MC'                                       .format(PART=lepton)),
             NTupleVariable('{pName}_recEffWeight'  .format(pName=pName), lambda event : multi_getattr(event, lepton + '.recEffWeight'          ), float, mcOnly=False, help='{PART} reconstruction ID + Iso weight'                                             .format(PART=lepton)),
           ]  

# muon
def fillMuon( pName, muon ):
    return fillLepton(pName, muon) + \
           [
             NTupleVariable('{pName}_looseId'.format(pName=pName), lambda event : multi_getattr(event, muon + '.looseId()'), float, mcOnly=False, help='Muon loose ID as defined in Heppy Muon'),
             NTupleVariable('{pName}_tightId'.format(pName=pName), lambda event : multi_getattr(event, muon + '.tightId()'), float, mcOnly=False, help='Muon tight ID as defined in Heppy Muon'),
           ] 
    # JAN FIXME - do we need the MVA iso and does it exist?
    # fill(tree, '{pName}_mvaIso'.format(pName=pName), muon.mvaIso() )

# electron
def fillEle( pName, ele ):
    return fillLepton(pName, ele) + \
           [
             NTupleVariable('{pName}_mvaTrigV0'          .format(pName=pName), lambda event : multi_getattr(event, ele + '.mvaTrigV0()'         ), float, mcOnly=False, help='Electron mvaTrigV0 ID as defined in Heppy Muon'         ),
             NTupleVariable('{pName}_mvaNonTrigV0'       .format(pName=pName), lambda event : multi_getattr(event, ele + '.mvaNonTrigV0()'      ), float, mcOnly=False, help='Electron mvaNonTrigV0 ID as defined in Heppy Muon'      ),
             NTupleVariable('{pName}_looseId'            .format(pName=pName), lambda event : multi_getattr(event, ele + '.looseIdForEleTau()'  ), float, mcOnly=False, help='Electron looseIdForEleTau ID as defined in Heppy Muon'  ),
             NTupleVariable('{pName}_tightId'            .format(pName=pName), lambda event : multi_getattr(event, ele + '.tightIdForEleTau()'  ), float, mcOnly=False, help='Electron tightIdForEleTau ID as defined in Heppy Muon'  ),
             NTupleVariable('{pName}_numberOfMissingHits'.format(pName=pName), lambda event : multi_getattr(event, ele + '.lostInner()'         ), float, mcOnly=False, help='Electron lostInner ID as defined in Heppy Muon'         ),
             NTupleVariable('{pName}_passConversionVeto' .format(pName=pName), lambda event : multi_getattr(event, ele + '.passConversionVeto()'), float, mcOnly=False, help='Electron passConversionVeto ID as defined in Heppy Muon'),
           ] 
    # JAN FIXME - do we need the MVA iso and does it exist?
    # fill(tree, '{pName}_mvaIso'.format(pName=pName), ele.mvaIso() )
    # fill(tree, '{pName}_mvaTrigV0'.format(pName=pName), ele.electronID('mvaTrigV0') )
    # fill(tree, '{pName}_mvaNonTrigV0'.format(pName=pName), ele.electronID('mvaNonTrigV0') )

# tau 
def fillTau( pName, tau ):

    ids = []
    for id in tauIDs :
      ids.append( NTupleVariable('{pName}_{tauID}'.format(pName=pName, tauID=id), lambda event : multi_getattr(event, tau + '.tauID(id)'), float, mcOnly=False, help='Tau {tauID} discriminator'.format(tauID=id)) )
    return fillLepton(pName, tau) + ids + \
           [
             NTupleVariable('{pName}_decayMode'.format(pName=pName), lambda event : multi_getattr(event, tau + '.decayMode()'), float, mcOnly=False, help='Tau decay mode'                                      ),
             NTupleVariable('{pName}_mass'     .format(pName=pName), lambda event : multi_getattr(event, tau + '.mass()'     ), float, mcOnly=False, help='Tau mass (RIC: isn\'t it the same as the p4.mass()?)'),
             NTupleVariable('{pName}_zImpact'  .format(pName=pName), lambda event : multi_getattr(event, tau + '.zImpact()'  ), float, mcOnly=False, help='Tau z impact parameter'                              ),
           ]

# jet and b-jet collections
cleanJetsCollection  = { 'cleanJets'  : NTupleCollection('Jets20', jetType, 4, help='Collection of clean jets with pt>20 and |eta|<4.7'),   # max 4 jets
                         'cleanJets30': NTupleCollection('Jets'  , jetType, 2, help='Collection of clean jets with pt>30 and |eta|<4.7'), } # max 2 jets
cleanBJetsCollection = { 'cleanBJets' : NTupleCollection('bJets' , jetType, 2, help='Collection of clean b-jets with pt>20 and |eta|<2.4 and medium CSV'), } # max 2 b-jets

# event global variables
evtVars = [
             NTupleVariable('rho'                  , lambda event : event.rho                  , float, mcOnly=False, help='rho'                  ),
             NTupleVariable('NUP'                  , lambda event : event.NUP                  , float, mcOnly=False, help='NUP'                  ),
             NTupleVariable('vertexWeight'         , lambda event : event.vertexWeight         , float, mcOnly=False, help='vertex weight'        ),
             NTupleVariable('nVert'                , lambda event : len(event.vertices)        , float, mcOnly=False, help='number of vertices'   ),
             NTupleVariable('embedWeight'          , lambda event : event.embedWeight          , float, mcOnly=False, help='embedWeight'          ),
             NTupleVariable('hqtWeight'            , lambda event : event.higgsPtWeight        , float, mcOnly=False, help='hqtWeight'            ),
             NTupleVariable('hqtWeightUp'          , lambda event : event.higgsPtWeightUp      , float, mcOnly=False, help='hqtWeightUp'          ),
             NTupleVariable('hqtWeightDown'        , lambda event : event.higgsPtWeightDown    , float, mcOnly=False, help='hqtWeightDown'        ),
             NTupleVariable('NJetWeight'           , lambda event : event.NJetWeight           , float, mcOnly=False, help='NJetWeight'           ),
             NTupleVariable('tauFakeRateWeightUp'  , lambda event : event.tauFakeRateWeightUp  , float, mcOnly=False, help='tauFakeRateWeightUp'  ),
             NTupleVariable('tauFakeRateWeightDown', lambda event : event.tauFakeRateWeightDown, float, mcOnly=False, help='tauFakeRateWeightDown'),
             NTupleVariable('tauFakeRateWeight'    , lambda event : event.tauFakeRateWeight    , float, mcOnly=False, help='tauFakeRateWeight'    ),
           ]

# common flat ntuple content
htt_globalVariables = diLeptonVars + ZttGenVars + VBFVars + evtVars
htt_globalObjects   = {}
htt_collections     = dict(cleanJetsCollection, **cleanBJetsCollection) # http://stackoverflow.com/questions/1551666/how-can-2-python-dictionaries-become-1/1551878#1551878




