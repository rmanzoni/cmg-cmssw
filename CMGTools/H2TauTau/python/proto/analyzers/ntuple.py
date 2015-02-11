#!/bin/env python

# from PhysicsTools.Heppy.analyzers.core.AutoFillTreeProducer  import AutoFillTreeProducer 
from PhysicsTools.Heppy.analyzers.core.autovars              import NTupleVariable 
from PhysicsTools.HeppyCore.utils.deltar                     import deltaR, deltaPhi
from CMGTools.H2TauTau.proto.analyzers.tauIDs                import tauIDs

def multi_getattr(obj, attr, **kw):
  '''
  Get a named attribute from an object; multi_getattr(x, 'a.b.c.d') is
  equivalent to x.a.b.c.d. When a default argument is given, it is
  returned when any attribute in the chain doesn't exist; without
  it, an exception is raised when a missing attribute is encountered.
  '''
  attributes = attr.split(".")
  for i in attributes:
    try:
      obj = getattr(obj, i)
      if callable(obj):
        obj = obj()
    except AttributeError:
      if kw.has_key('default'):
        return kw['default']
      else:
        raise
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

# vbf
VBFVars =  [
  NTupleVariable('mjj'      , lambda event : event.vbf.mjj             , float, mcOnly=False, help='VBF invariant mass of the two highest pt jets'                        ),
  NTupleVariable('deta'     , lambda event : event.vbf.deta            , float, mcOnly=False, help='VBF eta separation between the two highest pt jets'                   ),
  NTupleVariable('nCentral' , lambda event : len(event.vbf.centralJets), float, mcOnly=False, help='VBF number of jets with pt>30 in the eta gap between the tagging jets'),
  NTupleVariable('mva'      , lambda event : event.vbf.mva             , float, mcOnly=False, help='VBF MVA score'                                                        ),
  NTupleVariable('jdphi'    , lambda event : event.vbf.dphi            , float, mcOnly=False, help='VBF phi separation between the two highest pt jets'                   ),
  NTupleVariable('dijetpt'  , lambda event : event.vbf.dijetpt         , float, mcOnly=False, help='VBF pt of the di-jet system'                                          ),
  NTupleVariable('dijetphi' , lambda event : event.vbf.dijetphi        , float, mcOnly=False, help='VBF phi of the di-jet system'                                         ),
  NTupleVariable('hdijetphi', lambda event : event.vbf.dphidijethiggs  , float, mcOnly=False, help='VBF phi separation between the Higgs system and the di-jet system'    ),
  NTupleVariable('visjeteta', lambda event : event.vbf.visjeteta       , float, mcOnly=False, help='VBF visible jet eta'                                                  ),
  NTupleVariable('ptvis'    , lambda event : event.vbf.ptvis           , float, mcOnly=False, help='VBF visible pt'                                                       ),
]

# DY fakes
ZttGenVars =  [
  NTupleVariable('isZtt' , lambda event : event.isZtt , float, mcOnly=True, help='Z/H -> tautau, fully hadronic'),
  NTupleVariable('isZmt' , lambda event : event.isZmt , float, mcOnly=True, help='Z/H -> tautau, mutau'         ),
  NTupleVariable('isZet' , lambda event : event.isZet , float, mcOnly=True, help='Z/H -> tautau, etau'          ),
  NTupleVariable('isZee' , lambda event : event.isZee , float, mcOnly=True, help='Z/H -> tautau, ee'            ),
  NTupleVariable('isZmm' , lambda event : event.isZmm , float, mcOnly=True, help='Z/H -> tautau, mm'            ),
  NTupleVariable('isZem' , lambda event : event.isZem , float, mcOnly=True, help='Z/H -> tautau, em'            ),
  NTupleVariable('isZEE' , lambda event : event.isZEE , float, mcOnly=True, help='Z/H -> ee'                    ),
  NTupleVariable('isZMM' , lambda event : event.isZMM , float, mcOnly=True, help='Z/H -> mm'                    ),
  NTupleVariable('isZLL' , lambda event : event.isZLL , float, mcOnly=True, help='Z/H -> ll (ee or mm)'         ),
  NTupleVariable('isFake', lambda event : event.isFake, float, mcOnly=True, help='check DYJetsFakeAnalyzer'     ),
]


# simple particle
def fillParticle( pName, particle ):
    return [
             NTupleVariable('{pName}_pt'    .format(pName=pName), lambda event : multi_getattr(event, particle).pt()    , float, mcOnly=False, help='{PART} pt'    .format(PART=particle)),
             NTupleVariable('{pName}_eta'   .format(pName=pName), lambda event : multi_getattr(event, particle).eta()   , float, mcOnly=False, help='{PART} eta'   .format(PART=particle)),
             NTupleVariable('{pName}_phi'   .format(pName=pName), lambda event : multi_getattr(event, particle).phi()   , float, mcOnly=False, help='{PART} phi'   .format(PART=particle)),
             NTupleVariable('{pName}_charge'.format(pName=pName), lambda event : multi_getattr(event, particle).charge(), int  , mcOnly=False, help='{PART} charge'.format(PART=particle)),
           ]  
    
def fillGenParticle( pName, particle ):
    return fillParticle(pName, particle) + \
           [
             NTupleVariable('{pName}_mass'  .format(pName=pName), lambda event : multi_getattr(event, particle).mass() , float, mcOnly=True, help='{PART} gen mass'.format(PART=particle)),
             NTupleVariable('{pName}_pdgId' .format(pName=pName), lambda event : multi_getattr(event, particle).pdgId(), float, mcOnly=True, help='{PART} pdg ID'  .format(PART=particle)),
           ]  
   
# lepton
def fillLepton( pName, lepton ):
    return fillParticle(pName, lepton) + \
           [
             NTupleVariable('{pName}_relIso05'      .format(pName=pName), lambda event : multi_getattr(event, lepton).relIsoAllChargedDB05(), float, mcOnly=False, help='{PART} relative isolation, all charged particles, delta beta corrected, radius 0.5'.format(PART=lepton)),
             NTupleVariable('{pName}_dxy'           .format(pName=pName), lambda event : multi_getattr(event, lepton).dxy()                 , float, mcOnly=False, help='{PART} dxy wrt its own vertex (vertex has been assigned beforehand)'               .format(PART=lepton)),
             NTupleVariable('{pName}_dz'            .format(pName=pName), lambda event : multi_getattr(event, lepton).dz()                  , float, mcOnly=False, help='{PART} dz wrt its own vertex (vertex has been assigned beforehand)'                .format(PART=lepton)),
             NTupleVariable('{pName}_weight'        .format(pName=pName), lambda event : multi_getattr(event, lepton).weight                , float, mcOnly=False, help='{PART} weight. Product of all weights'                                             .format(PART=lepton)),
             NTupleVariable('{pName}_triggerWeight' .format(pName=pName), lambda event : multi_getattr(event, lepton).triggerWeight         , float, mcOnly=False, help='{PART} trigger weight computed as data efficiency divided by MC efficiency'        .format(PART=lepton)),
             NTupleVariable('{pName}_triggerEffData'.format(pName=pName), lambda event : multi_getattr(event, lepton).triggerEffData        , float, mcOnly=False, help='{PART} trigger efficiency as measured in data'                                     .format(PART=lepton)),
             NTupleVariable('{pName}_triggerEffMC'  .format(pName=pName), lambda event : multi_getattr(event, lepton).triggerEffMC          , float, mcOnly=False, help='{PART} trigger efficiency as measured in MC'                                       .format(PART=lepton)),
             NTupleVariable('{pName}_recEffWeight'  .format(pName=pName), lambda event : multi_getattr(event, lepton).recEffWeight          , float, mcOnly=False, help='{PART} reconstruction ID + Iso weight'                                             .format(PART=lepton)),
           ]  

# muon
def fillMuon( pName, muon ):
    return fillLepton(pName, muon) + \
           [
             NTupleVariable('{pName}_looseId'.format(pName=pName), lambda event : multi_getattr(event, muon).looseId(), float, mcOnly=False, help='Muon loose ID as defined in Heppy Muon'),
             NTupleVariable('{pName}_tightId'.format(pName=pName), lambda event : multi_getattr(event, muon).tightId(), float, mcOnly=False, help='Muon tight ID as defined in Heppy Muon'),
           ] 
    # JAN FIXME - do we need the MVA iso and does it exist?
    # fill(tree, '{pName}_mvaIso'.format(pName=pName), muon.mvaIso() )

# electron
def fillEle( pName, ele ):
    return fillLepton(pName, ele) + \
           [
             NTupleVariable('{pName}_mvaTrigV0'          .format(pName=pName), lambda event : multi_getattr(event, ele).mvaTrigV0()         , float, mcOnly=False, help='Electron mvaTrigV0 ID as defined in Heppy Muon'         ),
             NTupleVariable('{pName}_mvaNonTrigV0'       .format(pName=pName), lambda event : multi_getattr(event, ele).mvaNonTrigV0()      , float, mcOnly=False, help='Electron mvaNonTrigV0 ID as defined in Heppy Muon'      ),
             NTupleVariable('{pName}_looseId'            .format(pName=pName), lambda event : multi_getattr(event, ele).looseIdForEleTau()  , float, mcOnly=False, help='Electron looseIdForEleTau ID as defined in Heppy Muon'  ),
             NTupleVariable('{pName}_tightId'            .format(pName=pName), lambda event : multi_getattr(event, ele).tightIdForEleTau()  , float, mcOnly=False, help='Electron tightIdForEleTau ID as defined in Heppy Muon'  ),
             NTupleVariable('{pName}_numberOfMissingHits'.format(pName=pName), lambda event : multi_getattr(event, ele).lostInner()         , float, mcOnly=False, help='Electron lostInner ID as defined in Heppy Muon'         ),
             NTupleVariable('{pName}_passConversionVeto' .format(pName=pName), lambda event : multi_getattr(event, ele).passConversionVeto(), float, mcOnly=False, help='Electron passConversionVeto ID as defined in Heppy Muon'),
           ] 
    # JAN FIXME - do we need the MVA iso and does it exist?
    # fill(tree, '{pName}_mvaIso'.format(pName=pName), ele.mvaIso() )
    # fill(tree, '{pName}_mvaTrigV0'.format(pName=pName), ele.electronID('mvaTrigV0') )
    # fill(tree, '{pName}_mvaNonTrigV0'.format(pName=pName), ele.electronID('mvaNonTrigV0') )

# tau 
def fillTau( pName, tau ):

    ids = []
    for id in tauIDs :
      ids.append( NTupleVariable('{pName}_{tauID}'.format(pName=pName, tauID=id), lambda event : multi_getattr(event, tau).tauID(id), float, mcOnly=False, help='Tau {tauID} discriminator'.format(tauID=id)) )
    return fillLepton(pName, tau) + ids + \
           [
             NTupleVariable('{pName}_decayMode'.format(pName=pName), lambda event : multi_getattr(event, tau).decayMode(), float, mcOnly=False, help='Tau decay mode'                                      ),
             NTupleVariable('{pName}_mass'     .format(pName=pName), lambda event : multi_getattr(event, tau).mass()     , float, mcOnly=False, help='Tau mass (RIC: isn\'t it the same as the p4.mass()?)'),
             NTupleVariable('{pName}_zImpact'  .format(pName=pName), lambda event : multi_getattr(event, tau).zImpact()  , float, mcOnly=False, help='Tau z impact parameter'                              ),
           ]

# jet
def fillJet( pName, jet ):
    # JAN - only one PU mva working point, but we may want to add more
    # run in our skimming step
    # (for which Jet.py would have to be touched again)
    return fillParticle(pName, jet) + \
           [
             NTupleVariable('{pName}_puMva'        .format(pName=pName), lambda jet : jet.puMva('pileupJetIdFull:full53xDiscriminant'), float, mcOnly=False, help='Tau decay mode '),
             NTupleVariable('{pName}_looseJetId'   .format(pName=pName), lambda jet : jet.looseJetId()                                , float, mcOnly=False, help='Tau decay mode '),
             NTupleVariable('{pName}_btagMVA'      .format(pName=pName), lambda jet : jet.btagMVA                                     , float, mcOnly=False, help='Tau decay mode '),
             NTupleVariable('{pName}_area'         .format(pName=pName), lambda jet : jet.jetArea()                                   , float, mcOnly=False, help='Tau decay mode '),
             NTupleVariable('{pName}_genJetPt'     .format(pName=pName), lambda jet : jet.matchedGenJet.pt()                          , float, mcOnly=True , help='Tau decay mode '),
             NTupleVariable('{pName}_partonFlavour'.format(pName=pName), lambda jet : jet.partonFlavour()                             , float, mcOnly=False, help='Tau decay mode '),
           ]     