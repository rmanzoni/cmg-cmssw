from CMGTools.H2TauTau.proto.analyzers.ntuple import *

htt_tt_globalVariables = htt_globalVariables           + \
                         fillTau('l1','diLepton.leg1') + \
                         fillTau('l2','diLepton.leg2') + \
                         fillGenParticle('l1_gen','diLepton.leg1.genParticle()', isTau = True, genVisTau = 'diLepton.leg1.physObj.genJet()') + \
                         fillGenParticle('l2_gen','diLepton.leg2.genParticle()', isTau = True, genVisTau = 'diLepton.leg2.physObj.genJet()') 
                         
htt_tt_globalObjects = htt_globalObjects.copy()
htt_tt_globalObjects.update( {} )

htt_tt_collections = htt_collections.copy()
htt_tt_collections.update( {} )