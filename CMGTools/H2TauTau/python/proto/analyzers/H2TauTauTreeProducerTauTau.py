from CMGTools.H2TauTau.proto.analyzers.ntuple import *

htt_tt_globalVariables = diLeptonVars + ZttGenVars + fillTau('l1','diLepton.leg1') + fillTau('l2','diLepton.leg1')
htt_tt_globalObjects   = {}
htt_tt_collections     = {}

# 
# 
# class H2TauTauTreeProducerTauTau( TreeAnalyzerNumpy ):
#   '''Tree producer for the H->tau tau analysis'''
#   
#   def process(self, event):
#      
#     fill(self.tree, 'run'  , event.run    )
#     fill(self.tree, 'lumi' , event.lumi   )
#     fill(self.tree, 'event', event.eventId)
#         
#     fillDiLepton(self.tree, event.diLepton)  
#     
# #     import pdb ; pdb.set_trace()
# 
# #     fill(self.tree, 'metcov00', event.diLepton.metSig().significance()(0,0))
# #     fill(self.tree, 'metcov01', event.diLepton.metSig().significance()(0,1))
# #     fill(self.tree, 'metcov10', event.diLepton.metSig().significance()(1,0))
# #     fill(self.tree, 'metcov11', event.diLepton.metSig().significance()(1,1))
#     fill(self.tree, 'metPhi'  , event.diLepton.met().phi())
#     fill(self.tree, 'mex'     , event.diLepton.met().px() )
#     fill(self.tree, 'mey'     , event.diLepton.met().py() )
#     fill(self.tree, 'met'     , event.diLepton.met().pt() )
#     
#     fillTau(self.tree, 'l1', event.leg1 )
#     fillTau(self.tree, 'l2', event.leg2 )
# 
#     if hasattr(event,'genMass'):
#       fill(self.tree, 'genMass', event.genMass )
#     else:
#       fill(self.tree, 'genMass'         , -1 )
# 
#     nJets = len(event.cleanJets)
#     fill(self.tree, 'nJets', nJets )
#     if nJets >= 1 :
#       fillJet(self.tree, 'jet1'         , event.cleanJets[0]                                         )
#       fill   (self.tree, 'jet1Btag'     , event.cleanJets[0].btag('combinedSecondaryVertexBJetTags') )
# #       fill   (self.tree, 'jet1Bmatch'   , event.cleanJets[0].matchGenParton                          )
#     if nJets>=2:
#       fillJet(self.tree, 'jet2'         , event.cleanJets[1]                                         )
#       fill   (self.tree, 'jet2Btag'     , event.cleanJets[1].btag('combinedSecondaryVertexBJetTags') )
# #       fill   (self.tree, 'jet2Bmatch'   , event.cleanJets[1].matchGenParton                          )
#     
#     if hasattr(event, 'vbf') :
#       fillVBF(self.tree, 'ditau', event.vbf)
# 
#     nbJets = len(event.cleanBJets)
#     fill(self.tree, 'nbJets', nbJets )
#     if nbJets>=1:
#       fillJet(self.tree, 'bjet1', event.cleanBJets[0] )
#     if nbJets>=2:
#       fillJet(self.tree, 'bjet2', event.cleanBJets[1] )
# 
#     fill(self.tree, 'weight', event.eventWeight  )
#     fill(self.tree, 'nVert' , len(event.vertices)) 
# #     fill(self.tree, 'NUP'   , event.NUP          )
#     
#     fill(self.tree, 'isZtt'  , event.isZtt   )
#     fill(self.tree, 'isZee'  , event.isZee   )
#     fill(self.tree, 'isZmm'  , event.isZmm   )
#     fill(self.tree, 'isZj'   , event.isZj    )
#     fill(self.tree, 'isZttll', event.isZttll )
#     fill(self.tree, 'isZttj' , event.isZttj  )
#     
#     hasW = 0
#     if hasattr(event,'hasW') : hasW = event.hasW
#     fill(self.tree, 'hasW', hasW)
# 
#     hasZ = 0
#     if hasattr(event,'hasZ') : hasZ = event.hasZ
#     fill(self.tree, 'hasZ', hasZ)
# 
#     self.tree.tree.Fill() # funny
#     #self.tree.fill()
