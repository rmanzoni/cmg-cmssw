from itertools import product

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer
from PhysicsTools.Heppy.analyzers.core.AutoHandle import AutoHandle
from PhysicsTools.HeppyCore.utils.deltar import deltaR

import PhysicsTools.HeppyCore.framework.config as cfg


class L1TriggerAnalyzer(Analyzer):

    def __init__(self, *args, **kwargs):
        super(L1TriggerAnalyzer, self).__init__(*args, **kwargs)
        
        self.l1objDict = {
            'EmIsolated'   : 0,
            'EmNonIsolated': 1,
            'MET'          : 2,
            'MHT'          : 3,
            'HFRings'      : 4,
            'CentralJet'   : 5,
            'ForwardJet'   : 6,
            'IsoTau'       : 7,
            'Tau'          : 8,
            'Muon'         : 9
        }
    
    def declareHandles(self):
        super(L1TriggerAnalyzer, self).declareHandles()
        
        if (self.cfg_ana, 'label'):
            label = self.cfg_ana.label
        else:
            label = 'l1extraParticles'
            
        self.handles['EmIsolated'   ] = AutoHandle( (label, 'Isolated'   ), 'std::vector<l1extra::L1EmParticle>'    )
        self.handles['EmNonIsolated'] = AutoHandle( (label, 'NonIsolated'), 'std::vector<l1extra::L1EmParticle>'    )
        self.handles['MET'          ] = AutoHandle( (label, 'MET'        ), 'std::vector<l1extra::L1EtMissParticle>')
        self.handles['MHT'          ] = AutoHandle( (label, 'MHT'        ), 'std::vector<l1extra::L1EtMissParticle>')
        self.handles['HFRings'      ] = AutoHandle( (label, ''           ), 'std::vector<l1extra::L1HFRings>'       )
        self.handles['CentralJet'   ] = AutoHandle( (label, 'Central'    ), 'std::vector<l1extra::L1JetParticle>'   )
        self.handles['ForwardJet'   ] = AutoHandle( (label, 'Forward'    ), 'std::vector<l1extra::L1JetParticle>'   )
        self.handles['IsoTau'       ] = AutoHandle( (label, 'IsoTau'     ), 'std::vector<l1extra::L1JetParticle>'   )
        self.handles['Tau'          ] = AutoHandle( (label, 'Tau'        ), 'std::vector<l1extra::L1JetParticle>'   )
        self.handles['Muon'         ] = AutoHandle( (label, ''           ), 'std::vector<l1extra::L1MuonParticle>'  )
        
    def process(self, event):
        self.readCollections(event.input)
        
        # specify in the cfg what collections should be used to match.
        # Options are:
        #   0  EmIsolated 
        #   1  EmNonIsolated
        #   2  MET         
        #   3  MHT       
        #   4  HFRings    
        #   5  CentralJet   
        #   6  ForwardJet   
        #   7  IsoTau       
        #   8  Tau          
        #   9  Muon         

        collections = self.cfg_ana.collections
        
        dRmax = 0.5
        if hasattr(self.cfg_ana, 'dR'):
            dRmax = self.cfg_ana.dR
            
        legs = {event.diLepton.leg1():dRmax,
                event.diLepton.leg2():dRmax}    
                    
        for coll in collections:
     
            mycoll = self.handles[coll].product()
            for leg, l1 in product(legs.keys(), mycoll):
                dR = deltaR(l1.eta(), l1.phi(), leg.eta(), leg.phi())
                if dR < legs[leg]:
                    leg.L1 = l1
                    leg.L1flavour = self.l1objDict[coll]
                    legs[leg] = dR  

        return True

setattr(L1TriggerAnalyzer, 'defaultConfig', 
    cfg.Analyzer(
        class_object=L1TriggerAnalyzer,
        collections=['IsoTau', 'Tau', 'Muon'],
        dR=0.5
    )
)
