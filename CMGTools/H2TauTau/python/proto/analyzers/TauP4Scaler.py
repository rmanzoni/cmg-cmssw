import os
import random
import ROOT

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer

random.seed(29061986)

class TauP4Scaler(Analyzer):

    '''Calibrates tau four momentum'''

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(TauP4Scaler, self).__init__(cfg_ana, cfg_comp, looperName)

        fname1 = '/'.join( [os.environ['CMSSW_BASE'], 
                            'src', 
                            'CMGTools', 
                            'H2TauTau', 
                            'data', 
                            'tauEnergyResponse.root'] )

        fname2 = '/'.join( [os.environ['CMSSW_BASE'], 
                            'src', 
                            'CMGTools', 
                            'H2TauTau', 
                            'data', 
                            'tauEnergyResponsePeak.root'] )

        self.f1 = ROOT.TFile.Open(fname1, 'read')
        self.f2 = ROOT.TFile.Open(fname2, 'read')
        
    def process(self, event):
        '''
        '''
        if   self.cfg_ana.method == 'stat': self._scaleByStat(event)
        elif self.cfg_ana.method == 'mean': self._scaleByMean(event)
        elif self.cfg_ana.method == 'peak': self._scaleByPeak(event)

    def getTauProperties(self, event):
        tau = getattr(event, self.cfg_ana.leg)
        
        dm      = tau.decayMode
        if dm > 1:
            dm = 10 # treat 2 prongs as 3 prongs
        barrel  = abs(tau.eta()) < 1.479
        endcap  = abs(tau.eta()) >= 1.479
        pt      = tau.pt()
        pttuple = self.findPtBin(tau.pt())
        
        return tau, dm, barrel, endcap, pt, pttuple       

    def scaleP4(self, tau, scale, event):
       
        modifiedP4 = ROOT.TLorentzVector()
        modifiedP4.SetPtEtaPhiM(
            tau.pt() * scale,
            tau.eta(),
            tau.phi(),
            tau.mass() # do not scale mass
        )
        
        # I love ROOT
        modifiedP4LV = ROOT.LorentzVector(
            modifiedP4.Px(),
            modifiedP4.Py(),
            modifiedP4.Pz(),
            modifiedP4.E(),
        )
        
        if self.cfg_ana.scaleMET:
            metP4 = event.diLepton.met().p4()
            metP4 += modifiedP4LV
            metP4 -= tau.p4()
            event.diLepton.met().setP4(metP4)

        if self.cfg_ana.leg == 'leg1':
            event.leg1.setP4(modifiedP4LV)
        elif self.cfg_ana.leg == 'leg2':
            event.leg2.setP4(modifiedP4LV)
       
    def findPtBin(self, pt):
        ptbins = [20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 100, 140, 200, 250, 400, 600]
        jj = -1
        for i, bin in enumerate(sorted(ptbins)):
            if pt > bin:
                jj = i
        if jj >= 0 and jj < len(ptbins)-1: 
            return str(ptbins[jj]), str(ptbins[jj+1])
        elif jj < 0:
            return None
        else:
            return '600', 'Inf'
    
    def _scaleByMean(self, event):
        '''
        Scale the tau p4 by the quantity corresponding to the mean of the 
        generated to reco transfer function.
        '''
        scale = 1.
        
        tau, dm, barrel, endcap, pt, pttuple = self.getTauProperties(event)
        
        # get the appropriate transfer fucntion
        table = self.f1.Get('h1_inv_transfer_baseline_%s_dm%d_pt%sto%s' %('barrel'*barrel + 'endcap'*endcap, dm, pttuple[0], pttuple[1]))
        
        scale = table.GetMean()
        
        self.scaleP4(tau, scale, event)
        

    def _scaleByPeak(self, event):
        '''
        Scale the tau p4 by the quantity corresponding to the peak of the 
        generated to reco transfer function.
        '''
        scale = 1.
        
        tau, dm, barrel, endcap, pt, pttuple = self.getTauProperties(event)
        
        # get the appropriate transfer fucntion
        table = self.f2.Get('%s_DM_%d' %('barrel'*barrel + 'endcap'*endcap, dm))
        
        bin = table.FindBin(pt)
        scale = table.GetBinContent(bin)

        self.scaleP4(tau, scale, event)


#     def _scaleByStat(self, event):
#         '''
#         Scale the tau p4 by a random value distributed according to the
#         generated to reco transfer function PDF.
#         '''
#         scale = 1.
# 
#         tau, dm, barrel, endcap, pt, pttuple = self.getTauProperties(event)
#         
#         # get the appropriate transfer fucntion
#         tableGoverR = self.f1.Get('h1_inv_transfer_baseline_%s_dm%d_pt%sto%s' %('barrel'*barrel + 'endcap'*endcap, dm, pttuple[0], pttuple[1]))
#         # normalise the histogram to get a PDF
#         tableGoverR.Scale(1./tableGoverR.Integral()) 
# 
#         tableRoverG = self.f1.Get('h1_transfer_baseline_%s_dm%d_pt%sto%s' %('barrel'*barrel + 'endcap'*endcap, dm, pttuple[0], pttuple[1]))
#         # normalise the histogram to get a PDF
#         tableRoverG.Scale(1./tableRoverG.Integral()) 
# 
#         probMax  = 0.
#         bestPtGen = 1.
#         
#         for ptgen in range(1, 600, 1):
#             
#             ptgen = float(ptgen)
#             
#             probGoverR = tableGoverR.FindBin(ptgen/pt)
#             probRoverG = tableRoverG.FindBin(pt/ptgen)
#             
#             if probGoverR * probRoverG > probMax:
#                 bestPtGen = ptgen
#         
#         scale = bestPtGen / pt    
#                 
#         self.scaleP4(tau, scale, event)

        
    def _scaleByStat(self, event):
        '''
        Scale the tau p4 by a random value distributed according to the
        generated to reco transfer function PDF.
        '''
        scale = 1.

        tau, dm, barrel, endcap, pt, pttuple = self.getTauProperties(event)
        
        # get the appropriate transfer fucntion
        table = self.f1.Get('h1_inv_transfer_baseline_%s_dm%d_pt%sto%s' %('barrel'*barrel + 'endcap'*endcap, dm, pttuple[0], pttuple[1]))
        # normalise the histogram to get a PDF
        table.Scale(1./table.Integral()) 
        
        toss = True
        i = 0
        while toss:
            i += 1
            if self.cfg_ana.verbose and i > 100 and i%100 == 0:
                print 'WARNING: tossing many (%d) attemps to recalibrate taus' %i
            scale = random.uniform(0., 2.)
            bin = table.FindBin(scale)
            prob = table.GetBinContent(bin)
            if prob == 0.:
                continue
            if random.uniform(0., 1.) <= prob:
                toss = False
        
        if self.cfg_ana.verbose:
            print 'INFO: pt %f, eta %f, decay modes %d \nscale %f after %d attempts' %(tau.pt(), tau.eta(), tau.decayMode(), scale, i)
        
        self.scaleP4(tau, scale, event)
