import os
import random
import ROOT

from PhysicsTools.Heppy.analyzers.core.Analyzer import Analyzer

random.seed(29061986)

class TauP4Scaler(Analyzer):

    '''Calibrates tau four momentum'''

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(TauP4Scaler, self).__init__(cfg_ana, cfg_comp, looperName)

        fname = '/'.join( [os.environ['CMSSW_BASE'], 
                           'src', 
                           'CMGTools', 
                           'H2TauTau', 
                           'data', 
                           'tauEnergyResponse.root'] )

        self.f1 = ROOT.TFile.Open(fname, 'read')
        
    def process(self, event):
        tau = getattr(event, self.cfg_ana.leg)
        
        dm     = tau.decayMode
        if dm > 1:
            dm = 10 # treat 2 prongs as 3 prongs
        barrel = abs(tau.eta()) < 1.479
        endcap = abs(tau.eta()) >= 1.479
        pt     = self.findPtBin(tau.pt())
        
        # get the appropriate transfer fucntion
        table = self.f1.Get('h1_inv_transfer_baseline_%s_dm%d_pt%sto%s' %('barrel'*barrel + 'endcap'*endcap, dm, pt[0], pt[1]))
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
    

