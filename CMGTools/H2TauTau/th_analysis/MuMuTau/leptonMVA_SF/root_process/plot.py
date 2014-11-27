import ROOT, os, numpy
#from officialStyle import officialStyle
from ROOT import TLegend, TCanvas
from CMGTools.RootTools.DataMC.DataMCPlot import DataMCPlot
from CMGTools.H2TauTau.proto.plotter.officialStyle import officialStyle
from ROOT import TColor, kMagenta, kOrange, kRed, kBlue, kGray, kBlack

ROOT.gROOT.SetBatch(True)
officialStyle(ROOT.gStyle)
#ROOT.gStyle.SetOptTitle(1)
#ROOT.gStyle.SetPadLeftMargin(0.18)
#ROOT.gStyle.SetPadBottomMargin(0.15)

col_qcd = TColor.GetColor(250,202,255)
col_tt  = TColor.GetColor(155,152,204)
col_ttv  = TColor.GetColor(155,182,204)
col_ewk = TColor.GetColor(222,90,106)
col_zll = TColor.GetColor(100,182,232)
col_red = TColor.GetColor(248,206,104)

colours = [2, 3, 4, 6, 7, 8]
nbin=40

variables = {
    'muon_pt': {'nbin':nbin, 'xtitle':'tag muon pT (GeV)', 'xmin':0, 'xmax':200},
    'muon_eta': {'nbin':nbin, 'xtitle':'tag muon pT (GeV)', 'xmin':-2.5, 'xmax':2.5},
    'muon_neu_iso': {'nbin':nbin, 'xtitle':'tag muon neutral iso', 'xmin':0, 'xmax':0.3},
    'muon_ch_iso': {'nbin':nbin, 'xtitle':'tag muon charged iso', 'xmin':0, 'xmax':0.15},
    'muon_dB3D': {'nbin':nbin, 'xtitle':'log(abs(tag muon MVA dB3D))', 'xmin':-15, 'xmax':0.},
    'muon_dz': {'nbin':nbin, 'xtitle':'log(abs(tag muon dz))', 'xmin':-15, 'xmax':0},
    'muon_dxy': {'nbin':nbin, 'xtitle':'log(abs(tag muon dxy))', 'xmin':-15, 'xmax':0},
    'muon_csv': {'nbin':nbin, 'xtitle':'tag muon CSV', 'xmin':0, 'xmax':1},
    'muon_jet_dr': {'nbin':nbin, 'xtitle':'tag muon jet dR', 'xmin':0, 'xmax':0.4},
    'muon_ptratio': {'nbin':nbin, 'xtitle':'tag muon p_{T} ratio', 'xmin':0, 'xmax':1.5},
    'muon_mva': {'nbin':nbin, 'xtitle':'tag muon MVA', 'xmin':-1, 'xmax':1.},

    'smuon_pt': {'nbin':nbin, 'xtitle':'probe muon pT (GeV)', 'xmin':0, 'xmax':200},
    'smuon_eta': {'nbin':nbin, 'xtitle':'probe muon pT (GeV)', 'xmin':-2.5, 'xmax':2.5},
    'smuon_neu_iso': {'nbin':nbin, 'xtitle':'probe muon neutral iso', 'xmin':0, 'xmax':0.3},
    'smuon_ch_iso': {'nbin':nbin, 'xtitle':'probe muon charged iso', 'xmin':0, 'xmax':0.15},
    'smuon_dB3D': {'nbin':nbin, 'xtitle':'log(abs(probe muon MVA dB3D))', 'xmin':-15, 'xmax':0.},
    'smuon_dz': {'nbin':nbin, 'xtitle':'log(abs(probe muon dz))', 'xmin':-15, 'xmax':0},
    'smuon_dxy': {'nbin':nbin, 'xtitle':'log(abs(probe muon dxy))', 'xmin':-15, 'xmax':0},
    'smuon_csv': {'nbin':nbin, 'xtitle':'probe muon CSV', 'xmin':0, 'xmax':1},
    'smuon_jet_dr': {'nbin':nbin, 'xtitle':'probe muon jet dR', 'xmin':0, 'xmax':0.4},
    'smuon_ptratio': {'nbin':nbin, 'xtitle':'probe muon p_{T} ratio', 'xmin':0, 'xmax':1.5},
    'smuon_mva': {'nbin':nbin, 'xtitle':'probe muon MVA', 'xmin':-1, 'xmax':1.},

    'evt_Mmm': {'nbin':65, 'xtitle':'Invariant Mass (ll)', 'xmin':20, 'xmax':150},
    'evt_met': {'nbin':nbin, 'xtitle':'missing E_{T}', 'xmin':0, 'xmax':200},
}

selections = {
    'ttbar':{'selection':'evt_isMC==1 && evt_id >= 3 && evt_id <=4', 'col':col_tt, 'order':2},
    'DY':{'selection':'evt_isMC==1 && evt_id==6', 'col':col_zll, 'order':1},
    'W':{'selection':'evt_isMC==1 && evt_id==11', 'col':col_qcd, 'order':4},
    'data':{'selection':'evt_isMC==0', 'col':1, 'order':2999}
    }

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def makePlotsVars(tree, isSignal=False):
   
    for ivar, var in variables.items():

        c = TCanvas()
        h_stack = DataMCPlot('hs_' + ivar)
        h_stack.legendBorders = 0.65, 0.6, 0.85, 0.85
        
        for iprocess, process in selections.items():

            if not isSignal and iprocess=='signal':
                continue

            hist = ROOT.TH1F('h_'+ivar+'_'+iprocess, '', var['nbin'], var['xmin'], var['xmax'])
            hist.GetXaxis().SetTitle(var['xtitle'])
            hist.GetYaxis().SetNdivisions(507)
            hist.Sumw2()

            sel = '(' + process['selection'] + ')*evt_weight'

            tree.Project(hist.GetName(), ivar, sel)


            if iprocess=='data':
                hist.SetMarkerStyle(20)
            elif iprocess=='signal':
                hist.Scale(10.)
                hist.SetLineWidth(2)
                hist.SetLineColor(process['col'])
                hist.SetMarkerColor(process['col'])
                hist.SetMarkerSize(0)
            else:
                hist.SetFillStyle(1)
                hist.SetLineWidth(2)
                hist.SetFillColor(process['col'])
                hist.SetLineColor(process['col'])
                
            
            h_stack.AddHistogram(hist.GetName(), hist, process['order'])
            h_stack.Hist(hist.GetName()).legendLine = iprocess

            if iprocess=='data' or iprocess=='signal':
                h_stack.Hist(hist.GetName()).stack = False
            if iprocess=='signal':
                h_stack.Hist(hist.GetName()).legendLine = 'signal (x10)'

        print h_stack
        h_stack.DrawStack('HIST')
        
        ensureDir('plots')
        cname =  'plots/' + ivar + '.gif'
        print cname
        c.Print(cname)
        c.Print(cname.replace('gif','pdf'))

if __name__ == '__main__':

    tfile = ROOT.TFile('Zmumu.root')
    tree = tfile.Get('Tree')
    makePlotsVars(tree)
