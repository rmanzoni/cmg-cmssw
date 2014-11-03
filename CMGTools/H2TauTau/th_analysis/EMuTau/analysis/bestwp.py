from ROOT import TFile, TH3F, gStyle, TH1F, TCanvas, gROOT, TGraph
from officialStyle import officialStyle
import string

officialStyle(gStyle)
#gStyle.SetOptTitle(0)

def checkVal(file, str, n_ie, n_im):
    scan = TH1F("scan","scan",100,-2,2)
    cond = 'n_ie ==' + n_ie + ' && n_im ==' + n_im
    file.Tree.Draw(str + ' >> scan', cond)
    return str, ('%.3f' % scan.GetMean())
    
    
nsize = 51

file = TFile("Scan.root")
h = TH3F("h","h",nsize, 0, nsize, nsize, 0, nsize, 1000,0,1.)
g = TH3F("g","g",nsize, 0, nsize, nsize, 0, nsize, 1000,0,60.)
k = TH3F("k","k",nsize, 0, nsize, nsize, 0, nsize, 100,0,0.15)

file.Tree.Draw("nsig:n_im:n_ie >> h");
h_proj = h.Project3DProfile()

#file.Tree.Draw("significance:n_im:n_ie >> k");
file.Tree.Draw("nsig/TMath::Sqrt(nsig+nbkg+0.1*nbkg*nbkg):n_im:n_ie >> k");
k_proj = k.Project3DProfile()

file.Tree.Draw("nbkg:n_im:n_ie >> g");
g_proj = g.Project3DProfile()

xbin = []
ybin = []
nbkg = []
nsig = []

for ix in range(1, h_proj.GetXaxis().GetNbins()+1):
    
    max = 1000.
    max_Nsig = -1.

    _xbin_ = -1.
    _ybin_ = -1.
    _nbkg_ = -1.
    _nsig_ = -1.

    for iy in range(1, h_proj.GetYaxis().GetNbins()+1):

        Nsig = h_proj.GetBinContent(ix, iy)
        if Nsig < 0.44:
            continue
        
        if abs(Nsig-0.45) < max:
            max = abs(Nsig - 0.45)
            max_Nsig = Nsig
            _xbin_ = ix
            _ybin_ = iy
            _nsig_ = h_proj.GetBinContent(ix, iy)
            _nbkg_ = g_proj.GetBinContent(ix, iy)
            
    if max_Nsig != -1:
        xbin.append(_xbin_)
        ybin.append(_ybin_)
        nsig.append(_nsig_)
        nbkg.append(_nbkg_)

min = 1000.
min_x = -1.
min_y = -1.

for ii in range(len(xbin)):
    print xbin[ii], ybin[ii], nsig[ii], nbkg[ii]
    if nbkg[ii] < min:
        min = nbkg[ii]
        min_x = xbin[ii]
        min_y = ybin[ii]


print 'Results : Min. bkg = ', min, 'at (x, y) = ', min_x, min_y

graph = TGraph()

for a in range(len(xbin)):
    graph.SetPoint(a, xbin[a], ybin[a])

graph.SetMarkerSize(0)

c = TCanvas()
h_proj.SetTitle('Number of signal')
h_proj.GetYaxis().SetTitle('electron MVA step')
h_proj.GetXaxis().SetTitle('muon MVA step')
h_proj.Draw("colz")
graph.Draw("lsame")

c2 = TCanvas()
g_proj.SetTitle('Number of backgrounds')
g_proj.Draw("colz")
g_proj.GetYaxis().SetTitle('electron MVA step')
g_proj.GetXaxis().SetTitle('muon MVA step')
graph.Draw("lsame")

c3 = TCanvas()
k_proj.SetTitle('Number of signal')
k_proj.GetYaxis().SetTitle('electron MVA step')
k_proj.GetXaxis().SetTitle('muon MVA step')
k_proj.Draw("colz")

c3.SaveAs("significance.gif")
c3.SaveAs("significance.pdf")

e_th = str(min_x-1)
m_th = str(min_y-1)

print 'muon, barrel', checkVal(file, 'muon_b_mva', e_th, m_th)
print 'muon, endcap', checkVal(file, 'muon_e_mva', e_th, m_th)
print 'electron, barrel', checkVal(file, 'electron_b_mva', e_th, m_th)
print 'electron, endcap', checkVal(file, 'electron_e_mva', e_th, m_th)


