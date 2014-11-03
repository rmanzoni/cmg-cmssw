#for a in data WZ ZZ WW tt0l tt1l tt2l W1jet W2jet W3jet W4jet DY DY1 DY2 DY3 DY4
for a in data WZ ZZ TTW TTZ
#for a in TTW TTZ
#for a in WW W1jet W2jet W3jet W4jet DY DY1 DY2 DY3 DY4
#for a in data WZ ZZ tt1l tt2l
  do
  python wjet_control.py --channel muon --phys $a&
  python wjet_control.py --channel electron --phys $a&
done

