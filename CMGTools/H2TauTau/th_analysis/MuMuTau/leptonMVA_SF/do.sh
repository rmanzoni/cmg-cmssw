## for sample generation
#for a in Wjet tW tbW WZ ZZ data WW DY tt1l tt2l
for a in Wjet data DY tt1l tt2l
  do
  python sync.py --phys $a &
done

