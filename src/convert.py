

src = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI"]
dst = ["chr"+str(i) for i in range(1,17)]

map = dict(zip(src,dst))

with open("data/Eaton_2010_ORC_ACS_V64.bed") as f:
    for line in f:
        line = line.strip()
        if line.startswith("chr"):
            i = line.index("\t")
            curr = line[:i]
            print map[curr]+line[i:]
        else:
            print line
