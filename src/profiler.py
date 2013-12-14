from collections import Counter
from Bio import SeqIO
import sys

class Profiler:
    def __init__(self):
        self.acs = self.read_acs()
        self.ori = self.read_ori()
        sys.stderr.write( "reading genome...\n")
        self.genome = self.read_genome()
        sys.stderr.write( str(self.genome.keys()) )

        sys.stderr.write( "reading done\n")
        self.trim_cnt = 0
        self.profile_half_len = 10000
        self.half_num_samples = 300
        self.profile_lens = []

    def read_genome(self):
        result = dict()
        for chromo in ["chr"+str(i) for i in range(1,17)]:
            with open("data/sacCer1/"+chromo+".fa", "r") as f:
                fasta = SeqIO.parse(f, "fasta")
                contig = fasta.next()
                result[chromo] = str(contig.seq).upper()
        return(result)

    def read_acs(self):
        result = dict()
        with open("data/Eaton_2010_ORC_ACS_V64_converted.bed", "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("chr"):
                    cols = line.split("\t")
                    chromo = cols[0]
                    rec = (int(cols[1]), int(cols[2]), cols[3])
                    if chromo in result:
                        result[chromo].append(rec)
                    else:
                        result[chromo] = [rec]
        return result

    def read_ori(self):
        result = dict()
        with open("data/all_available_ori_scerevisiae.txt", "r") as f:
            cnt = 0
            for line in f:
                if cnt >0:
                    line = line.strip()
                    cols = line.split("\t")
                    if cols[5]=="Confirmed":
                        chromo = "chr"+cols[0]
                        rec = (int(cols[1]), int(cols[2]), cols[3])
                        if chromo in result:
                            result[chromo].append(rec)
                        else:
                            result[chromo] = [rec]
                cnt += 1
        return result

    def get_prev_ars(self, chromo, acs_center):
        lst = self.ori[chromo]
        prev = None
        for o in lst:
            if o[1] > acs_center:
                break
            prev = o
        return(prev)

    def get_prev_distance(self, acs_center, ars):
        if ars==None:
            result = acs_center/2
        else:
            result = (acs_center-ars[1])/2
        return(result)

    def get_next_ars(self, chromo, acs_center):
        lst = self.ori[chromo]
        lst = [i for i in reversed(lst)]
        prev = None
        for o in lst:
            if o[0] < acs_center:
                break
            prev = o
        return prev

    def get_next_distance(self, acs_center, ars, chrlen):
        if ars==None:
            result = (chrlen-acs_center)/2
        else:
            result = (ars[0]-acs_center)/2
        return(result)

    def getNuc(self, chromo, pos):
        return(self.genome[chromo][pos])

    def run(self):
        result = []
        for chromo in ["chr"+str(i) for i in range(1,17)]:
            sys.stderr.write("doing "+chromo+"\n")
            for acs_rec in self.acs[chromo]:
                acs_id = acs_rec[2]
                acs_center = acs_rec[0]+17

                cnt = Counter()
                prev_ars = self.get_prev_ars(chromo, acs_center)
                dist = self.get_prev_distance(acs_center, prev_ars)

                if dist>self.profile_half_len:
                    dist = self.profile_half_len
                elif dist<self.half_num_samples:
                    sys.stderr.write("skipping "+acs_id+"\n")
                    continue

                next_ars = self.get_next_ars(chromo, acs_center)
                tmpdist = self.get_next_distance(acs_center, next_ars, len(self.genome[chromo]))
                if tmpdist<self.half_num_samples:
                    sys.stderr.write("skipping "+acs_id+"\n")
                    continue

                profile_start = acs_center-dist

                step = 1.0*dist/self.half_num_samples
                curr = 1.0*profile_start
                sample_num = 0
                while(sample_num<self.half_num_samples):
                    cnt.update(self.getNuc(chromo, int(curr)))
                    result.append((acs_id, int(curr), cnt["A"], cnt["T"], cnt["C"], cnt["G"]))
                    sample_num += 1
                    curr += step


                dist = self.get_next_distance(acs_center, next_ars, len(self.genome[chromo]))
                if dist > self.profile_half_len:
                    dist = self.profile_half_len
                elif dist<self.half_num_samples:
                    sys.stderr.write("skipping "+acs_id+"\n")
                    continue

                step = 1.0*dist/self.half_num_samples
                curr = 1.0*acs_center
                while(sample_num<2*self.half_num_samples):
                    cnt.update(self.getNuc(chromo, int(curr)))
                    result.append((acs_id, int(curr), cnt["A"], cnt["T"], cnt["C"], cnt["G"]))
                    sample_num += 1
                    curr += step

        return result

    def debug_neighbours(self):
        print("acs_id, chromo, pos, prev_ars, pstart, pend, next_ars, nstart, nend")
        for chromo in ["chr"+str(i) for i in range(1,17)]:
            for acs_rec in self.acs[chromo]:
                acs_id = acs_rec[2]
                acs_center = acs_rec[0]+17

                prev_ars = self.get_prev_ars(chromo, acs_center)
                prev_id = "None"
                prev_start = -1
                prev_end=-1
                if prev_ars!=None:
                    prev_id = prev_ars[2]
                    prev_start = prev_ars[0]
                    prev_end = prev_ars[1]
                prev_distance = self.get_prev_distance(acs_center, prev_ars)

                next_ars = self.get_next_ars(chromo, acs_center)
                next_id = "None"
                next_start = -1
                next_end = -1
                if next_ars != None:
                    next_id = next_ars[2]
                    next_start = next_ars[0]
                    next_end = next_ars[1]
                next_distance = self.get_next_distance(acs_center, next_ars, len(self.genome[chromo]))
                print("%s, %s, %d, %s, %d, %d, %d, %s, %d, %d, %d"%(acs_id, chromo, acs_center,
                                                    prev_id, prev_start, prev_end, prev_distance,
                                                    next_id, next_start, next_end, next_distance,
                                                    ))

def main():
    prof()
    # dbg()

# results/neighbours.txt
def dbg():
    profiler = Profiler()
    profiler.debug_neighbours()

# results/acs_new.txt
def prof():
    profiler = Profiler()
    result = profiler.run()
    print "id,pos,cntA,cntT,cntC,cntG"
    for rec in result:
        print "%s,%d,%d,%d,%d,%d"%tuple(rec)

if __name__ == "__main__":
    main()
