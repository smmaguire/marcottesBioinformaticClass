#microtus avpr1 simulation

# this is to parse the fasta file. Probably can be greately improved but w/e it works. Only need to use it once to read in reference genome
import random
import string
def parse_fasta(fasta):
    fullSeq=[]
    seq = list(line.rstrip("\n").split(None,1)[0] for line in open(fasta))
    for line in seq:
        if line[0] != '>':
            fullSeq.append(line)
    fullSeq=string.join(fullSeq,"")
    return fullSeq

#here we are loading in the same txt seq for both alleals we'll need to make something that
#changes one or both of the alleals to have the SNPs that we're looking for. Should be a function of some kind.
alleal_1=parse_fasta("avpr1_cds.txt")
alleal_2=parse_fasta("avpr1_cds.txt")

# function that returns a random read from a randomly selected alleal (skipping reverse compliment stuff... Does this matter? should we add?)
def make_read(alleal1,alleal2):
    READLEN=50 #hard coded. Probably we should sample read lengths from a normal distribution centered at the mean of the real read lengths. 
    which=random.randint(1,2)
    if(which == 1):
        alleal=alleal1
    else:
        alleal=alleal2
    # pick a random starting position between 0 and the length of the gene, minus the read length
    pos = random.randint(0, len(alleal) - READLEN) 
    # extract the read sequence
    read = alleal[pos:pos + READLEN]
    return read

# error rate hardcoded to be 1 in 50.
# note -- about 1 error every 2 reads (should be higher/ lower?)
# the error rate is set by the random.randint function call, below -- random.randint returns a number between
# the first and second parameter, inclusive.  For a 10% error rate you could use 1, 10.

# should probably make a function that simulates methylation at some constant level (say 50%). This could run
# inside of the error function or seperately.

def make_substitution_error(read):
    for i in range(len(read)):
        if random.randint(1, 50) == 1:
            new_nt = random.choice("acgt") # these are lower case to make the identifiable, don't know if that is a problem for the sequence aligner
            read = read[:i] + new_nt + read[i + 1:]
    return read

#test
r = make_read(alleal_1,alleal_2)
r2 = make_substitution_error(r)
print r2

def make_sim_seq_set(animalID,numRead,alleal1,alleal2):
    outfp2 = open(animalID, 'w')
    for i in range(numRead):
        read = make_read(alleal1,alleal2)
        read = make_substitution_error(read)
        outfp2.write('>read%d\n%s\n' % (i, read))
    outfp2.close()

# how many reads per animal?
make_sim_seq_set("test.fa",1000,alleal_1,alleal_2)

