# Calculates nucleotide frequencies and GC content in user-specified FASTA DNA file

import string

# Gets a file object and returns a dictionary of the frequencies of each nucleotide and of gc content
def bases_freq(dna_seq_file):
    title = dna_seq_file.readline()
    # Initialise dictionary
    nuclFreq = {'a':0, 'c':0, 'g':0, 't':0, 'gc':0, 'n':0}
    # Read line by line, because of large files
    for line in dna_seq_file:
        #If the sequence is in capital letters, convert it to lower letters (saves doing some actions twice)
        if line.isupper():
            line = line.lower()
        nuclFreq['a'] += line.count("a")
        nuclFreq['c'] += line.count("c")
        nuclFreq['g'] += line.count("g")
        nuclFreq['t'] += line.count("t")
        nuclFreq['gc'] += line.count("c") + line.count("g")
        nuclFreq['n'] += line.count("n") # To reduce Ns from final count
    # Total length of sequence
    seqLen = nuclFreq.get('a') + nuclFreq.get('c') + nuclFreq.get('g') + nuclFreq.get('t') - nuclFreq.get('n')
    dna_seq_file.close()
    # Convert all values to percentages
    for val in nuclFreq:
        nuclFreq[val] = float(nuclFreq[val]) / ((seqLen) / 100)
    return nuclFreq

fileName = input("Enter the file name:\n")
# Ensures file name contains suffix
if "." not in fileName:
    fileName += ".txt"
DNAfile = open(fileName, "r")
freq = bases_freq(DNAfile)
print ("a\tc\tg\tt\tcg")
print ("%.1f    %.1f    %.1f    %.1f    %.1f" %(freq['a'], freq['c'], freq['g'], freq['t'], freq['gc']))

'''
Example run:

Enter the file name:
sample_dna_1.txt
a	c	g	t	cg
27.4    22.7    21.9    28.1    44.5
'''
