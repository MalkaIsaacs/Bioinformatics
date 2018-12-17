# This program gets a motive and a file name with sequences in FASTA format and prints a list of
# gi numbers, the length of the sequences and if the given motive exists in the sequences.
# The program will then state the total number of sequences.

import string

fileName = input("Enter the name of the file, with the correct suffix.\n")
# If file name does not contain a suffix, default '.txt' is assumed
if "." not in fileName:
    fileName = fileName + '.txt'
FASTAfile = open(fileName, 'r')
motive = input("Enter the motive to look for:\n")
# Ensure upper case
if motive.islower():
    motive = motive.upper()
print (">gi number     Length   Motive\n\
>-------------------------------------------\n")
numOfSeq = 0
for sequence in FASTAfile.read().split("\n\n"):
    if sequence.startswith(">gi"):
        # Count number of sequences
        numOfSeq += 1
        # Split by lines and then join all lines but the first, i.e. the actual sequence without the title
        seq = "".join(sequence.split("\n")[1 : ])
        # Just in case there's an extra space at the end...
        seq.rstrip()
        # Ensure uppercase
        if seq.islower():
            seq = seq.upper()
        # If specified motive exists in sequence
        if motive in seq:
            mot = "+\n"
        else:
            mot = "-\n"
    print (sequence[0 : 12], "\t", len(seq), "\t  ", mot)
FASTAfile.close()
print ("Total of %d sequences\n" % numOfSeq)
