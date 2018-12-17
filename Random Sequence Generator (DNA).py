#Generates random nucleotide sequence (useful for testing code, comparing frequencies, etc)

import random

dna = ['A', 'C', 'G', 'T']
size = int(input("How many nucleotides would you like the random sequence to contain?\n"))
i = 0
seq = ""
while i < size:
    seq += random.choice(dna)
    i += 1
fileName = "randomDNAsequence" + str(size) + ".txt"
outputFile = open(fileName, "w")
outputFile.write("Randomly-generated sequence:\n")
outputFile.write(seq)
outputFile.close()
print ("Successfully saved generated sequence to \'" + fileName + "\'\n")

'''
Example run:

How many nucleotides would you like the random sequence to contain?
20000
Successfully saved generated sequence to 'randomDNAsequence20000.txt'
'''
