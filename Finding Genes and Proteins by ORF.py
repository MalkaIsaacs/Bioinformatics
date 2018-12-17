# Scans the 6 possible reading frames to find an ORF (Open Reading Frame) resembling a gene, of greater length
# than a user-speficied size. Then translates the gene into protein.

import string


# Gets a gene and translates it into a protein
def translateDNAtoProtein(gene):
    code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"s", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
    "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    i = 0
    protein = ""
    while(i < len(gene) - 1):
        # Get the next 3 nucleotides
        aminoAcid = gene[i] + gene[i+1] + gene[i+2]
        # Get the amino acid coded by the 3 nucleotides
        for aa in code.keys():
          if aa == aminoAcid:
              # Add amino acid to the protein, unless it's a STOP codon
              if code[aa] != "STOP":
                  protein += code[aa]
        i += 3
    print ("Protein sequence:\n" + protein + "\n\n")


# Looks for ORFs for genes (of minSize nucleotides or more) in DNA sequence it receives.
# If finds genes, prints them, sends them to be translated to protein and prints that too
def findGene(dnaSeq):
    # Set to beginning of sequence
    i = 0
    # Set to end of sequence
    j = len(dnaSeq) - 1
    # Counts nucleotides from start codon (for later use to check if gene is long enough)
    count = 0
    index = 0
    orfExists = 0
    # For the entire sequence
    while((j - i) > 1):
        # Start with the first 3 nucleotides for this search
        start = dnaSeq[i] + dnaSeq[i + 1] + dnaSeq[i + 2]
        index = i + 3
        flag = 1
        # If the first 3 nucleotides represent a start codon
        if start == "ATG":
            # Look at the next 3 nucleotides
            nextNuc = dnaSeq[index] + dnaSeq[index + 1] + dnaSeq[index + 2]
            count = 3
            # So that we know where to look for stop codon
            end = index
            gene = start + nextNuc
            stop = ""
            # Look for a stop codon (until one is found or the end of the sequence is reached)
            while(((j - end) > 1) and (flag > 0)):
                stop = dnaSeq[end] + dnaSeq[end + 1] + dnaSeq[end + 2]
                # Copies the 3 nucleotides into the gene sequence
                gene += stop
                # Sets to the next 3 nucleotides
                end += 3
                # Add 3 to the nucleotide count
                count += 3
                #If an appropriate stop codon is found
                if stop == "TAA" or stop == "TAG" or stop == "TGA":
                    # If there are at least minSize nucleotides in the ORF, print the details and call the function that translates from DNA to protein to do so
                    if count > minSize - 1:
                        orfExists = 1
                        secondLastCodon = dnaSeq[end - 6] + dnaSeq[end - 5] + dnaSeq[end - 4]
                        print ("start: %d\tend: %d\tlength: %d\t" % ((index - 1), (end + 1), ((end) - (index - 2))) + start + nextNuc + "\t" + secondLastCodon + stop + "\n")
                        # Print gene sequence
                        print ("Gene sequence:\n" + gene + "\n")
                        translateDNAtoProtein(gene)
                    flag = -1
                i = end
        # If start codon has not yet been found, keep checking nucleotides in sets of 3 (start the search process from 3 nucleotides later)
        else:
            i += 3     
    return orfExists


# Program:
fileName = input("Enter file name:\n")
# If no suffix entered, add it
if "." not in fileName:
    fileName += ".txt"
originalFile = open(fileName, "r")
# Extract the actual DNA sequence from the read file and capitalise it
dnaSeq = "".join(originalFile.read().split("\n")[1:]).strip().upper()
originalFile.close()

minSize = int(input("Enter minimum size for ORF:\n"))

print ("\nFrame 1 (+1):\n")
# Send the first frame to the function that looks for a gene in it. If no appropriate ORF was found, say so...
if findGene(dnaSeq) == 0:
    print ("No ORF of at least " + str(minSize) + " nucleotides has been found.\n\n")
# Save the sequence for later use (when reversing the complement)
dna = dnaSeq

# Start from the second nucleotide
dnaSeq = dnaSeq[1:]
print ("Frame 2 (+2):\n")
# Send the second frame to the function that looks for a gene in it. If no appropriate ORF was found, say so...
if findGene(dnaSeq) == 0:
    print ("No ORF of at least " + str(minSize) + " nucleotides has been found.\n\n")

# Start from the third nucleotide
dnaSeq = dnaSeq[1:]
print ("Frame 3 (+3):\n")
#Send the third frame to the function that looks for a gene in it. If no appropriate ORF was found, say so...
if findGene(dnaSeq) == 0:
    print ("No ORF of at least " + str(minSize) + " nucleotides has been found.\n\n")

complementDNA = ""
i = 0
# Convert the entire sequence to its complement
while(i < (len(dna))):
    if dna[i] == "G":
       complementDNA += "C"
    if dna[i] == "C":
       complementDNA += "G"
    if dna[i] == "A":
       complementDNA += "T"
    if dna[i] == "T":
       complementDNA += "A"
    i += 1
 
reverseDNA = ""
i = 0
j = len(complementDNA)
# Reverse the entire complement sequence
while(i < j):
    reverseDNA += complementDNA[j - 1]
    j -= 1

print ("Frame 4 (-1):\n")
# Send the fourth frame to the function that looks for a gene in it. If no appropriate ORF was found, say so...
if findGene(reverseDNA) == 0:
    print ("No ORF of at least " + str(minSize) + " nucleotides has been found.\n\n")

# Start from the second reverse nucleotide
reverseDNA = reverseDNA[1:]
print ("Frame 5 (-2):\n")
# Send the fifth frame to the function that looks for a gene in it. If no appropriate ORF was found, say so...
if findGene(reverseDNA) == 0:
    print ("No ORF of at least " + str(minSize) + " nucleotides has been found.\n\n")

# Start from the third reverse nucleotide
reverseDNA = reverseDNA[1:]
print ("Frame 6 (-3):\n")
# Send the sixth frame to the function that looks for a gene in it. If no appropriate ORF was found, say so...
if findGene(reverseDNA) == 0:
    print ("No ORF of at least " + str(minSize) + " nucleotides has been found.\n\n")


'''
Example run:

Enter file name:
plasmid pACYC184
Enter minimum size for ORF:
300

Frame 1 (+1):

No ORF of at least 300 nucleotides has been found.


Frame 2 (+2):

start: 1754	end: 2137	length: 383	ATGCGC	GCATAA

Gene sequence:
ATGCGCCGCACCCGTTCTCGGAGCACTGTCCGACCGCTTTGGCCGCCGCCCAGTCCTGCTCGCTTCGCTACTTGGAGCCACTATCGACTACGCGATCATGGCGACCACACCCGTCCTGTGGATCCTCTACGCCGGACGCATCGTGGCCGGCATCACCGGCGCCACAGGTGCGGTTGCTGGCGCCTATATCGCCGACATCACCGATGGGGAAGATCGGGCTCGCCACTTCGGGCTCATGAGCGCTTGTTTCGGCGTGGGTATGGTGGCAGGCCCCGTGGCCGGGGGACTGTTGGGCGCCATCTCCTTGCATGCACCATTCCTTGCGGCGGCGGTGCTCAACGGCCTCAACCTACTACTGGGCTGCTTCCTAATGCAGGAGTCGCATAA

Protein sequence:
MRRTRSRSTVRPLWPPPSPARFATWSHYRLRDHGDHTRPVDPLRRTHRGRHHRRHRCGCWRLYRRHHRWGRSGSPLRAHERLFRRGYGGRPRGRGTVGRHLLACTIPCGGGAQRPQPTTGLLPNAGVA


start: 2675	end: 3013	length: 338	ATGGAA	GGTTAG

Gene sequence:
ATGGAAGAACGGGTTGGCATGGATTGTAGGCGCCGCCCTATACCTTGTCTGCCTCCCCGCGTTGCGTCGCGGTGCATGGAGCCGGGCCACCTCGACCTGAATGGAAGCCGGCGGCACCTCGCTAACGGATTCACCACTCCAAGAATTGGAGCCAATCAATTCTTGCGGAGAACTGTGAATGCGCAAACCAACCCTTGGCAGAACATATCCATCGCGTCCGCCATCTCCAGCAGCCGCACGCGGCGCATCTCGGGCAGCGTTGGGTCCTGGCCACGGGTGCGCATGATCGTGCTCCTGTCGTTGAGGACCCGGCTAGGCTGGCGGGGTTGCCTTACTGGTTAG

Protein sequence:
MEERVGMDCRRRPIPCLPPRVASRCMEPGHLDLNGSRRHLANGFTTPRIGANQFLRRTVNAQTNPWQNIsIAsAIsSSRTRRISGSVGsWPRVRMIVLLSLRTRLGWRGCLTG


Frame 3 (+3):

start: 1580	end: 2770	length: 1190	ATGAAA	ACCTGA

Gene sequence:
ATGAAAAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTCACCCTGGATGCTGTAGGCATAGGCTTGGTTATGCCGGTACTGCCGGGCCTCTTGCGGGATATCGTCCATTCCGACAGCATCGCCAGTCACTATGGCGTGCTGCTAGCGCTATATGCGTTGATGCAATTTCTATGCGCACCCGTTCTCGGAGCACTGTCCGACCGCTTTGGCCGCCGCCCAGTCCTGCTCGCTTCGCTACTTGGAGCCACTATCGACTACGCGATCATGGCGACCACACCCGTCCTGTGGATCCTCTACGCCGGACGCATCGTGGCCGGCATCACCGGCGCCACAGGTGCGGTTGCTGGCGCCTATATCGCCGACATCACCGATGGGGAAGATCGGGCTCGCCACTTCGGGCTCATGAGCGCTTGTTTCGGCGTGGGTATGGTGGCAGGCCCCGTGGCCGGGGGACTGTTGGGCGCCATCTCCTTGCATGCACCATTCCTTGCGGCGGCGGTGCTCAACGGCCTCAACCTACTACTGGGCTGCTTCCTAATGCAGGAGTCGCATAAGGGAGAGCGTCGACCGATGCCCTTGAGAGCCTTCAACCCAGTCAGCTCCTTCCGGTGGGCGCGGGGCATGACTATCGTCGCCGCACTTATGACTGTCTTCTTTATCATGCAACTCGTAGGACAGGTGCCGGCAGCGCTCTGGGTCATTTTCGGCGAGGACCGCTTTCGCTGGAGCGCGACGATGATCGGCCTGTCGCTTGCGGTATTCGGAATCTTGCACGCCCTCGCTCAAGCCTTCGTCACTGGTCCCGCCACCAAACGTTTCGGCGAGAAGCAGGCCATTATCGCCGGCATGGCGGCCGACGCGCTGGGCTACGTCTTGCTGGCGTTCGCGACGCGAGGCTGGATGGCCTTCCCCATTATGATTCTTCTCGCTTCCGGCGGCATCGGGATGCCCGCGTTGCAGGCCATGCTGTCCAGGCAGGTAGATGACGACCATCAGGGACAGCTTCAAGGATCGCTCGCGGCTCTTACCAGCCTAACTTCGATCACTGGACCGCTGATCGTCACGGCGATTTATGCCGCCTCGGCGAGCACATGGAACGGGTTGGCATGGATTGTAGGCGCCGCCCTATACCTTGTCTGCCTCCCCGCGTTGCGTCGCGGTGCATGGAGCCGGGCCACCTCGACCTGA

Protein sequence:
MKKSNNALIVILGTVTLDAVGIGLVMPVLPGLLRDIVHsDSIASHYGVLLALYALMQFLCAPVLGALsDRFGRRPVLLASLLGATIDYAIMATTPVLWILYAGRIVAGITGATGAVAGAYIADITDGEDRARHFGLMSACFGVGMVAGPVAGGLLGAIsLHAPFLAAAVLNGLNLLLGCFLMQESHKGERRPMPLRAFNPVSsFRWARGMTIVAALMTVFFIMQLVGQVPAALWVIFGEDRFRWSATMIGLSLAVFGILHALAQAFVTGPATKRFGEKQAIIAGMAADALGYVLLAFATRGWMAFPIMILLAsGGIGMPALQAMLsRQVDDDHQGQLQGSLAALTSLTSITGPLIVTAIYAASASTWNGLAWIVGAALYLVCLPALRRGAWSRATST


Frame 4 (-1):

start: 5	end: 442	length: 437	ATGGCA	GCGTAA

Gene sequence:
ATGGCAGCAATGAAAGACGGTGAGCTGGTGATATGGGATAGTGTTCACCCTTGTTACACCGTTTTCCATGAGCAAACTGAAACGTTTTCATCGCTCTGGAGTGAATACCACGACGATTTCCGGCAGTTTCTACACATATATTCGCAAGATGTGGCGTGTTACGGTGAAAACCTGGCCTATTTCCCTAAAGGGTTTATTGAGAATATGTTTTTCGTCTCAGCCAATCCCTGGGTGAGTTTCACCAGTTTTGATTTAAACGTGGCCAATATGGACAACTTCTTCGCCCCCGTTTTCACCATGGGCAAATATTATACGCAAGGCGACAAGGTGCTGATGCCGCTGGCGATTCAGGTTCATCATGCCGTCTGTGATGGCTTCCATGTCGGCAGAATGCTTAATGAATTACAACAGTACTGCGATGAGTGGCAGGGCGGGGCGTAA

Protein sequence:
MAAMKDGELVIWDSVHPCYTVFHEQTETFSSLWSEYHDDFRQFLHIYSQDVACYGENLAYFPKGFIENMFFVSANPWVSFTSFDLNVANMDNFFAPVFTMGKYYTQGDKVLMPLAIQVHHAVCDGFHVGRMLNELQQYCDEWQGGA


Frame 5 (-2):

start: 1670	end: 2128	length: 458	ATGGTC	CATTAG

Gene sequence:
ATGGTCGTCGTCATCTACCTGCCTGGACAGCATGGCCTGCAACGCGGGCATCCCGATGCCGCCGGAAGCGAGAAGAATCATAATGGGGAAGGCCATCCAGCCTCGCGTCGCGAACGCCAGCAAGACGTAGCCCAGCGCGTCGGCCGCCATGCCGGCGATAATGGCCTGCTTCTCGCCGAAACGTTTGGTGGCGGGACCAGTGACGAAGGCTTGAGCGAGGGCGTGCAAGATTCCGAATACCGCAAGCGACAGGCCGATCATCGTCGCGCTCCAGCGAAAGCGGTCCTCGCCGAAAATGACCCAGAGCGCTGCCGGCACCTGTCCTACGAGTTGCATGATAAAGAAGACAGTCATAAGTGCGGCGACGATAGTCATGCCCCGCGCCCACCGGAAGGAGCTGACTGGGTTGAAGGCTCTCAAGGGCATCGGTCGACGCTCTCCCTTATGCGACTCCTGCATTAG

Protein sequence:
MVVVIYLPGQHGLQRGHPDAAGSEKNHNGEGHPASRRERQQDVAQRVGRHAGDNGLLLAETFGGGTSDEGLSEGVQDsEYRKRQADHRRAPAKAVLAENDPERCRHLsYELHDKEDSHKCGDDSHAPRPPEGADWVEGSQGHRSTLsLMRLLH


Frame 6 (-3):

start: 1961	end: 2311	length: 350	ATGACC	ATATAG

Gene sequence:
ATGACCACCCAGAGCGCTGCCGGCACCTGTCCTACGAGTTGCATGATAAAGAAGACAGTCATAAGTGCGGCGACGATAGTCATGCCCCGCGCCCACCGGAAGGAGCTGACTGGGTTGAAGGCTCTCAAGGGCATCGGTCGACGCTCTCCCTTATGCGACTCCTGCATTAGGAAGCAGCCCAGTAGTAGGTTGAGGCCGTTGAGCACCGCCGCCGCAAGGAATGGTGCATGCAAGGAGATGGCGCCCAACAGTCCCCCGGCCACGGGGCCTGCCACCATACCCACGCCGAAACAAGCGCTCATGAGCCCGAAGTGGCGAGCCCGATCTTCCCCATCGGTGATGTCGGCGATATAG

Protein sequence:
MTTQSAAGTCPTSCMIKKTVISAATIVMPRAHRKELTGLKALKGIGRRSPLCDsCIRKQPSSRLRPLSTAAARNGACKEMAPNSPPATGPATIPTPKQALMSPKWRARSsPSVMSAI
'''
