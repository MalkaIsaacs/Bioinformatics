# Assembles a large sequence out of many smaller overlapping sequences, by matching their overlaps

import string


# Gets a file name and returns a dictionary with fragments of sequences
def readDataFromFile(fileName):
    fragmentsFile = open(fileName, "r")
    fragments = {}
    for line in fragmentsFile.readlines():
        key, value = line.strip().split(" ")
        fragments[key] = value.upper()
    fragmentsFile.close()
    return fragments


# Gets a file name and returns the average length of the sequences it contains
def meanLength(fileName):
    fragments = readDataFromFile(fileName)
    totalLength = 0
    totalNum = 0
    for key in fragments.keys():
        totalLength += fragments[key]
        totalNum += 1
    # Returns the average length of a fragment
    return float(totalLength) / totalNum


# Gets 2 fragments and returns the length of their largest real overlap left to right
def getOverlap(left, right):
    overlap = ""
    leftLength = len(left)
    # Compares the end of the left fragment with the begining of right one
    for i in range(leftLength):
        if left[i] == right[0]:
            j = leftLength - i
            if left[i : i + j] == right[: j]:
                overlap += right[: j]
                return overlap
    return overlap


# Returns a dictionary of dictionaries of all possible overlaps of all the fragments
def getAllOverlaps(reads):
    overlaps = {}
    # Make a dictionary with all overlaps for each key
    for key1 in reads.keys():
        overlaps[key1] = {}
        for key2 in reads.keys():
            # Ensure fragment is not overlapped with itself...
            if key1 != key2:
                # Sets the right place in the dictionary to the length of the overlap
                overlaps[key1][key2] = len(getOverlap(reads[key1], reads[key2]))
    return overlaps


# Prints the overlaps in an ordered table, where the rows represent the left fragment and the columns represent the right fragment
def prettyPrint(overlaps):
    print ("Overlaps matrix:\n")
    for key1 in range(len(overlaps) + 1):
        if key1 == 0:
            print ("   ", end = " ")
            for key2 in range(1, len(overlaps) + 1):
                print (" ", key2, end = " ")
            print ("\n", end = " ")
        else:
            print ("% 3d" % key1, end = " ")
            for key3 in range(1, len(overlaps) + 1):
                    if key1 == key3:
                        print (" ", "-", end = " ")
                    else:
                        print ("% 3d" % overlaps[str(key1)][str(key3)], end = " ")
            print ("\n", end = " ")


# Gets the overlaps dictionary of dictionaries and returns the leftmost fragment, with which to start the chain
def findFirstRead(overlaps):
    length = len(overlaps) - 1
    for key1 in range(1, len(overlaps) + 1):
        counter = 0
        for key2 in range(1, len(overlaps) + 1):
            if key1 != key2:
                if overlaps[str(key2)][str(key1)] < 2:
                    counter += 1
                if counter == length:
                    return key1


# Finds largest value in a given dictionary
def findKeyForLargestValue(d):
    maxValue = max(d.values())
    for value in range(1, len(d) + 2):
        if str(value) in d:
            if d[str(value)] == maxValue:
                return value


# Recursive function that returns a list of fragment names in the order in which they make up the sequence
def findOrder(name, overlaps):
    if max(overlaps[str(name)].values()) < 3:
        return [name]
    else:
        nextName = findKeyForLargestValue(overlaps[str(name)])
        return [name] + findOrder(str(nextName), overlaps)


# Assembles the entire sequence/genome
def assembleGenome(readOrder, reads, overlaps):
    seq = ""
    seq += reads[readOrder[0]]
    for fragment in range(1, len(reads)):
        overlap = max(overlaps[readOrder[fragment]].values())
        seq += reads[readOrder[fragment]][overlap:]
    return seq


# Program
fileName = input("Enter name of file:\n")
# If no suffix entered, add it
if "." not in fileName:
    fileName += ".txt"
# Gets the dictionary of reads
fragments = readDataFromFile(fileName)
# Print the raw fragments
print ("\nRaw fragments:\n")
for line in fragments:
    print (line + " " + fragments[line])
print ("\n")
# Gets the dictionary of dictionaries of overlaps
overlaps = getAllOverlaps(fragments)
# Prints the overlaps matrix#
prettyPrint(overlaps)
# Finds leftmost fragment
leftmost = findFirstRead(overlaps)
# Finds with which fragment the largest overlap is
key = findKeyForLargestValue(overlaps[str(leftmost)])
order = findOrder(str(leftmost), overlaps)
print ("\nOrder of fragments:\n", order, "\n")
# Assemble the genome
sequence = assembleGenome(order, fragments, overlaps)
print ("Final joined sequence:\n" + sequence + "\n")


'''
Example run:

Enter name of file:
genome_assembly

Raw fragments:

1 GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGC
2 CTTTACCCGGAAGAGCGGGACGCTGCCCTGCGCGATTCCAGGCTCCCCACGGG
3 GTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTCGTGAACACATCAGT
4 TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCG
5 CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC
6 TGACAGTAGATCTCGTCCAGACCCCTAGCTGGTACGTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGT


Overlaps matrix:

      1   2   3   4   5   6 
   1   -   1   0   0   1  29 
   2  13   -   1   0  21   0 
   3   0   0   -   1   0   1 
   4   1  17   1   -   2   0 
   5  39   1   0   0   -  14 
   6   0   0  43   1   0   - 
 
Order of fragments:
 ['4', '2', '5', '1', '6', '3'] 

Final joined sequence:
TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCGGCTGCCCTGCGCGATTCCAGGCTCCCCACGGGTAGATCTCA
GTAGATCTCGTCCAGACCCCTAGCTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTTCTTCAGTAGAAAATTGTTTTTTT
CTTCCAAGAGGTCGGAGTCGTGAACACATCAGT
'''
