#This program calculates the most commonly occuring k-mers in the file specified by the user and their frequency in it

import string

fileName = input("Enter the name of the file, with the correct suffix.\n")
# If file name does not contain a suffix, default (.txt) is assumed
if "." not in fileName:
    fileName = fileName + '.txt'
file = open(fileName, 'r')
seq = file.read().split("\n")[1]
k = int(input("Enter \"word\" size:\n"))

words = {}
# Find all possible k-mers and count the frequency of each
for char in range(len(seq) - k + 1):
    word = seq[char : char + k]
    if word not in words:
        words[word] = seq.count(word)
        
# Assign 0 to all words that appear less than the most common.
for w in words:
    for val in words:
        if words[val] < words[w]:
            words[val] = 0

# Print all words that appear the most
for w in words:
    if words[w] != 0:
        print (w)
        num = words[w]
print ("Each of these \"words\" appears %d times in the sequence." % num)

'''
Example run:
Enter the name of the file, with the correct suffix.
dna.txt
Enter "word" size
4
atca
atga
tgat
Each of these "words" appears 11 times in the sequence.
'''
