# Generates data for creation of GC Content graph, using sliding windows

import string

# Recieves file and window length and returns lists of x and y coordinates for GC content graph
def slidingwindowplot(dna_file, window_length):
    newFile = open("graphData.txt", "w+")
    # Get the actual sequence from the read file and capitalise it
    seq = "".join(dna_file.read().split("\n")[1:]).strip().upper()
    newFile.write(seq)
    newFile.seek(0,0)
    count = 1
    xVals = []
    yVals = []
    i = 1
    while(i > 0):
        # Reads the 'window' from the file
        w = newFile.read(window_length)
        if(not w):
            break
        gc = w.count("G") + w.count("C") # Count frequency of g and c in the window
        gc = "%.2f" % (float(gc) / (len(w) - w.count("N"))) # Frequency as fraction but minus the Ns
        # Stores x values and ensures that the end is shown appropriately
        if len(w) == window_length:
            xVals.append(count * window_length)
            xVals.append("\n")
        #if there is less than window size in the last reading
        else:
            xVals.append(((count - 1) * window_length) + (len(w)))
            xVals.append("\n")
        # Store y values
        yVals.append(gc)
        yVals.append("\n")
        count += 1
    newFile.close()
    dna_file.close()
    # Return lists of x and y values for the graph
    return xVals, yVals


windowSize = int(input("Enter window size:\n"))
fileName = input("Enter file name:\n")
if "." not in fileName:
    fileName += ".txt"
readFile = open(fileName, "r")
xVals, yVals = slidingwindowplot(readFile, windowSize)
# Convert lists to string
xValsStr = "".join(str(x) for x in xVals)
yValsStr = "".join(str(y) for y in yVals)
newFile = open("FinalGraphData.txt", "w+")
newFile.write("X values (window size intervals):\n" + xValsStr + "\n\nY values (GC content):\n")
newFile.write(yValsStr)
newFile.close()
print ("Data for graph has been successfully saved in \'FinalGraphData.txt\' file.")

'''
Example run:

Enter window size:
200
Enter file name:
sample_dna_1
Data for graph has been successfully saved in 'FinalGraphData.txt' file.
'''

