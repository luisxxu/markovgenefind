import math
import matplotlib.pyplot as plt

# convert fasta file into string
db = ""
with open("/Users/luisxu/Documents/CSE 182/A4/chrA.fasta", "r") as fastafile:
        filecontent = fastafile.readlines()
        for line in filecontent:
            if line.startswith(">"):
                continue
            db += line.strip()

# find island starts and ends

islands = []
starts = []
ends = []

with open("/Users/luisxu/Documents/CSE 182/A4/chrA.islands", "r") as islandfile:
    filecontent = islandfile.readlines()
    for line in filecontent:
        start = ""
        end = ""
        pos = 0
        while line[pos] != ' ':
            start += line[pos]
            pos += 1
        pos += 1
        while line[pos] != '\n':
            end += line[pos]
            pos += 1
        starts.append(int(start))
        ends.append(int(end))
        islands.append(int(start))
        islands.append(int(end))


# split fastafile string into cpg and notcpg strings
cpg = ""
notcpg = ""

j = 0
island = False
for i, char in enumerate(db):
    if j < len(islands):
        if i != islands[j]:
            
            if island:
                cpg += char
            else:
                notcpg += char
        else:
            j += 1
            if island:
                island = False
                notcpg += char
            else:
                island = True
                cpg += char

# compute dinucleotide frequencies for cpg/notcpg strings
cpgfreqs = [0] * 16
notcpgfreqs = [0] * 16
difreqnames = ["AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"]
for i in range(len(cpg) - 1):
    index = 0
    curr = cpg[i:i+2]
    if curr[0] == "A":
        index += 0
    elif curr[0] == "C":
        index += 4
    elif curr[0] == "G":
        index += 8
    else:
        index += 12
    if curr[1] == "A":
        index += 0
    elif curr[1] == "C":
        index += 1
    elif curr[1] == "G":
        index += 2
    else:
        index += 3
    
    cpgfreqs[index] += 1

for i in range(len(notcpg) - 1):
    index = 0
    curr = notcpg[i:i+2]
    if curr[0] == "A":
        index += 0
    elif curr[0] == "C":
        index += 4
    elif curr[0] == "G":
        index += 8
    else:
        index += 12
    if curr[1] == "A":
        index += 0
    elif curr[1] == "C":
        index += 1
    elif curr[1] == "G":
        index += 2
    else:
        index += 3
    
    notcpgfreqs[index] += 1

total = 0
for i in range(16):
    cpgfreqs[i] = cpgfreqs[i] / (len(cpg) - 1)
    total += cpgfreqs[i]
    notcpgfreqs[i] = notcpgfreqs[i] / (len(notcpg) - 1)

w = 400

potential = []

for window in range(len(db) - w):
    cpgtotal = 0
    notcpgtotal = 0
    for ind in range(w):
        index = 0
        curr = db[window + ind: window + ind + 2]
        if curr[0] == "A":
            index += 0
        elif curr[0] == "C":
            index += 4
        elif curr[0] == "G":
            index += 8
        else:
            index += 12
        if curr[1] == "A":
            index += 0
        elif curr[1] == "C":
            index += 1
        elif curr[1] == "G":
            index += 2
        else:
            index += 3
        cpgtotal += math.log(cpgfreqs[index])
        notcpgtotal += math.log(notcpgfreqs[index])

    potential.append(cpgtotal - notcpgtotal)

pcpgstarts = []
pcpgends = []
pncpgstarts = []
pncpgends = []

island = False
pottotal = 0
potnum = 0
pncpgstarts.append(0)
for i, pot in enumerate(potential):
    if pot > 40 and island == False:
        island = True
        pcpgstarts.append(i)
        pncpgends.append(i)
    elif pot < 40 and island == True:
        island = False
        pcpgends.append(i)
        pncpgstarts.append(i)

if len(pncpgstarts) > len(pncpgends):
    pncpgends.append(len(db))

pcpgresult = []

tp = 0
fp = 0
fn = 0
tn = 0

for j, pstart in enumerate(pcpgstarts):
    overlap = 10000000000
    for i, iend in enumerate(ends):
        if (starts[i] < pcpgends[j] and pstart < iend):
            if abs(iend-pstart) < overlap:
                overlap = abs(iend - pstart)
    if (4 * overlap) > (pcpgends[j] - pstart) and overlap != 10000000000:
        pcpgresult.append("true positive")
        tp += 1
    else:
        pcpgresult.append("false positive")
        fp += 1

pncpgresult = []
for j, pstart in enumerate(pncpgstarts):
    overlap = 10000000000
    for i, iend in enumerate(ends):
        if (starts[i] < pncpgends[j] and pstart < iend):
            if abs(iend-pstart) < overlap:
                overlap = abs(iend - pstart)
    if (4 * overlap) > (pncpgends[j] - pstart) and overlap != 10000000000:
        pncpgresult.append("false negative")
        fn += 1
    else:
        pncpgresult.append("true negative")
        tn += 1




with open("/Users/luisxu/Documents/CSE 182/A4/predictedregions.txt", "w") as outputfile:
    for i in range(len(pcpgstarts)):
        outputfile.write(str(pncpgstarts[i]) + "\t" + str(pncpgends[i]) + "\t" + pncpgresult[i] + "\n")
        outputfile.write(str(pcpgstarts[i]) + "\t" + str(pcpgends[i]) + "\t" + pcpgresult[i] + "\n")
    outputfile.write(str(pncpgstarts[len(pcpgstarts)]) + "\t" + str(pncpgends[len(pcpgstarts)]) + "\t" + pncpgresult[len(pcpgstarts)] + "\n")

print("Statistics of predicted islands:\nTrue Positive Rate: " + str(tp / (len(pcpgresult) + (len(pncpgresult)))))
print("False Positive Rate: " + str(fp / (len(pcpgresult) + (len(pncpgresult)))))
print("True Negative Rate: " + str(tn / (len(pcpgresult) + (len(pncpgresult)))))
print("False Negative Rate: " + str(fn / (len(pcpgresult) + (len(pncpgresult)))))