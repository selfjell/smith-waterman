import sys

def read_score_matrix(path):
    score_matrix = open(path)
    indexes = {}
    scores = []
    for line in score_matrix.readlines():
        #skip metadata
        if line[0] == "#":
            continue
        #put together lookup table
        elif line[0] == " ":
            #exclude newline char
            lst = line[:(len(line)-1)]
            #remove whitespace
            lst = lst.split(" ")
            #remove empty strings
            lst = [x for x in lst if x != ""]
            for i in range(len(lst)):
                indexes[lst[i]] = i
        #append to score matrix
        else:
            #remove whitespace
            lst = line.split(" ")
            #remove letter
            lst = lst[1:]
            scores.append([int(x) for x in lst if x != "" and x != "\n"])

    return indexes, scores

def format_and_print(a,b,x, a_start, b_start, a_end, b_end, width = 50):
    done = False
    #print sequences line by line in a certain width
    while not done:
        if len(a) < width:
            #last line
            print(f'A:\t{a_start}\t{a}\t{a_end}')
            print(f'\t\t{x}')
            print(f'B:\t{b_start}\t{b}\t{b_end}')
            done = True
            break
        #Grab next part of the sequences to be printed
        line_a = a[:width]
        line_b = b[:width]
        line_x = x[:width]
        #Calculate indices at the ends of the lines
        line_a_gaps = len([x for x in line_a if x == "-"])
        line_b_gaps = len([x for x in line_b if x == "-"])
        line_a_end = a_start+(width-line_a_gaps)
        line_b_end = b_start+(width-line_b_gaps)
        
        #Print the line with indexes
        print(f'A:\t{a_start}\t{line_a}\t{line_a_end}')
        print(f'\t\t{line_x}')
        print(f'B:\t{b_start}\t{line_b}\t{line_b_end}')

        a_start = line_a_end+1
        b_start = line_b_end+1
        
        #Discard the parts already printed
        a = a[width:]
        b = b[width:]
        x = x[width:]

if len(sys.argv) < 4:
    raise Exception("Missing arguments")

print("Reading files..")
seq1_file = open(sys.argv[1])
seq2_file = open(sys.argv[2])
header_seq1 = seq1_file.readline()
header_seq2 = seq2_file.readline()
seq1 = seq1_file.read()
seq1 = seq1.replace("\n","")
seq2 = seq2_file.read()
seq2 = seq2.replace("\n","")


indexes, scores = read_score_matrix(sys.argv[3])
gap_opening = 11
gap_extention = 1
#parse optional arguments
for arg in sys.argv:
    if arg.startswith("go"):
        gap_opening = int(arg.split("=")[1])
    elif arg.startswith("ge"):
        gap_extention = int(arg.split("=")[1])
print(f'Penalty: Affine\nGap-opening penalty: {gap_opening}\nGap-extention penalty: {gap_extention}')
print(f'Score matrix:{sys.argv[3]}')


#smith-waterman
#initialization
matrix = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
dir_matrix = [[[] for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
#Aligning matrix
G = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
#Insertion matrix
E = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
#Deletion matrix
F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

#Set first row and column of indel matrices
for i in range(1,len(seq1)+1):
    E[i][0] = -(2*gap_opening) - ((i+1)*gap_extention)
    E[0][i] = float("-inf")
    F[i][0] = float("-inf")
    F[0][i] = -(2*gap_opening) - ((i+1)*gap_extention)

#coordinates for the highest found score, backtrack from there
highest_score = [0,0,0]
for i in range(1,len(seq1)+1):
    for j in range(1,len(seq2)+1):
        G[i][j] = matrix[i-1][j-1] + scores[indexes[seq1[i-1]]][indexes[seq2[j-1]]]
        F[i][j] = max(F[i-1][j]-gap_extention,matrix[i-1][j]-gap_opening-gap_extention)
        E[i][j] = max(E[i][j-1]-gap_extention,matrix[i][j-1]-gap_opening-gap_extention)
        
        matrix[i][j] = max(G[i][j],F[i][j],E[i][j],0)
        if matrix[i][j] == G[i][j]:
            dir_matrix[i][j].append("D")
        if matrix[i][j] == E[i][j]:
            dir_matrix[i][j].append("L")
        if matrix[i][j] == F[i][j]:
            dir_matrix[i][j].append("U")
        #If it is zero do nothing, it's not part of our local alignment
        #Check if this is the highest score so far
        if matrix[i][j] > highest_score[2]:
            highest_score = [i,j,matrix[i][j]]

#alignments to be printed
align1, align2, alignX = '','',''
i,j = highest_score[0],highest_score[1]


seq1_start = 1
seq2_start = 1
prev_i, prev_j = 0,0
#find the alignment using direction matrix, stops when it hits a zero
while i > 0 or j > 0:
    if matrix[i][j] == 0:
        #terminate
        seq2_start = prev_j
        seq1_start = prev_i
        break
    dir = dir_matrix[i][j]
    
    #Keep track of previous in case next one is a zero
    prev_i, prev_j = i, j
    if "L" in dir:
        align1 = "-" + align1
        align2 = seq2[j-1] + align2
        alignX  = ' ' + alignX
        j-=1
    elif "D" in dir:
        align1 = seq1[i-1] + align1
        align2 = seq2[j-1] + align2
        if seq1[i-1] == seq2[j-1]:
            alignX = "|" + alignX
        elif scores[indexes[seq1[i-1]]][indexes[seq2[j-1]]] > 0:
            alignX = "+" + alignX
        else:
            alignX =  " " + alignX
        i -= 1
        j -= 1
    elif "U" in dir:
        align1  = seq1[i-1] + align1
        align2 = "-" + align2
        alignX = " " + alignX
        i -= 1

seq2_end = highest_score[1]
seq1_end = highest_score[0]

#printing results 
print(f'Sequence A:\n{header_seq1}')
print(f'Sequence B:\n{header_seq2}')
print(f'Score: {highest_score[2]}\nAlignment:')
format_and_print(align1,align2,alignX, seq1_start, seq2_start, seq1_end, seq2_end, width=50)

