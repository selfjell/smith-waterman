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
        line_a_end = a_start+(width-line_a_gaps-1)
        line_b_end = b_start+(width-line_b_gaps-1)
        
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
debug = False
#parse optional arguments
for arg in sys.argv:
    if arg.startswith("go"):
        gap_opening = int(arg.split("=")[1])
    elif arg.startswith("ge"):
        gap_extention = int(arg.split("=")[1])
    elif arg.startswith("d"):
        debug = True
print(f'Penalty: Affine\nGap-opening penalty: {gap_opening}\nGap-extention penalty: {gap_extention}')
print(f'Score matrix:{sys.argv[3]}')


#smith-waterman
#initialization
#Aligning matrix
matrix = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
#Insertion matrix
E = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
#Deletion matrix
F = [[0 for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
#Direction matrices
dir_matrix = [["" for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
E_dir = [["" for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]
F_dir = [["" for j in range(len(seq2)+1)] for i in range(len(seq1)+1)]

#coordinates for the highest found score, backtrack from there
highest_score = [0,0,0]

for i in range(1,len(seq1)+1):
    for j in range(1,len(seq2)+1):
        #Match/mismatch
        diag = matrix[i-1][j-1] + scores[indexes[seq1[i-1]]][indexes[seq2[j-1]]]
        #Extending or opening a new gap
        F[i][j] = max(F[i-1][j]-gap_extention,matrix[i-1][j]-gap_opening-gap_extention,0)
        E[i][j] = max(E[i][j-1]-gap_extention,matrix[i][j-1]-gap_opening-gap_extention,0) 
        #Assuming no insertion right after deletions
        if F[i][j] == F[i-1][j]-gap_extention:
            F_dir[i][j] = "U"
        elif F[i][j] == matrix[i-1][j]-gap_opening-gap_extention:
            F_dir[i][j] = "D"

        if E[i][j] == E[i][j-1]-gap_extention:
            E_dir[i][j] = "L"  
        elif E[i][j] == matrix[i][j-1]-gap_opening-gap_extention:
            E_dir[i][j] = "D"

        matrix[i][j] = max(diag, F[i][j],E[i][j],0)
        if matrix[i][j] == diag:
            dir_matrix[i][j] = "D"        
        elif matrix[i][j] == F[i][j]:
            dir_matrix[i][j] = "U"           
        elif matrix[i][j] == E[i][j]:
            dir_matrix[i][j] = "L"
        
        #If it is zero do nothing, it's not part of our local alignment
        #Check if this is the highest score so far
        if matrix[i][j] > highest_score[2]:
            highest_score = [i,j,matrix[i][j]]

#alignments to be printed
align1, align2, alignX = '','',''
i,j = highest_score[0],highest_score[1]

#trim first row and column of arrays
E = [row[1:] for row in E[1:]]
E_dir = [row[1:] for row in E_dir[1:]]
F = [row[1:] for row in F[1:]]
F_dir = [row[1:] for row in F_dir[1:]]

#print all the matrices for debugging
if debug:
    print("Score matrix (main)")
    for row in matrix:
        print(row)
    print("Direction matrix (main)")
    for row in dir_matrix:
        print(row)
    print("E scores")
    for row in E:
        print(row)
    print("E direction:")
    for row in E_dir:
        print(row)
    print("F scores")
    for row in F:
        print(row)
    print("F direction:")
    for row in F_dir:
        print(row)


seq1_start = 1
seq2_start = 1
prev_i, prev_j = 0,0
cur_dir = dir_matrix
score_matrix = matrix
#find the alignment using direction matrix, stops when it hits a zero
while i >= 0 or j >= 0:
    if debug:
        print(f'I: {i} J: {j}')
    if score_matrix[i][j] <= 0:
        #terminate
        seq2_start = prev_j
        seq1_start = prev_i
        break
    dir = cur_dir[i][j]
    prev_i, prev_j = i, j
    #move based on which matrix we are in
    if cur_dir == dir_matrix:
        i -= 1
        j -= 1
    elif cur_dir == E_dir:
        j-=1
    elif cur_dir == F_dir:
        i-=1
    if debug:
        print(f'Elements in list: {dir}')
    #deletion
    if "L" == dir:
        if debug:
            print(f'Deletion for residue:{seq2[j]}')
        align1 = "-" + align1
        align2 = seq2[j] + align2
        alignX  = ' ' + alignX
        cur_dir = E_dir
        score_matrix = E
    #insertion
    elif "U" == dir:
        if debug:
            print(f'Insertion for residue:{seq1[i]}')
        align1  = seq1[i] + align1
        align2 = "-" + align2
        alignX = " " + alignX
        cur_dir = F_dir
        score_matrix = F
    #match/mismatch
    elif "D" == dir:
        if debug:
            print(f'Alignment for: {seq1[i]} and {seq2[j]}')
        align1 = seq1[i] + align1
        align2 = seq2[j] + align2
        if seq1[i] == seq2[j]:
            alignX = "|" + alignX
        elif scores[indexes[seq1[i]]][indexes[seq2[j]]] > 0:
            alignX = "+" + alignX
        else:
            alignX =  " " + alignX
        cur_dir = dir_matrix
        score_matrix = matrix
    

seq2_end = highest_score[1]
seq1_end = highest_score[0]

#print alignment 
print(f'Sequence A:\n{header_seq1}{seq1}')
print(f'Sequence B:\n{header_seq2}{seq2}')
print(f'Score: {highest_score[2]}\nAlignment:')
format_and_print(align1,align2,alignX, seq1_start, seq2_start, seq1_end, seq2_end, width=50)

#Verify the score based on printed alignment
if debug:
    actual_score = 0
    inGap = False
    for i in range(len(align1)):
        if align1[i] == "-" or align2[i] == "-":
            if inGap:
                actual_score -= gap_extention
            else:
                actual_score = (actual_score) - gap_opening - gap_extention
                inGap = True
        else:
            inGap = False
            actual_score += scores[indexes[align1[i]]][indexes[align2[i]]]
    print(f'Actual score from alignment: {actual_score}')

