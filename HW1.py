#!/usr/bin/env python 
# -*- coding:utf-8 -*-

#!/usr/bin/python
__author__ = "Gaoqianxue Liu"
__email__ = "qianxue.liu@yale.edu"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "1.0.0"
### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse
### This is one way to read in arguments in Python.
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False,
default=-2)
parser.add_argument('-e', '--extgap', help='extension gap',
required=False, default=-1)
args = parser.parse_args()


# gap = '' -> we define a gap penalty function at the end of the code
# (so the program is suitable for affine penalties too)


def runSW(inputFile, scoreFile, openGap, extGap):
    #read in files
    Substi = []

    for row in open(scoreFile):
        Substi.append([x for x in row.split()])
    Substi.pop(0)
    Substi.pop()
    for i in range(len(Substi) - 1):
        Substi[i].pop(0)

    with open(scoreFile) as f:
        alphabet: str = f.readline().strip().replace(" ","")

    with open(inputFile) as f:
        seqA = f.readline().strip('\n')
        seqB = f.readline().strip('\n')

    #print sequences
    print("-----------")
    print("|Sequences|")
    print("-----------")
    print("sequence1")
    print(seqA)
    print("sequence2")
    print(seqB)

    #print score matrix
    print("--------------")
    print("|Score Matrix|")
    print("--------------")
    rows = len(seqA) + 1
    cols = len(seqB) + 1

    score_matrix, start_pos, max_score, path_matrix = create_score_matrix(rows, cols, seqA, seqB, Substi, openGap, extGap, alphabet)
    score_matrix = list(map(list, zip(*score_matrix)))
    score_matrix.insert(0, list(seqA))
    for i in range(2, len(seqB)+2):
        score_matrix[i].insert(0, list(seqB)[i - 2])
    print_matrix(score_matrix)

    print("----------------------")
    print("|Best Local Alignment|")
    print("----------------------")

    #Printing the alignements
    print('Alignment Score:{}'.format(max_score))
    print("Alignment Results:")
    seqA_aligned, seqB_aligned = traceback(seqA, seqB, path_matrix, start_pos, score_matrix)
    assert len(seqA_aligned) == len(seqB_aligned), 'aligned strings are not the same size'
    print_alignment(seqA, seqB, seqA_aligned, seqB_aligned)
    return (0)
### Run your Smith-Waterman Algorithm


def traceback(seqA, seqB, path_matrix, start_pos, score_matrix):
    '''The function does the traceback based on the starting position
    and on the instructions contained in the path_matrix,
    it also displays the move done at each step '''

    x, y = start_pos
    aligned_seqA = []
    aligned_seqB = []

    while path_matrix[x][y] != [0, 'NULL']:
        d, direction = path_matrix[x][y][0], path_matrix[x][y][1]
        if direction == 'DIAG':
            assert d == 1, 'path_matrix wrongly constructed !'
            aligned_seqA.append(seqA[x - 1])
            aligned_seqB.append(seqB[y - 1])
            x -= 1
            y -= 1
        elif direction == 'UP':
            for c in range(d):
                aligned_seqA.append(seqA[x - 1])
                aligned_seqB.append('-')
                x -= 1
        elif direction == 'LEFT':
            for c in range(d):
                aligned_seqA.append('-')
                aligned_seqB.append(seqB[y - 1])
                y -= 1

    return ''.join(reversed(aligned_seqA)), ''.join(reversed(aligned_seqB))


def create_score_matrix(rows, cols, seqA, seqB, Substi, openGap, extGap, alphabet):
    # create the score_matrix and the path_matrix
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]
    path_matrix = [[[0, 'NULL'] for col in range(cols)] for row in range(rows)]

    max_score = 0
    max_pos = None
    for i in range(1, rows):
        for j in range(1, cols):
            score, antecedent = calc_score(score_matrix, i, j, seqA, seqB, Substi, openGap, extGap, alphabet)
            if score > max_score:
                max_score = score
                max_pos = (i, j)

            score_matrix[i][j], path_matrix[i][j] = score, antecedent

    assert max_pos is not None, 'No maximum found'

    return score_matrix, max_pos, max_score, path_matrix



def calc_score(score_matrix, x, y, seqA, seqB, Substi, openGap, extGap, alphabet):
    # score_matrix construction,


    similarity = Substitution_score(Substi, x, y, seqA, seqB, alphabet)

    same_row = [(score_matrix[x][y - l] + openGap + gap_penalty(l-1, extGap)) for l in range(1, x + 1)]
    same_col = [(score_matrix[x - k][y] + openGap + gap_penalty(k-1, extGap)) for k in range(1, x + 1)]

    up_score = max(same_col)
    left_score = max(same_row)

    diag_score = score_matrix[x - 1][y - 1] + similarity
    pos_max_up = first_pos_max(same_col)
    pos_max_left = first_pos_max(same_row)

    score = max(0, diag_score, up_score, left_score)

    if score == diag_score:
        antecedent = [1, 'DIAG']
        return score, antecedent
    elif score == up_score:
        antecedent = [pos_max_up + 1, 'UP']
        return score, antecedent
    elif score == left_score:
        antecedent = [pos_max_left + 1, 'LEFT']
        return score, antecedent
    else:
        return score, [0, 'NULL']


def Substitution_score(Substi, x, y, seqA, seqB, alphabet):
    # score between seqA[x-1] and seqB[y-1]
    a_i = alphabet.index(seqA[x - 1])
    b_j = alphabet.index(seqB[y - 1])

    return int(Substi[a_i][b_j])

def alignment_string(aligned_seqA, aligned_seqB):
    # join aligned string
    alignment_string = []

    for base1, base2 in zip(aligned_seqA, aligned_seqB):  # loop that runs over both sequences simultaneously
        if base1 == base2:
            alignment_string.append('|')
        else:
            alignment_string.append(' ')

    return ''.join(alignment_string)



def print_alignment(seqA, seqB, seqA_aligned, seqB_aligned):
    #print alignment
    alignment_str = alignment_string(seqA_aligned, seqB_aligned)
    seqaa = seqA_aligned.replace("-","")
    if seqaa in seqA:
        stra = seqA.replace(seqaa, '*')
    seqba = seqB_aligned.replace("-", "")
    if seqba in seqB:
        strb = seqB.replace(seqba, '*')
    stra = stra.split("*")
    strb = strb.split("*")

    m = max(len(stra[0]), len(strb[0]))
    x = max(len(stra[1]), len(strb[1]))
    print('{n:>{m}}({p:^s}){q:<{x}}'.format(n=stra[0], m=m, p=seqA_aligned, q=stra[1], x=x))
    print('{n:>{m}}{p:^s}{q:<{x}}'.format(n="", m=m + 1, p=alignment_str, q="", x=x+1))
    print('{n:>{m}}({p:^s}){q:<{x}}'.format(n=strb[0], m=m, p=seqB_aligned, q=strb[1], x=x))



def print_matrix(matrix):  # print score_matrix
    print("\t\t"+"\t".join(matrix[0]))
    print("\t"+'\t'.join([str(item) for item in matrix[1]]))
    print('\n'.join(['\t'.join([str(item) for item in row]) for row in matrix[2:]]))


def first_pos_max(list): #first position
    maxi = max(list)
    return [i for i, j in enumerate(list) if j == maxi][0]


def gap_penalty(k, extGap): #gap extend penalty
    return extGap * k

runSW(args.input, args.score, args.opengap, args.extgap)