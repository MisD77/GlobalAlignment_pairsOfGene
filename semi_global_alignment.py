#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 00:41:55 2018

@author: dikshya
"""

# Dikshya Acharya
# CS-327 Computation Bilology Project
# Semi - global Alignment of Pair of Genes

import numpy as np

gap = -1
mismatch = 0
match = 1


#initialization
def semi_global_algo(s1_seq, s2_seq):
    m = len(s2_seq)
    n = len(s1_seq)
    matrix = np.zeros((n + 1, m + 1), dtype=int)
    #Build an alignment matrix
    matrix[0][0] = 0
    
    # Eliminating the terminal gaps from the matrix

    #Building the scoring matrix
    
    for i in range(1, n+1):
        for j in range(1, m+1):
            if s1_seq[i-1] == s2_seq[j-1]:
                score1 = matrix[i-1][j-1] + match
            else:
                score1 = matrix[i-1][j-1] + mismatch

            score2 = matrix[i][j-1] + gap
            score3 = matrix[i-1][j] + gap

            matrix[i][j] = max(score1, score2, score3)

    # create the directional strings
    d_string = build_directional_strings(matrix, n, m)

    # Build alignments using directional strings
    s1_Pos = n - 1 # index of last element in seq1
    s2_Pos = m - 1 # index of last element in seq2
    s1_alignment = ''
    s2_alignment = ''
    dir_pos = 0
    while dir_pos < len(d_string):
        if d_string[dir_pos] == 'D':
            s1_alignment += s1_seq[s1_Pos]
            s2_alignment += s2_seq[s2_Pos]
            s1_Pos -= 1
            s2_Pos -= 1

        elif d_string[dir_pos] == 'V':
            s1_alignment += s1_seq[s1_Pos]
            s2_alignment += '-'
            s1_Pos -= 1
        elif d_string[dir_pos] == 'H':
            s1_alignment += '-'
            s2_alignment += s2_seq[s2_Pos]
            s2_Pos -= 1
            
        dir_pos += 1
        
    s1_alignment = s1_alignment[::-1]
    s2_alignment = s2_alignment[::-1]

    return s1_alignment, s2_alignment


# Finding the paths through the matrix
def build_directional_strings(matrix, n, m):
    d_string = ''
    i = n
    j = m

    while i != 0 or j != 0:
        if i == 0:
            d_string += 'H'
            j -= 1

        elif j == 0:
            d_string += 'V'
            i -= 1

        elif matrix[i][j] == matrix[i][j-1] + gap:
            d_string += 'H'
            j -= 1

        elif matrix[i][j] == matrix[i-1][j] + gap:
            d_string += 'V'
            i -= 1

        else:
            d_string += 'D'
            i -= 1
            j -= 1

    return d_string


def check(seq1, seq2):
    seq1_align, seq2_align = semi_global_algo(seq1, seq2)
    print("The semi global alignment of seq1 is ", seq1_align)
    print("The semi global alignment of seq2 is ", seq2_align)
    print()
    
    
def read_input(filename):
    f = open(filename, 'r')
    pandemic_influenza = f. readline().rstrip('/n')
    M1_gene = f.readline().rstrip('/n')
    s1_align, s2_align = semi_global_algo(pandemic_influenza, M1_gene)
    check(s1_align, s2_align)


def main():
    check('CGCTATAG','CTA')
    check('LILLAVAV','AVAVFTTA')
    check('TCT', 'AATCTATA')
    read_input('global_M1_CodingRegion_test.txt')
    

if __name__ == "__main__":
    main()