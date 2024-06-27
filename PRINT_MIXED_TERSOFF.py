#====================================================
# Prints all Tersoff parameters for a system of N
# elements.
#
# You have to write an input file with the pure element
# Tersoff parameters (e.g., if your system is SiGe you
# need to specify the pure Tersoff parameters for Si
# and Ge) and include the mixing parameter Chi if you
# have it (e.g., Chi_{Si-Ge} = 1.00061)
#
# This program takes that input file and outputs all
# mixing parameters in the 3-atom bond described in the 
# Tersoff potential.
# 
# Example of input file for Si-Ge-H system:
#==========================================
#   H         3.0 1.0 0.0 0.       1.     1.       1.      4.         1.98    43.531    0.90    0.10   3.7879  86.7120
#   Si(D)     3.0 1.0 0.0 1.0039e5 16.217 -0.59825 0.78734 1.1000e-6  1.7322  471.18    2.85    0.15   2.4799  1830.8
#   Ge        3.0 1.0 0.0 1.0643e5 15.652 -0.43884 0.75627 9.0166e-7  1.7047  419.23    2.95    0.15   2.4451  1769.0
#
#   Si(D) Ge 1.00061
#   Si(D) H 0.78
#   Ge H 0.76
#==========================================
# ^^^ The last three lines are the Chi values which has the format:
#   (element 1) (element 2) Chi_{element1-element2} (dont worry about which element you list first)
# The first three lines are the Tersoff parameters for the pure elements in LAMMPS format (notice the
# "1.0 3.0" in the beginning; this specifies a Tersoff_2 style potential (see LAMMPS documentation on Tersoff).
#
#
# Usage:
# python PRINT_MIXED_TERSOFF.py --input example1
# Author: Jesper Kristensen, Cornall University 2013
#====================================================
#!/bin/python3
import argparse
import re, os
import numpy as np

#HELPER PARSER FUNCTION
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#PARSES INPUT FILE
def parse_xyz_file(inputfile):

    with open(inputfile, 'r') as f:
        symbols = []
        tersoff_params = [] #convert xyz file to a matrix
        print('0000 Parsing Tersoff parameters...')
        for line in f.readlines():
            # Ignore comments and empty lines:
            li = line.strip()
            if li.startswith("#") or not li:
                continue
            extracted = [re.sub(r'\n', r'', line_all) for line_all in filter(None, line.split(' '))]
            if not is_number(extracted[1]):
                continue
            symbols.append(extracted[0])

            if len(extracted[1:]) == 12:
                print('EEEE     You did not supply enough Tersoff parameters? Did you remember the initial two "3.0 1.0" numbers? Read the LAMMPS doc for the format.')
                quit()
            if len(extracted[1:]) != 14:
                print('EEEE     You did not supply enough Tersoff parameters? Read the LAMMPS doc for the format.')
                quit()
            else:
                # Number of Tersoff parameters OK:
                tersoff_params.append(extracted[1:])
        print('0000 done.')

        # NOW GET CHI-VALUES (MIXING PARAMETER):
        print('0000 Parsing CHI values (mixing parameter)...')
        n = len(symbols)
        f.seek(0) # reset reading of file
        CHIS_IND = []
        CHIS = []
        for line in f.readlines():
            # Ignore comments and empty lines:
            li = line.strip()
            if li.startswith("#") or not li:
                continue
            extracted = [re.sub(r'\n', r'', line_all) for line_all in filter(None, line.split(' '))]
            if not is_number(extracted[1]):
                # CHI VALUE:
                ind1 = symbols.index(extracted[0]) # element 1
                ind2 = symbols.index(extracted[1]) # element 2
                CHIS_IND.append([min(ind1, ind2), max(ind1, ind2)]) # order: lowest index first
                CHIS.append(float(extracted[2]))

        CHI_MATRIX = np.ones(shape=(n, n))
        for i in range(0, n):
            for j in range(i + 1, n):
                if [i, j] in CHIS_IND: # we have ordered "CHIS"
                    CHI_MATRIX[i][j] = CHIS[CHIS_IND.index([i, j])]
                    CHI_MATRIX[j][i] = CHI_MATRIX[i][j] # Symmetric
                    print('0000     Chi-value of ', str(CHI_MATRIX[i][j]), ' found between element ', str(symbols[i]), ' and ', str(symbols[j]))
                else:
                    print('0000     ****Note: I do not have any Chi-parameter for the (', str(symbols[i]), ' and ', str(symbols[j]), ') combination: Assuming it to be 1.')
        print('0000 done.')

    return symbols, tersoff_params, CHI_MATRIX

def compute_RD_mixing(R_lammps1, R_lammps2, D1, D2):
    # element 1
    Rp1 = R_lammps1 - D1
    Sp1 = R_lammps1 + D1
    # element 2
    Rp2 = R_lammps2 - D2
    Sp2 = R_lammps2 + D2
    # Use rule of mixing from Tersoff 1989 paper:
    Rptot = np.sqrt(Rp1 * Rp2)
    Sptot = np.sqrt(Sp1 * Sp2)
    # Convert back to LAMMPS parameters:
    R_lammps = (Sptot + Rptot) / 2.0
    D = (Sptot - Rptot) / 2.0
    return R_lammps, D

def compute_two_body_params(tersoff_params, CHI_MATRIX, i, j):
    RES = tersoff_params[i]
    RES = RES[0:8]
    RES1 = [float(x) for x in tersoff_params[i]] # element 1
    RES2 = [float(x) for x in tersoff_params[j]] # element 2

    # Naming constants in LAMMPS notation
    # comments to the far right gives the variable name in terms of the Tersoff 1989 paper:
    lambda2 = (RES1[8] + RES2[8]) / 2.0 # in Tersoff 1989 paper: mu
    B = CHI_MATRIX[i][j] * np.sqrt(RES1[9] * RES2[9]) # in Tersoff 1989 paper: B
    R_lammps, D = compute_RD_mixing(RES1[10], RES2[10], RES1[11], RES2[11]) # in Tersoff 1989 paper: He calls R_lammps for R and D for S
    lambda1 = (RES1[12] + RES2[12]) / 2.0 # in Tersoff 1989 paper: lambda
    A = np.sqrt(RES1[13] * RES2[13]) # in Tersoff 1989 paper: A

    RES.extend([str(lambda2), str(B), str(R_lammps), str(D), str(lambda1), str(A)])
    return RES

def compute_three_body_params(tersoff_params, CHI_MATRIX, i, k):
    # relevant: elements i and k
    RES = tersoff_params[i]
    RES = RES[0:7]
    RES1 = [float(x) for x in tersoff_params[i]] # element 1
    RES2 = [float(x) for x in tersoff_params[k]] # element 2

    beta = 0.0
    lambda2 = 0.0
    B = 0.0
    R_lammps, D = compute_RD_mixing(RES1[10], RES2[10], RES1[11], RES2[11]) # in Tersoff 1989 paper: He calls R_lammps for R and D for S
    lambda1 = 0.0
    A = 0.0
    RES.extend([str(beta), str(lambda2), str(B), str(R_lammps), str(D), str(lambda1), str(A)])
    return RES

#==========
# Main
#=========
def main():
    #==== Parse command line: ====
    parser = argparse.ArgumentParser(description='@Jesper Kristensen, Cornell 2013: Compute and output Tersoff mixing parameters.')
    parser.add_argument('--input', dest='inputfile', type=str, help='input file with pure element Tersoff parameters (and chi mixing values)')
    args = parser.parse_args()

    infile = args.inputfile
    if not os.path.isfile(infile):
        print('')
        print('EEEE Input the file containing pure Tersoff coefficients')
        print('')
        quit()
    #==== 
    
    #==== EXTRACT SYMBOL NAMES AND THEIR COEFFICIENTS:
    print('')
    print('0000 Parsing input file: ', str(infile))
    symbols, tersoff_params, CHI_MATRIX = parse_xyz_file(infile)
    print('0000 Done parsing input file')
    print('')
    print('0000 Outputting all interactions in LAMMPS Tersoff_2 format:')
    n = len(symbols)
    outfile = ''
    for sym in symbols:
        outfile += str(sym)
    outfile += '.tersoff'
    with open(outfile, 'w') as f:
        #==== interactions:
        # think of 3 bins which can each take n values:
        # [n]x[n]x[n]
        # for example:
        # n=1:  [1]x[1]x[1] = 1 combination (e.g., SiSiSi)
        # n=2:  [2]x[2]x[2] = 8 combinations (e.g., {SiSiSi,SiSiC,SiCSi,...})
        # n=3:  [3]x[3]x[3] = 27 combinations (e.g., {SiSiSi, SiSiC, ..., SiSiGe, ...}
        # etc.
        for i in range(0, n): # bin 1
            for j in range(0, n): # bin 2
                for k in range(0, n): # bin 3
                    if i == j:
                        if j == k:
                            # PURE ATOM (all indices equal)
                            to_print = str(symbols[i]) + ' ' + str(symbols[j]) + ' ' + str(symbols[k]) + '     ' + ' '.join(tersoff_params[i])
                            print(to_print)
                            f.write('%s\n' % str(to_print))
                        else:
                            params = compute_three_body_params(tersoff_params, CHI_MATRIX, i, k)
                            to_print = str(symbols[i]) + ' ' + str(symbols[j]) + ' ' + str(symbols[k]) + '     ' + ' '.join(params)
                            print(to_print)
                            f.write('%s\n' % str(to_print))
                    else:
                        if j == k:
                            # TWO-BODY INTERACTION (second index repeated but different from first index)
                            params = compute_two_body_params(tersoff_params, CHI_MATRIX, i, j)
                            to_print = str(symbols[i]) + ' ' + str(symbols[j]) + ' ' + str(symbols[k]) + '     ' + ' '.join(params)
                            print(to_print)
                            f.write('%s\n' % str(to_print))
                        else:
                            params = compute_three_body_params(tersoff_params, CHI_MATRIX, i, k)
                            to_print = str(symbols[i]) + ' ' + str(symbols[j]) + ' ' + str(symbols[k]) + '     ' + ' '.join(params)
                            print(to_print)
                            f.write('%s\n' % str(to_print))

    print('')
    print('')
    print('0000 result in file ', str(outfile))

if __name__ == "__main__":
    main()
