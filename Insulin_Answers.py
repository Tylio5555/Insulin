# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:22:35 2020

@author: Mickael
"""

import Insulin_Functions as I_F

seq = "AGCCCTCCAGGACAGGCTGCATCAGAAGAGGCCATCAAGCAGGTCTGTTCCAAGGGCCTTTGCGTCAGATCACTGTCCTTCTGCCATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCCGCCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC"


###########################################
# 1. How long is the proteinic sequence ? #
###########################################
coding_seq = I_F.get_coding_seq(seq)
coding_seq_len = I_F.protein_seq_length(coding_seq)

print ("1. How long is the proteinic sequence ?")
print ("The protein is made of "+str(int(coding_seq_len))+" aminoacid.")
print ("")


##################################################
# 2. How many Isoleucin in the coding sequence ? #
##################################################


ile_nb = I_F.get_nb_isoleucin(coding_seq)
#ile_nb_bis = get_nb_isoleucin_biopython(ps)
print ("2. How many Isoleucin in the coding sequence ?")
print ("There is "+str(ile_nb)+" Isoleucine.")
print ("")


##################################################################
# 3. Find the largest sequence of nucleotids which appears twice #
##################################################################

doublon = I_F.largest_sequence_efficient(seq)
print ("3. Find the largest sequence of nucleotids which appears twice")
print ("The largest sequence is: " + doublon)
