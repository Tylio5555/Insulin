# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:22:35 2020

@author: Mickael
"""
from Bio.Seq import Seq

###########################################
# 1. How long is the proteinic sequence ? #
###########################################


def get_coding_seq(seq):
    """
    Function to extract the coding sequence
    """
    #Finding the start of the coding sequence
    i=0
    seq_start = None 
    while i < len(seq)-2:
        if seq[i : i+3] == "ATG":
            seq_start = i
            break
        i+=1
        
    if not seq_start:
        return "No starting codon found"

    #Finding the end of the coding sequence
    seq_end = None 
    j=i
    while j < len(seq)-2:
        if seq[j : j+3] in ["TAA", "TGA", "TAG"]:
            seq_end = j+3
            break
        j+=3
    
    if not seq_end:
        return "No Stop codon found"

    return seq[seq_start : seq_end]


def protein_seq_length(seq):
    """
    Function to divide the length by 3 since 
    1 aminoacid is made of 3 nucleotides.
    """
    return len(seq)/3


##################################################
# 2. How many Isoleucin in the coding sequence ? #
##################################################


def get_nb_isoleucin(seq):
    """
    Function to get the number of Isoleucin if the 
    sequence is not translated
    """
    ile_count = 0
    for elt in ["ATT", "ATC", "ATA"]:
        ile_count += seq.count(elt)

    return ile_count


def get_nb_isoleucin_biopython(n_seq):
    """
    Function to get the number of Isoleucin if the sequence
    is translated using 1 single letter aminoacid code.
    """
    aa_seq = str(Seq(n_seq).translate())
    return aa_seq.count("I")


##################################################################
# 3. Find the largest sequence of nucleotids which appears twice #
##################################################################


def largest_sequence(seq, progress=True):
    """
    Brute force algorithm to find the largest double.
    The function return a list if multiple double are
    found having the same size.
    """
    seen_list = []
    d_list = []
    
    i = 0
    while i < len(seq):
        #Indication of function progress
        if progress:
            if i%25 == 0:
                print("Percent done:"+str(round((i/len(seq))*100))+"%")

        j=i+1
        while j < len(seq)+1:
            seq_splice = seq[i:j]
            
            if seq_splice not in seen_list:
                seen_list.append(seq_splice)
                
                #add only if counted twice
                if seq.count(seq_splice) == 2:
                    d_list.append(seq_splice)

            j+=1
        i+=1

    d_list = sorted(d_list, reverse=True, key=len)

    #Case of multiple sequences appearing twice while having the same length
    if len(d_list[1]) == len(d_list[0]):
        out_list = [d_list[0]]
        
        for elt in d_list[1:]:
            if len(elt) ==len(d_list[0]):
                out_list.append(elt)
            else:
                return out_list
        return out_list

    #case of one sequence
    else:
        return d_list[0]


def largest_sequence_efficient(seq):
    """
    Function to find the largest double.
    The function scan for a double from 
    the largest substring size to the smallest.
    """

    frame = len(seq)-1
    while frame > 0:
        i = 0
        while i<len(seq)-frame:
            seq_splice = seq[i:i+frame]
            if seq.count(seq_splice)==2:
                return seq_splice
            i+=1
        frame -= 1
    return "No double found"
    
    
    
    
