
#! /usr/bin/env python3

import sys

complement_dict = {'A':'U', 'G':'C', 'C':'G', 'U':'A'}
bases = ['A', 'a', 'U', 'u', 'C', 'c', 'G', 'g']

def translate_sequence(rna_sequence, genetic_code):

#   Translates a sequence of RNA into a sequence of amino acids.
#
#    Translates `rna_sequence` into string of amino acids, according to the
#    `genetic_code` given as a dict. Translation begins at the first position of
#    the `rna_sequence` and continues until the first stop codon is encountered
#    or the end of `rna_sequence` is reached.
#
#    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
#    an empty string is returned.
#
#    Parameters
#    ----------
#    sequence : str
#        A string representing an RNA sequence (upper or lower-case).
#
#    genetic_code : dict
#        A dictionary mapping all 64 codons (strings of three RNA bases) to
#        amino acids (string of single-letter amino acid abbreviation). Stop
#        codons should be represented with asterisks ('*').
#
#    Returns
#    -------
#    str
#        A string of the translated amino acids


    rna_sequence = rna_sequence.upper()
    protein = []
    if len(rna_sequence) < 3:
        return ''
    elif len(rna_sequence) % 3 != 0:
        while len(rna_sequence) % 3 != 0:
            rna_sequence = rna_sequence[:-1]

    for x in range(0, len(rna_sequence), 3):
        if genetic_code[rna_sequence[x:x+3]] == '*':
            break

        codon = rna_sequence[x:x+3]
        protein += genetic_code[codon]
    return ''.join(protein)

#Made a generator that is able to chunk the sequence into sections by 3. The if and elif
#statements look for sequences less than 3 or sequences ending with extra neucleotides
#that don't code for anything. If there are 1 or 2 extra nucleotides, they are chopped
#off as they are not needed.

def get_all_translations(rna_sequence, genetic_code):
#Get a list of all amino acid sequences encoded by an RNA sequence.
#
#    All three reading frames of `rna_sequence` are scanned from 'left' to
#    'right', and the generation of a sequence of amino acids is started
#    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
#    be in the correct orientation (i.e., no reverse and/or complement of the
#    sequence is explored).
#
#    The function returns a list of all possible amino acid sequences that
#    are encoded by `rna_sequence`.
#
#    If no amino acids can be translated from `rna_sequence`, an empty list is
#    returned.
#
#    Parameters
#    ----------
#    rna_sequence : str
#        A string representing an RNA sequence (upper or lower-case).
#
#    genetic_code : dict
#        A dictionary mapping all 64 codons (strings of three RNA bases) to
#        amino acids (string of single-letter amino acid abbreviation). Stop
#        codons should be represented with asterisks ('*').
#
#    Returns
#    -------
#    list
#        A list of strings; each string is an sequence of amino acids encoded by
#        `rna_sequence`.

    rna_sequence = rna_sequence.upper()
    protein = []
    if len(rna_sequence) < 3:
        return ''
    for x in range(len(rna_sequence)):
        codon = rna_sequence[x:x+3]
        if codon == "AUG":
            tran_seq = translate_sequence(rna_sequence[x:], genetic_code)
            protein.append(tran_seq)
    return protein

#Generator that looks at each codon without translating it. Searches for AUG in the
#RNA sequence; once found, calls the other function previously defined, but specifies
#that is must start at position x. This goes on until another AUG is encountered. Append
#is used since a list is specified for.


def get_reverse(sequence):
#    Reverse orientation of `sequence`.
#
#    Returns a string with `sequence` in the reverse order.
#
#    If `sequence` is empty, an empty string is returned.
#
#    Examples
#    --------
#    >>> get_reverse('AUGC')
#    'CGUA'
#    >>> get_reverse('ATGC')
#    'CGTA'
#
    return (sequence[::-1]).upper()

#Lists the sequence backwards

def get_complement(sequence):
#Get the complement of a `sequence` of nucleotides.
#
#    Returns a string with the complementary sequence of `sequence`.
#
#    If `sequence` is empty, an empty string is returned.
#
#    Examples
#    --------
#    >>> get_reverse('AUGC')
#    'UACG'
#    >>> get_reverse('ATGC')
#    'TACG'

    seq_list = list(sequence.upper())
    seq_list = [complement_dict[base] for base in seq_list]
    return ''.join(seq_list)

#Used list comprehension to create a new list. Start by making the sequence input into
#a list, all uppercase. Take that same list and find its value according to the dictionary
#defined earlier, then make that into a new list. Finally join the new list together
#with no separation in between letters to get the complement.""

def reverse_and_complement(sequence):
#Get the reversed and complemented form of a `sequence` of nucleotides.
#
#    Returns a string that is the reversed and complemented sequence
#    of `sequence`.
#
#    If `sequence` is empty, an empty string is returned.
#
#    Examples
#    --------
#    >>> reverse_and_complement('AUGC')
#    'GCAU'
#    >>> reverse_and_complement('ATGC')
#    'GCAT'

    seq = get_reverse(sequence)
    seq = get_complement(seq)
    return seq

#Takes the previous functions defined and uses them again here

def get_longest_peptide(rna_sequence, genetic_code):
#Get the longest peptide encoded by an RNA sequence.
#
#    Explore six reading frames of `rna_sequence` (the three reading frames of
#    `rna_sequence`, and the three reading frames of the reverse and complement
#    of `rna_sequence`) and return (as a string) the longest sequence of amino
#    acids that it encodes, according to the `genetic_code`
#    If no amino acids can be translated from `rna_sequence` nor its reverse and
#    complement, an empty string is returned.
#
#    Parameters
#    ----------
#    rna_sequence : str
#        A string representing an RNA sequence (upper or lower-case).
#
#    genetic_code : dict
#        A dictionary mapping all 64 codons (strings of three RNA bases) to
#        amino acids (string of single-letter amino acid abbreviation). Stop
#        codons should be represented with asterisks ('*').
#
#    Returns
#    -------
#    str
#        A string of the longest sequence of amino acids encoded by
#        `rna_sequence`.

    protein = get_all_translations(rna_sequence, genetic_code)
    reverse_com = reverse_and_complement(rna_sequence)
    reverse_trans = get_all_translations(reverse_com, genetic_code)

    if protein:
        for x in protein:
            max_length = 0
            peptide = 0
            if len(x) > max_length:
                peptide = x
                max_length = len(x)
        return peptide

    elif reverse_trans:
        for y in reverse_trans:
            max_length = 0
            peptide = 0
            if len(y) > max_length:
                peptide = y
                max_length = len(y)
        return max(reverse_trans, key=len)

    elif not protein or reverse_trans:
        return ''

#If loop that checks for the length of the peptide that is provided. If the peptide is
#longer than the previous one, the value should overwrite the smaller value and
#be stored as peptide. For the reverse_trans function, I'm not sure why I have to
#include max(...,key=len), as otherwise it doesn't work and gives me the shortest
#peptide. If I type min(reverse_trans), the tests come out OK. Possibly don't understand
#max and min fully.

if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
