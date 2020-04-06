#! /usr/bin/env python3

import argparse
import find_orf
import translate

parser = argparse.ArgumentParser()

parser.add_argument('sequence', type = str, help='RNA sequence that will be'
                    ' translated if there is an open reading frame')

parser.add_argument('-p', '--path', action = 'store_true', help = 'Path to'
                    ' a file containing sequence to be translated. Sequence'
                    ' positional argument must be path to file.')

parser.add_argument('-s', '--start-codons', type = str, nargs = '+',
                    default = ['AUG'], help = 'One or more possible'
                    ' start codons.')

parser.add_argument('-x', '--stop-codons', type = str, nargs = '+',
                    default = ['UAA', 'UAG', 'UGA'], help = 'One or'
                    ' more possible stop codons.')

args = parser.parse_args()

genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}

if args.path:
    sequence = find_orf.parse_sequence_from_path(args.sequence)
else:
    sequence = args.sequence

orf = find_orf.find_first_orf(sequence, start_codons = args.start_codons,
                    stop_codons = args.stop_codons)

translation = translate.translate_sequence(orf, genetic_code)

print('Your (specified) start codon(s): ' + str(args.start_codons))
print('Your (specified) stop codon(s): ' + str(args.stop_codons))
print('Your sequence was: ' + sequence)
print('The open reading frame found was: ' + orf)
print('The translation of the reading frame is:\n' + translation)

