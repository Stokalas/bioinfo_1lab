from math import sqrt
import os
from Bio.Seq import Seq

SMTH = ['T', 'C', 'A', 'G']

codons_temp = []

for first in SMTH:
    for second in SMTH:
        for third in SMTH:
            codons_temp.append(f"{first}{second}{third}")

CODONS = codons_temp

# CODONS = [
#     "TTT", "TTC", "TTA", "TTG",
#     "CTT", "CTC", "CTA", "CTG",
#     "ATT", "ATC", "ATA", "ATG",
#     "GTT", "GTC", "GTA", "GTG",

#     "TCT", "TCC", "TCA", "TCG",
#     "CCT", "CCC", "CCA", "CCG",
#     "ACT", "ACC", "ACA", "ACG",
#     "GCT", "GCC", "GCA", "GCG",

#     "TAT", "TAC", "TAA", "TAG",
#     "CAT", "CAC", "CAA", "CAG",
#     "AAT", "AAC", "AAA", "AAG",
#     "GAT", "GAC", "GAA", "GAG",

#     "TGT", "TGC", "TGA", "TGG",
#     "CGT", "CGC", "CGA", "CGG",
#     "AGT", "AGC", "AGA", "AGG",
#     "GGT", "GGC", "GGA", "GGG",
# ]

def file_reader(filename):
    file = open(f"./data/{filename}", 'r')

    first = True

    name = ''
    strand = ''
    for line in file:
        if first:
            name = line.split('>')[1].split(' ')[0].replace('\n', '')
            first = False
            continue
        else:
            strand += line
    
    strand = strand.replace('\n', '')

    file.close()

    return (name, strand)

def get_reverse_compl(strand):
    seq = Seq(strand)
    rev_com = seq.reverse_complement()

    return (strand, rev_com)

def find_start_stop(strand):
    START = "ATG"
    STOP = ["TAA", "TAG", "TGA"]
    length = len(strand)

    rez = []

    for i in range(0,3):
        j = i
        temp = []
        
        while j < length-2:
            codon = strand[j:j+3]
            if codon == START:
                temp.append(('START', j))
            elif codon in STOP:
                temp.append(('STOP', j))
            
            j += 3
    

        searching_stop = False
        start_index = 0
        for codon in temp:
            type, pos = codon
            if searching_stop == False and type == 'START':
                searching_stop = True
                start_index = pos
            elif searching_stop == True and type == 'STOP':
                rez.append((start_index, pos))
                searching_stop = False

    return rez

def filter(pairs):
    rez = []
    for pair in pairs:
        start, stop = pair
        length = stop - start + 2 + 1
        if length % 3 != 0:
            print("err")
        if (stop - start + 2 + 1) >= 100:
            rez.append(pair)
    
    return rez

def calculate_codon_freq(seq, pairs, freq_dict = None):
    if freq_dict == None:
        # print('Set to None!')
        freq_dict = {}
        for codon in CODONS:
            freq_dict[codon] = 0

    for pair in pairs:
        st, en = pair
        st += 3 #ignoring start codon
        while st < en: #ignoring stop codon
            temp = seq[st:st+3]
            freq_dict[temp] += 1

            st += 3

    return freq_dict

def calculate_dicodon_freq(seq, pairs, freq_dict = None):
    if freq_dict == None:
        # print('Set to None!')
        freq_dict = {}
        for codon1 in CODONS:
            for codon2 in CODONS:
                freq_dict[codon1 + codon2] = 0

    for pair in pairs:
        st, en = pair
        st += 3
        while st < (en - 3): #ignoring stop codon
            temp = seq[st:st+6]
            freq_dict[temp] += 1

            st += 3

    return freq_dict

def get_calculated_data_from_file(filename):

    name, data = file_reader(filename)
    (strand, reverse_compl) = get_reverse_compl(data)

    kodons = filter(find_start_stop(strand))
    r_kodons = filter(find_start_stop(reverse_compl))

    c_freq_codon = calculate_codon_freq(strand, kodons)
    d_freq_dicodon = calculate_dicodon_freq(strand, kodons)

    c_freq = calculate_codon_freq(reverse_compl, r_kodons, c_freq_codon)
    d_freq = calculate_dicodon_freq(reverse_compl, r_kodons, d_freq_dicodon)

    return (name, c_freq, d_freq)

most_var = {}

def calc_dist(freq1, freq2):
    if len(freq1) != len(freq2):
        raise Exception('Frequencies\' lists are not of the same length!')
    # if len(d_freq1) != len(d_freq2):
    #     raise Exception('Di-codon frequencies\' lists are not of the same length!')
    
    result = 0
    max = 0
    maxName = ''

    for codon in freq1:
        diff = freq2[codon] - freq1[codon]
        if diff > max:
            max = diff
            maxName = codon


        result += pow(diff, 2)

    # for dicodon in d_freq1:
    #     result += pow(d_freq2[dicodon] - d_freq1[dicodon], 2)

    result = sqrt(result)
    if maxName in most_var:
        most_var[maxName] += max
    else:
        most_var[maxName] = max

    return result

def write_results(filename, dist_matrix):
    result_file = open(filename, 'w')
    result_file.write(f"{len(dist_matrix)}\n")

    for line in dist_matrix:
        line_str = str(line)[1:-1].replace('\'', '').replace(',', '')
        result_file.write(f"{line_str}\n")

    result_file.close()

files = os.listdir("./data")
files.sort()

print(files)


freq_dict = {}
c_dist_matrix = []
d_dist_matrix = []

for f in files:
    name, c_freq, d_freq = get_calculated_data_from_file(f)
    freq_dict[name] = (c_freq, d_freq)

for item in freq_dict:
    codon, dicodon = freq_dict[item]
    temp_c = [item]
    temp_d = [item]

    for item2 in freq_dict:
        codon2, dicodon2 = freq_dict[item2]
        dist_c = "{:.3f}".format(calc_dist(codon, codon2))
        dist_d = "{:.3f}".format(calc_dist(dicodon, dicodon2))
        temp_c.append(dist_c)
        temp_d.append(dist_d)

    c_dist_matrix.append(temp_c)
    d_dist_matrix.append(temp_d)    

write_results('results_codon.txt', c_dist_matrix)
write_results('results_dicodon.txt', d_dist_matrix)

print({key: val for key, val in sorted(most_var.items(), key = lambda ele: ele[1])})



