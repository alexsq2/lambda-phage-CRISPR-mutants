#this script will take a text file with separate lines and extract DNA sequences
#then it will match these extracted sequences to a full sequence to determine if it is a deletion
#deletions will be listen in an excel sheet with length of the deletion and any homology at the deletion sites



import xlsxwriter
import glob
wt_region = 'TCTTCAGAACAGGCATTCGCGTCTGAATATCCTTTGGTTCCCATACCGTATAACCATTTGGCTGTCCAAGCTCCGGGTTGATATCAACCTGCAATACGGTGAGCGGTATATCCCAGAACTTCACAACTTCCCTGACAAACCGATATGTCATTGGATGTTCACAACCTGTATCCATGAAAACGTAATGCACGTCTTTACCTGCCCGTCGCTTTTGCTCCATTAGCCAGAGCAAATATGCTGACGTCCTGCCACCGGAGAAACTAACGACATTTATCATGCAGCCCTGTCTCCCCATCTCGCTTTCCACTCCAGAGCCAGTCTCGCTTCGTCTGACCACTTAACGCCACGCTCTGTACCGAATGCCTGTATAAGCTCTAATAGCTCCGCAAATTCGCCCCTACACGCATCCTGCTGGTTGACTGGCCTATTACCACAAAGCCATTCCCGGCAAGGTTAGGAACAACATCCTGCTGCTTTAATGCTGCGGTAAACACACACTTCCAGCTTTCTGCATCCAGCCAGCGACCATGCCATTCAACCTGACGAGAGACGTCACCTAAGCAGGCCCATAGCTTCCTGTTTTGGTCTAAGCTGCGGTTGCGTTCCTGAATGGTTACTACGATTGGTTTGGTTGGGTCTGGAAGGATTTGCTGTACTGCGTGAATAGCGTTTTGCTGATGTGCTGGAGATCGAATTTCAAAGGTTAGTTTTTTCATGACTTCCCTCTCCCCAAATAAAAAGGCCTGCGATTACCAGCAGGCCTGTTATTAGCTCAGTAATGTAGATGGTCATCTTTTAACTCCATATACCGCCAATACCCGTTTCATCGCGGCACTCTGGCGACACTCCTTAAAAACCAGGTTCGTGCTCATCTTTCCTTCCCGTTCTTCCCTGGTAGCAAACCGGTAATACACCGTTCGCCAGACCTTACCTTCGATAACCAGAAGACCTGCCCGTGCCATTTTAGCCGCGGCCTGATTTATGCTGGTTACTGTTGCGCCTGTTAGCGCGGCAACGTCCGGCGCACAGAAGCTATTATGA'
#wt_region is the sequence in which you are looking for deletions

#outputs the reverse complement of a DNA sequence
def reverse_comp(in_sequence):
    complement_list = [['A', 'T'], ['T', 'A'], ['C', 'G'], ['G', 'C']]
    out_sequence = ''
    for x in range(len(in_sequence)):
        for y in range(len(complement_list)):
            if in_sequence[len(in_sequence) -x - 1 ] == complement_list[y][0]:
                out_sequence += complement_list[y][1]
    # print(len(in_sequence),len(out_sequence))
    if len(in_sequence) == len(out_sequence):
        return out_sequence
    else:
        return '0'

#this function will take in a longer sequence 'wt_seq' and a shorter one "candidate' and find the longest stretch
#of the shorter sequence in the longer one
def largest_match(wt_seq,candidate):
    seed_length = 5 #find a match of length 5 somewhere in wt_seq and match the rest of the sequence from that point
    index = 0
    found_positions = [] # the places where matches are found of length "seed length"
    while index + seed_length < len(wt_seq):
        if wt_seq.find(candidate[0:seed_length],index) not in found_positions and wt_seq.find(candidate[0:seed_length],index) !=-1 :
            found_positions.append(wt_seq.find(candidate[0:seed_length],index))
        index += 1

    longest_match_list = [0,0] #length, position in wt sequence
    for position in found_positions:
        # print(position)
        mismatch_found = False
        index2 = 0
        while mismatch_found == False and index2 < len(candidate) and position + index2 < len(wt_seq):
            if candidate[index2] != wt_seq[position + index2]:
                mismatch_found = True
                # print('match length is',index2, 'at position', position)
                if index2 > longest_match_list[0]: #updating the final list if the match is longest one found
                    longest_match_list = [index2,position]
            index2 += 1
        if mismatch_found == False:
            # print('complete match at', position)
            longest_match_list = [len(candidate),position]

    return longest_match_list

#this function uses largest match to identify deletions by finding the largest match in the first half of the
#candidate sequence, cutting it, then finding the largest match in the latter half
#it returns two coordinates which are the coordinates in wt_region of the deleted bases
#requires "largest match function" to be defined above
def deletion_no_start_end (wt_seq, candidate):
    if candidate[0:50] in wt_seq:

        # print(largest_match(wt_seq, candidate)[0])
        wt_second_half = wt_seq[largest_match(wt_seq, candidate)[0] + largest_match(wt_seq, candidate)[1]:]
        candidate_second_half = candidate[largest_match(wt_seq, candidate)[0]:]
        first_coordinate = len(wt_seq) - len(wt_second_half)
        second_coordinate = largest_match(wt_second_half, candidate_second_half)[1] + len(wt_seq) - len(wt_second_half)
        if second_coordinate - first_coordinate > 10 and candidate_second_half[0:50] in wt_second_half and len(candidate_second_half ) >20:
        # if second_coordinate - first_coordinate > 10 :
            return [first_coordinate,second_coordinate]
        else:
            return [0,0]


complement_list = [['A', 'T'], ['T', 'A'], ['C', 'G'], ['G', 'C']]
# perfect_target = 'TTTGTCATTTATGGAGCGTGAGGA' #FORWARD orientation


target_sequence = 'TTTAACGTGCGCTAACTGCGGTCA' #This is the CRISPR target sequence, not required

to_find1_R1 = 'TACCGGAACATCTCGGTAACTGCA'  #extract the target sequence between these two sequences from the R1 file
to_find2_R1 = 'CAAGCAATGCGTGGTGTGCAACC'

rev_to_find1_R1 = 'GGTTGCACACCACGCATTGCTTG'
rev_to_find2_R1 = 'TGCAGTTACCGAGATGTTCCGGTA'

# to_find1_R2 = 'TTTACCCAT'  #extract the target sequence between these two sequences from the R2 file
# to_find2_R2 = 'ACCGCCCG'

# R1_file = open('demultiplex.bc1002--bc1002.hifi_reads.fastq')
# lines_R1 = R1_file.readlines()
# print(len(lines_R1))

seq_list = [] # this list is added to seq_count_list when
total_sequences = 0


total_lines = 0 #number of lines in the file
stranded_sequences = 0 #total number of valid DNA sequences found in the file
good_sequences = 0 #valid deletion sequences
wt_region_rev = reverse_comp(wt_region) #this line will depend on the orientation of the sequence you are looking for
seq_count_list = []  #[sequence,count,read length, deletion length, contains target]
# this will contain the sequences and data this written in the excel sheet

with open("demultiplex.bc1001--bc1001.hifi_reads.fastq") as file1: # enter the file you wish to find the deletion sequences in here
    for line in file1:
        out_list = []
        total_lines += 1
        # print(total_lines)
        if len(line) > 20 and '~' not in line and len(line) < 1000 and '/' not in line:
            # print(line)
            stranded_sequences += 1
            temp_sequence = str(line)[:-2]
            # print(temp_sequence)
            temp_rev_sequence = reverse_comp(temp_sequence)
            # print(temp_rev_sequence)
            deletion_forward = deletion_no_start_end(wt_region,temp_sequence)
            deletion_reverse = deletion_no_start_end(wt_region,temp_rev_sequence)
            out_list.append(deletion_forward)
            out_list.append(deletion_reverse)

            if out_list != [None, None] and out_list != [[0, 0], None] and out_list !=  [None, [0, 0]] and out_list != [[0, 0], [0, 0]]:
                for item in out_list:
                    if item == None or item == [0,0]:
                        out_list.remove(item)
                out_list.append(temp_sequence)
                print(out_list)
                good_sequences += 1
                #find homology between the two deletion sites
                mismatch_found_2 = False
                counter = 0
                micro_homology = 'value not assigned'
                while mismatch_found_2 == False:
                    if wt_region[out_list[0][0] - counter - 1 ] != wt_region[out_list[0][1] - counter -1]:
                        mismatch_found_2 = True
                        if counter == 0:
                            micro_homology = 'none'
                        if counter > 0:
                            micro_homology = wt_region[out_list[0][1] - counter: out_list[0][1]  ]
                    counter += 1

                print('homology length', counter, micro_homology)
                in_list = False
                for item in seq_count_list:
                            if item[2] == temp_sequence:
                                item[0] += 1
                                in_list = True
                if in_list == False:
                    seq_count_list.append([1,out_list[0][0],out_list[0][1], out_list[1],micro_homology])
        # if len(seq_count_list) > 5: #set break for testing small list
        #     break

print(stranded_sequences)
print(total_lines)
print(good_sequences)

#making excel sheet with seq_count_list information
seq_count_list.sort()
seq_count_list.reverse()
print(seq_count_list)


workbook = xlsxwriter.Workbook('myworkbook.xlsx') #name your excel workbook here
worksheet = workbook.add_worksheet('myworksheet') #name the worksheet here, there will only be one worksheet anyway

row = 1
column = 0
#writing the data about each deletion from seq_count_list into the excel sheet
for count, coordinates1, coordinates2 ,sequence1, homology in seq_count_list:
    worksheet.write(row, column, count)
    worksheet.write(row, column + 1, coordinates1)
    worksheet.write(row, column + 2, coordinates2)
    worksheet.write(row, column + 3, sequence1)
    worksheet.write(row, column + 4, homology)



    row += 1


worksheet.write(0,0,'counts')
worksheet.write(0,1,'coordinate 1')
worksheet.write(0,2,'coordinate 2')
worksheet.write(0,3,'sequence')
worksheet.write(0,4,'homology')


workbook.close()

