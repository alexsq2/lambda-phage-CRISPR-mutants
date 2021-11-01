#script only works for gene J targets (398) unless "to_find" sequences are switched out
# skeleton of script that checks R1 and R2 files to see if they match

# take two sequences that are supposed to be reverse complement of each other
import xlsxwriter
import glob

def reverse_comp(in_sequence): #to do - implement this function to simplify
    complement_list = [['A', 'T'], ['T', 'A'], ['C', 'G'], ['G', 'C']]
    out_sequence = ''
    for x in range(len(in_sequence)):
        for y in range(len(complement_list)):
            if in_sequence[len(in_sequence) -x - 1 ] == complement_list[y][0]:
                out_sequence += complement_list[y][1]
    if len(in_sequence) == len(out_sequence):
        return out_sequence
    else:
        print('invalid input sequence')


complement_list = [['A', 'T'], ['T', 'A'], ['C', 'G'], ['G', 'C']]
perfect_target = 'CTGCGGGCGGTTTTGTCATTTATGGAGCGTGAGGAATGGGTAA' #FORWARD orientation



perfect_target_list = []
for x in range(len(perfect_target)):
    perfect_target_list += perfect_target[x]
print(perfect_target_list)


to_find1_R1 = 'CT'  #extract the target sequence between these two sequences from the R1 file
to_find2_R1 = 'AA'

to_find1_R2 = 'TT'  #extract the target sequence between these two sequences from the R2 file
to_find2_R2 = 'AG'

workbook = xlsxwriter.Workbook('cas12a_escapers_J1_R1R2 3.xlsx')  # naming the output excel sheet

# selecting files to analyze, here it's all the Cas12a ones and the control
total_file_list = []
total_file_list += glob.glob('*MASJ*')

# total_file_list += glob.glob('*lys*')
print(total_file_list)
print(len(total_file_list) , 'files gathered')

file_list_length_half = int(len(total_file_list) / 2)
print(file_list_length_half)

total_file_list_pairs = [] #this list will contain pairs of files ordered by their order in total file list
#total file list must be in order or the wrong files will be matched!!

for x in range(file_list_length_half):
    total_file_list_pairs.append([total_file_list[2*x]] + [total_file_list[2*x + 1]])

print(total_file_list_pairs)

for file_pair in total_file_list_pairs:


    mismatch_list = []  # list for making the mismatch distribution for each file
    for a in range(len(perfect_target)):
        mismatch_list.append([0, 0, 0, 0, 0, 0])  # mutation/deletion types: A,T,C,G,deletion,multiple mismatches where A means previous base mutated to A

    sequence_lines = 0    #this represents the total number of valid sequences found, valid = R1R2 match
    wild_type_sequences = 0
    single_mm_sequences = 0
    indel_sequences = 0
    multi_mm_sequences = 0

    # make a list with values for each file
    # total counts, WT, 1MM, 1deletion, locations of mm/del, 0,1,2,3,4,5... ,23

    seq_count_list = []   #defining this list now for use much later on

    R1_file = open(file_pair[0])
    R2_file = open(file_pair[1])

    lines_R1 = R1_file.readlines()
    lines_R2 = R2_file.readlines()

    print(len(lines_R1))
    print(len(lines_R2))

    line_count_r1 = 0
    line_count_r2 = 0

    same_length_sequences = 0
    matching_sequences = 0

    for line_num in range(len(lines_R1)):  # lines_R1 and R2 SHOULD have the same length
        temp_sequence_R1 = ''
        temp_sequence_R2 = ''
        if lines_R1[line_num].find(to_find1_R1) != -1 and lines_R1[line_num].find(to_find2_R1) != -1:
            temp_sequence_R1 = lines_R1[line_num][lines_R1[line_num].find(to_find1_R1,30,35) : lines_R1[line_num].find(to_find2_R1,69) +2 ]   #72 here because there are many occcurences of "to_find_R1" in the line
            line_count_r1 += 1
            # print('true')

            # print(lines_R1[line_num].find(to_find1_R1))
            # print(lines_R1[line_num].find(to_find2_R1, 72) )

        if lines_R2[line_num].find(to_find1_R2) != -1 and lines_R2[line_num].find(to_find2_R2) != -1:
            temp_sequence_R2 = lines_R2[line_num][lines_R2[line_num].find(to_find1_R2,31) : lines_R2[line_num].find(to_find2_R2,70) + 2]
            # print(temp_sequence_R2)
            # print(temp_sequence_R2)
            line_count_r2 += 1

        # print(len(temp_sequence_R1),len(temp_sequence_R2))
        # print(temp_sequence_R2)

        if len(temp_sequence_R1) == len(temp_sequence_R2) and len(temp_sequence_R1) > 20 and len(temp_sequence_R2) > 20:  # if they are same length, rev complement R2, match with R1
            same_length_sequences += 1
            R2_rev_complement = reverse_comp(temp_sequence_R2)

            if R2_rev_complement == temp_sequence_R1:
                matching_sequences += 1
                sequence_lines += 1
                # print(temp_sequence_R1)
                temp_sequence = temp_sequence_R1  # the final sequence used for analysis and comparison
                # print(temp_sequence)

                if len(temp_sequence) == len(perfect_target) - 1 :  # counting single indels, multiple deletions are ignored
                    for x in range(len(temp_sequence)):
                        if temp_sequence[x] != perfect_target[x]:
                            mismatch_list[x][4] += 1  #4th position is number of deletions in "mismatch list"
                            indel_sequences += 1
                            break


                all_matching = True
                complement_base_pos_dict = {'T': 1, 'A': 0, 'G': 3, 'C': 2}

                if len(temp_sequence) == len(perfect_target):
                    mm_count_temp = 0
                    for x in range(len(temp_sequence)):  #counting total mismatches
                        if temp_sequence[x] != perfect_target[x]:
                            mm_count_temp += 1

                    for x in range(len(temp_sequence)):
                        if temp_sequence[x] != perfect_target[x]:
                            if mm_count_temp == 1: #single mismatches
                                mismatch_list[x][complement_base_pos_dict[temp_sequence[x]]] += 1
                                single_mm_sequences += 1
                            if mm_count_temp > 1: #multiple mismatches
                                mismatch_list[x][5] += 1
                                multi_mm_sequences += 1
                            break
                    if mm_count_temp == 0 :
                        wild_type_sequences += 1


                # making seq count list - to compare to control file to see if mutated sequence was there already
                in_list = False
                for x in range(len(seq_count_list)):
                    if temp_sequence in seq_count_list[x][0]:
                        seq_count_list[x][1] += 1
                        in_list = True
                        break
                if in_list == False:
                    seq_count_list.append([temp_sequence,0,0,0,0,0,0])  #sequence, count, control count, status, base before, base after, position

    print(sequence_lines, 'sequence lines')
    seq_count_list = [item for item in seq_count_list if item[1] > 0 and item[0] != perfect_target]  # WT Cas12a protospacer


    # match sequences in list to control file
    # control_file = open('lys-J1-Fn-563-1_S228_L001_R1_001.fastq')
    # for linenumber, line in enumerate(control_file):
    #     if line.find(to_find1_R2) != -1 and line.find(to_find2_R2) != -1:
    #
    #         temp_sequence = line[line.find(to_find1_R2) + len(to_find1): line.find(to_find2_R2)]
    #         for sequence2 in seq_count_list:
    #             if sequence2[0] == temp_sequence:
    #                 sequence2[2] += 1

    # classifying sequences by mismatch or deletion, then mismatch type (M to N)
    for sequence4 in seq_count_list:
        for x in range(len(perfect_target)):
            if sequence4[0][x] != perfect_target[x]:
                sequence4[6] = x + 1  #sixth position is mismatch/indel location
                sequence4[4] = perfect_target[x]
                sequence4[5] = sequence4[0][x]
                break
        if len(sequence4[0]) < len(perfect_target):
            sequence4[3] = 'deletion'
        if len(sequence4[0]) > len(perfect_target):
            sequence4[3] = 'extension'
        if len(sequence4[0]) == len(perfect_target):
            sequence4[3] = 'mutation'


    print(seq_count_list)


    for x in range(len(mismatch_list)): #making a fraction of mutated sequences column for heat maps
        sum = 0
        for y in range(len(mismatch_list[x])):
            sum += mismatch_list[x][y]

        mismatch_list[x].append(sum/sequence_lines)
    print(mismatch_list)

    if sequence_lines != 0 :  # shouldn't actually need this line if you glob files correctly above

        file_name = file_pair[0][0:22]
        print(file_name)
        worksheet = workbook.add_worksheet(file_name)

        # Creating chart object
        chart = workbook.add_chart({'type': 'column', 'subtype': 'stacked'})
        # worksheet.write_column('A1', mismatch_list)

        row = 2
        column = 2
        for A0, T1, C2, G3, Deletion4, multi_mm_5, fractions in mismatch_list:
            worksheet.write(row, column, A0)
            worksheet.write(row, column + 1, T1)
            worksheet.write(row, column + 2, C2)
            worksheet.write(row, column + 3, G3)
            worksheet.write(row, column + 4, Deletion4)
            worksheet.write(row, column + 5, multi_mm_5)
            worksheet.write(row, column + 6, fractions)
            row += 1

        worksheet.write_column('B1', perfect_target_list)
        name_list = ['A', 'T', 'C', 'G', 'Deletion', 'multi_mm']
        for type in range(len(name_list)):
            chart.add_series({'values': [file_name, 2, 2 + type, 2 + len(perfect_target) - 1, 2 + type],
                              # [sheetname, first_row, first_col, last_row, last_col]
                              'categories': [file_name, 0, 1, 42, 1],
                              'name': name_list[type]
                              })  # values strings specifies what values in the excel sheet to make the chart from

        chart.set_x_axis({'name': 'mutation position'})
        chart.set_y_axis({'name': 'sequence counts'})

        worksheet.insert_chart('N3', chart)

        worksheet.write(49, 2, 'mutated sequence')
        worksheet.write(49, 3, 'count')
        worksheet.write(49, 4, 'count in control ')

        if seq_count_list != []:  # writing counts of sequences found in control file
            row = 50
            column = 2
            for sequence1, count, control_count, status, base_before, base_after, position in seq_count_list:
                worksheet.write(row, column, sequence1)
                worksheet.write(row, column + 1, count)
                worksheet.write(row, column + 2, control_count)
                worksheet.write(row, column + 3, status)
                worksheet.write(row, column + 4, position)
                worksheet.write(row, column + 5, base_before)
                worksheet.write(row, column + 6, 'to')
                worksheet.write(row, column + 7, base_after)
                row += 1

        # writing checks to make sure all sequences are accounted for
        worksheet.write(2, 9, str(sequence_lines) )
        worksheet.write(2, 10, str("total sequences") )
        worksheet.write(3, 9, str(single_mm_sequences + indel_sequences) + " single mutation or indel sequences")
        worksheet.write(4, 9, str(wild_type_sequences) + " wild-type sequences")
        worksheet.write(5, 9, str(sequence_lines - (indel_sequences + single_mm_sequences + wild_type_sequences + multi_mm_sequences)) + " other sequences")
        worksheet.write(6, 9, "fraction mutated is " )
        worksheet.write(6, 10, str(round(((indel_sequences + single_mm_sequences + multi_mm_sequences) / sequence_lines), 2)))

        # values_string = "=" + file_name + "!$A$1:$A$24"
        # print(values_string)
        options = {'fill': {'none': True}, 'line': {'none': True}, 'font': {'color': 'black',
                                                                            'size': 14},
                   'font': {'underline': True}}
        worksheet.insert_textbox('R5', 'Mutated Base', options)

workbook.close()

