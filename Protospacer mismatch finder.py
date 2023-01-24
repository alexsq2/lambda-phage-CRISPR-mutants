import xlsxwriter
import glob

def reverse_comp(in_sequence): #faster way lifted from stackoverflow
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    seq = in_sequence
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return  reverse_complement


workbook = xlsxwriter.Workbook('(top) general Cas12a escapers.xlsx')  # naming the output excel sheet

# selecting files to analyze, here it's all the Cas12a ones and the control
total_file_list = []
total_file_list += glob.glob('*Fn_J_red_441*')


print(total_file_list)
print(len(total_file_list) , 'files gathered')

file_list_length_half = int(len(total_file_list) / 2)
print(file_list_length_half)

total_file_list_pairs = [] #this list will contain pairs of files ordered by their order in total file list
#total file list must be in order or the wrong files will be matched!!

for x in range(file_list_length_half):
    total_file_list_pairs.append([total_file_list[2*x]] + [total_file_list[2*x + 1]])

print(total_file_list_pairs)

all_sample_info = [] # this will be written on the last sheet containing info for all of the samples processed
worksheet_bar_graphs = workbook.add_worksheet('Bar graphs')
worksheet_all = workbook.add_worksheet('All info')
chart_coord1 = 1
chart_coord2 = 1
chart_rep = 0
chart_rep_9 = 0
replicate_File_count = 0 # 3 pairs of files will be analyzed at once


perfect_target = 'TTTGTCATTTATGGAGCGTGAGGA' #FORWARD orientation

to_find1_R1 = 'TGCGGGCGGT'  #extract the target sequence between these two sequences from the R1 file
to_find2_R1 = 'ATGGGT'

to_find1_R2 = 'TTTACCCAT'  #extract the target sequence between these two sequences from the R2 file
to_find2_R2 = 'ACCGCCCG'


def file_pair_processing(file_pair,perfect_target,to_find1_R1,to_find2_R1,to_find1_R2,to_find2_R2):
    r1_file = open(file_pair[0])
    r2_file = open(file_pair[1])
    line_count_r1 = 0

    end_of_lines = False
    r1_temp_seq_dict = {}
    r2_temp_seq_dict = {}
    while not end_of_lines:  # lines_R1 and R2 SHOULD have the same length
        temp_sequence_r1 = ''
        temp_sequence_r2 = ''
        r1_line = r1_file.readline()
        r2_line = r2_file.readline()

        if r1_line.find(to_find1_R1) != -1 and r1_line.find(to_find2_R1) != -1 and r2_line.find(
                to_find1_R2) != -1 and r2_line.find(to_find2_R2) != -1:


            temp_sequence_r1 = r1_line[r1_line.find(to_find1_R1) + len(to_find1_R1): r1_line.find(
                to_find2_R1)]  # 72 here because there are many occcurences of "to_find_R1" in the line
            temp_sequence_r2 = r2_line[r2_line.find(to_find1_R2) + len(to_find1_R2): r2_line.find(to_find2_R2)]
            # print(temp_sequence_R1)
            # print(temp_sequence_R2)
            line_count_r1 += 1
            if temp_sequence_r1 not in r1_temp_seq_dict:
                r1_temp_seq_dict[temp_sequence_r1] = {'count': 1}
            else:
                r1_temp_seq_dict[temp_sequence_r1]['count'] += 1

            if temp_sequence_r2 not in r2_temp_seq_dict:
                r2_temp_seq_dict[temp_sequence_r2] = {'count': 1}
            else:
                r2_temp_seq_dict[temp_sequence_r2]['count'] += 1

        if not r1_line:
            end_of_lines = True
            r1_file.close()
            r2_file.close()
            print(r1_temp_seq_dict)
            print(r2_temp_seq_dict)
    temp_seq_dict = {}
    for item1 in r1_temp_seq_dict:
        if reverse_comp(item1) in r2_temp_seq_dict:
            temp_seq_dict[item1] = {
                'count': min(r1_temp_seq_dict[item1]['count'], r2_temp_seq_dict[reverse_comp(item1)]['count'])}


    valid_sequences = 0
    for sequence in temp_seq_dict:
        valid_sequences += temp_seq_dict[sequence]['count']
    for sequence in temp_seq_dict:
        mm_count_temp = 0
        mm_pos_1 = 'NA'
        mm_pos_2 = 'NA'

        for position in range(min(len(sequence),len(perfect_target))):  # counting total mismatches
            if sequence[position] != perfect_target[position]:
                mm_count_temp += 1
                if mm_count_temp == 1:
                    mm_pos_1 = position
                if mm_count_temp == 2:
                    mm_pos_2 = position
        temp_seq_dict[sequence]['deletion'] = False
        temp_seq_dict[sequence]['mismatches'] = mm_count_temp
        temp_seq_dict[sequence]['mismatch1'] = mm_pos_1
        temp_seq_dict[sequence]['mismatch2'] = mm_pos_2
        temp_seq_dict[sequence]['fraction'] = temp_seq_dict[sequence]['count']/valid_sequences
        if len(sequence) == len(perfect_target) - 1:
            temp_seq_dict[sequence]['deletion'] = True


    print(temp_seq_dict)
    print(len(temp_seq_dict))
    return temp_seq_dict

perfect_target_list = []
for x in range(len(perfect_target)):
    perfect_target_list += perfect_target[x]

for file_pair in total_file_list_pairs:
    file_dict = file_pair_processing(file_pair,perfect_target,to_find1_R1,to_find2_R1,to_find1_R2,to_find2_R2)
    print(file_dict)
    info_list_one_sample = [] # contains the info for only one sample: sample name, % mutated, number of samples over 1% of all reads

    file_name = file_pair[0][0:22]
    worksheet = workbook.add_worksheet(file_name)
    mismatch_pos_dict = {'A':0, 'T':1, 'C':2, 'G':3, 'deletion':4, 'multi':5}
    row = 2
    column = 2
    multi_mismatch_array = []
    for x in range(len(perfect_target)):
        multi_mismatch_array.append(0)


    for sequence in file_dict:
        col_write = 'unassigned'
        if file_dict[sequence]['mismatches'] == 1 and file_dict[sequence]['deletion'] == False :
            col_write = mismatch_pos_dict[sequence[  file_dict[sequence]['mismatch1']       ]]

        if file_dict[sequence]['mismatches'] > 1 and file_dict[sequence]['deletion'] == True:
            col_write = 4
        if file_dict[sequence]['mismatches'] > 1 and file_dict[sequence]['deletion'] == False:
            multi_mismatch_array[file_dict[sequence]['mismatch1']] += file_dict[sequence]['fraction']


        if col_write != 'unassigned':
            worksheet.write(file_dict[sequence]['mismatch1'] + row, col_write + column, file_dict[sequence]['fraction'])


    chart = workbook.add_chart({'type': 'column', 'subtype': 'stacked'})

    worksheet.write_column('B3', perfect_target_list)
    name_list = ['A', 'T', 'C', 'G', 'Deletion', 'multi_mm']
    for type in range(len(name_list)):
        chart.add_series({'values': [file_name, 2, 2 + type, 2 + len(perfect_target) - 1, 2 + type],
                          # [sheetname, first_row, first_col, last_row, last_col]
                          'categories': [file_name, 0, 1, 23, 1],
                          'name': name_list[type]
                          })  # values strings specifies what values in the excel sheet to make the chart from

    chart.set_x_axis({'name': 'mutation position'})
    chart.set_y_axis({'name': 'sequence counts'})
    chart.set_title({'name': file_name})
    chart.set_size({'width': 800, 'height': 300})

    worksheet.insert_chart('N3', chart)  # uncomment if you want graphs on separate sheets

workbook.close()
