program 1
#########################################################
### reading MultiFASTA flile containing DNA sequences ###
#########################################################
# opening the files containing sequences in FASTA format to be analyesd
filename = "example file"
file_contents = open(filename, 'r')
# two functions used to split the headers and the sequences into seperate unordered lists
(headers[] and sequences[])
# the header seqeunces is identified by the start character '>' and the lines of sequence
following that are joined and appended to the list sequences[]
# meanwhile the header is appended to the list headers[]
def FastaFileSequences_into_List(filename):
# opening file and defining the empty lists to be filled
 file_contents = open(filename, 'r')
 sequences = []
 seqFragments = []
# identifying new entries in the multifasta file by the header sequence that starts with the
'>' character
# if a line is identified as a header, the seqFragments list is kept emtpy
 for line in file_contents:
 if line.startswith('>'):
 if seqFragments:
 sequence = ''.join(seqFragments)
 sequences.append(sequence)
 seqFragments = []
# each line of sequences in the FASTA format contains an end of line character at the end and
they are hence interpreted as seperate lines
# sequence fragments are joined together to create a single sequence
 else:
 seq = line.rstrip()
 seqFragments.append(seq)
 if seqFragments:
# lines of nucleotide sequences are identified and appended to a list
# gaps are eradicated to link each sequence
 sequence = ''.join(seqFragments)
 sequences.append(sequence)
 file_contents.close()
 return sequences
def FastaFileHeaders_into_List(filename):
 file_contents = open(filename, 'r')
 headers = []
#if the firts character of the line is '>', the line is added to the headers list
 for line in file_contents:
 if line[0] == ">":
 headers.append(line)
 return headers
###headers and sequences are defined and called outside of their resepective functions
headers = FastaFileHeaders_into_List(filename)
sequences = FastaFileSequences_into_List(filename)
number_of_sequences = len(sequences)
## Writing files containing retrieved information
with open("SEQUENCES", "a") as output:
 output.write(str(sequences))
with open("HEADERS", "a") as output:
 output.write(str(headers))
with open("NUMBER OF SEQUENCES", "a") as output:
 output.write(str(number_of_sequences))





Program 2:
## Importing the regular expression and numpy (for calculations) modules
import re
import numpy
sequences = [ENTER SEQUENCE HERE]
########################################
# identifying repeats of any unit size #
########################################
# establishing the empty lists that will store data sets about the nature of the repeats
unit = []
start_of_repeat = []
end_of_repeat = []
length = []
length1 = []
sequence_length = []
sequence_number_containing_repeat = []
sequence_containing_repeat = []
repeat_info = []
full_repeat = []
relative_position_of_start = []
relative_position_of_end = []
middle_of_repeat = []
relative_position_of_middle = []
sequences_containing_repeats = 0
total_repeats_found = 0
## reading each sequence individually
for seq in sequences:
## printing which sequence number the programm is working on
 print("Working on sequence number " + str(sequences.index(seq)) + ":")
## if a repeat is present in a sequence, the whole sequence is searched
## the regular expression contains two variables that determine unit size and echo:
## ([atcg]{a,b})(\1){x,y}(?!\1)?
## a and b define the upper and lower limits of unit size
## x and y define the upper and lower limits of echo
## a,b,x and y may be left blank if either upper or lower limits are not required
## the variables used here search for repeats of unit size 1+ repeated 3+ times
## the repeat variable refers to the times the unit is repeated additionally to the one
identified
## therefore, the echo will be x+1 (in this case 4)
 if re.search(r"([atcg]{1,})(\1){3,}(?!\1)?", seq):
 repeat = re.finditer(r"([atcg]{1,})(\1){3,}(?!\1)?", seq)
## notification regarding whether or not a repeat was found
 print("Repeat found!")
## addition of one to the tally
 sequences_containing_repeats += 1
## collecting data
 for r in repeat:
 total_repeats_found += 1
 sequence_number_containing_repeat.append(sequences.index(seq))
 sequence_containing_repeat.append(seq)
 repeat_info.append(r)
 start_of_repeat.append(r.start())
 end_of_repeat.append(r.end())
 full_repeat.append(r.group(0))
 unit.append(r.group(1))
 sequence_length.append(len(seq))
## notification regarding whether or not a repeat was found
 else:
 print("No repeat found!")
for i in full_repeat:
 length.append(len(i))
##print(length)
## calculation using start and end position of the repeat within the sequence to find the
position of the middle
array = (numpy.add(start_of_repeat, end_of_repeat))/2
for i in array:
 middle_of_repeat.append(i)
33
## using the length of the respective sequence the repeat was found in, the relative position
of the repeat is calculated
position_m = (numpy.divide(middle_of_repeat, sequence_length))*100
for i in position_m:
 relative_position_of_middle.append(i)
with open('POSITION OF MIDDLE OF REPEAT', 'a') as output:
 output.write(str(relative_position_of_middle))
print(relative_position_of_middle)
#################################
#### COLLECTING ALL LENGTHS #####
#################################
## all repeats of all lengths are seperated into groups of total repeats of each length
def all_lengths(l):
 lengths = 0
 for i in length:
 if i == l:
 lengths +=1
 return lengths
l = []
for i in range(4,21,1):
 l.append(all_lengths(i))
print(l)
## for easier interpretation of the list created above, it is zipped with the numbers 4-20
r = list(range(4,21,1))
both = list(zip(r,l))
print(both)
#########################
##### UNIT FINDER #######
#########################
print('yes')
r_length = []
full_repeat = []
unit = []
sequence_length = []
repeat_length = []
### in order to find repeats of larger unit sizes but shorter length, I need to change the
margins of the second variable ####
## the function 'unit_finder' uses two variables in order to identify repeats based on their
unit size and length
## the regular expression identifies a homomeric repeat for different unit sizes depending on
what is searched for
## e.g. a repeat of 12 'A's can be interpreted as a repeat of 1x12, 2x6 or 3x4
## to avoid that the program returns false repeats, the list of units is analysed for
repeating nucleotides and removes those from the 'units' list
def unit_finder(unit_size, length):
## list and values are established
 full_repeat = []
 units = []
 false_unit = 0
 all_units = 0
## print("Unit size: " + str(unit_size))
 for seq in sequences:
## print("Working on sequence number " + str(sequences.index(seq)) + ":")
## the regular expression '([atcg]{a,b})(\1){x,y}(?!\1)?' is utilised and the two variables
are integrated
 pattern = r"([atcg]{" + re.escape(str(unit_size)) + r"})(\1){" +
re.escape(str(length)) + r",}(?!\1)?"
## pattern = r"([atcg]{" + re.escape(str(unit_size)) + r"})(\1){3,}(?!\1)?"
## searching for repeats in each sequence
## first it is checked if a repeat is present
 if re.search(pattern, seq):
## if a repeat is present, the whole sequence is searched for more repeats
 repeat = re.finditer(pattern, seq)
 for r in repeat:
34
## data collection
 unit = r.group(1)
 units.append(unit)
 all_units += 1
 full_repeat.append(r.group(0))
## the following two for loops identify homomeric repeats and remove them from the two lists
## these need to be commented out when searching for homomeric repeats/ repeats of unit size
1
 for homo in full_repeat:
## a second regular expression is used to identify those homomeric repeats
 pattern_2 = r"([atgc])\1+"
 if re.search(pattern_2, homo):
 false_unit +=1
 while homo in full_repeat:
 full_repeat.remove(homo)
## for homo in units:
## pattern_2 = r"([atgc])\1+"
## if re.search(pattern_2, homo):
## while homo in units:
## units.remove(homo)
## calculation for the number of repeats containing 'true' repeats, not homomeric repeats
## using the values collected in the programme
 total_sequences_containing_repeats = all_units - false_unit
 return units, total_sequences_containing_repeats, full_repeat
## the following commands ustilise the 'unit_'finder' function for unit sizes 1-8 and lengths
4+
## the function is called individually - while searching for one unit size, the other
commands need to be commented out
## identifying repeats of unit size 1
##print('Results for Unit Size 1')
##unit_finder_1 = unit_finder(1,3)
##repeat_len_1 = []
##repeat_1 = unit_finder_1[2]
##for i in repeat_1:
## repeat_len_1.append(len(i))
####print('Length 1')
####print(repeat_len_1)
##def all_lengths(l):
## lengths = 0
## for i in repeat_len_1:
## if i == l:
## lengths +=1
## return lengths
##l = []
##for i in range(4,21,1):
## l.append(all_lengths(i))
##print(l)
##total_sequences_containing_repeats_1 = unit_finder_1[1]
##print('TotalSeq 1')
##print(total_sequences_containing_repeats_1)
####units_found_1 = unit_finder_1[0]
####print(units_found_3)
####print(len(units_found_3))
##print('Total 1')
##print(len(repeat_len_1))
## identifying repeats of unit size 2
print('Results for Unit Size 2')
unit_finder_2 = unit_finder(2,1)
repeat_len_2 = []
repeat_2 = unit_finder_2[2]
for i in repeat_2:
 repeat_len_2.append(len(i))
##print('Length 2')
##print(repeat_len_2)
def all_lengths(l):
 lengths = 0
 for i in repeat_len_2:
 if i == l:
 lengths +=1
 return lengths
l = []
35
for i in range(4,21,1):
 l.append(all_lengths(i))
print(l)
total_sequences_containing_repeats_2 = unit_finder_2[1]
print('TotalSeq 2')
print(total_sequences_containing_repeats_2)
##units_found_2 = unit_finder_2[0]
##print(units_found_3)
##print(len(units_found_3))
print('Total 2')
print(len(repeat_len_2))
## identifying repeats of unit size 3
print('Results for Unit Size 3')
unit_finder_3 = unit_finder(3,1)
repeat_len_3 = []
repeat_3 = unit_finder_3[2]
##print(repeat_3)
for i in repeat_3:
 repeat_len_3.append(len(i))
##print('Length 3')
##print(repeat_len_3)
def all_lengths(l):
 lengths = 0
 for i in repeat_len_3:
 if i == l:
 lengths +=1
 return lengths
l = []
for i in range(4,21,1):
 l.append(all_lengths(i))
print(l)
total_sequences_containing_repeats_3 = unit_finder_3[1]
print('TotalSeq 3')
print(total_sequences_containing_repeats_3)
##units_found_3 = unit_finder_3[0]
##print(units_found_3)
##print(len(units_found_3))
print('Total 3')
print(len(repeat_len_3))
print('Results for Unit Size 4')
unit_finder_4 = unit_finder(4,1)
repeat_len_4 = []
repeat_4 = unit_finder_4[2]
##print(repeat_4)
for i in repeat_4:
 repeat_len_4.append(len(i))
##print('Length 4')
##print(repeat_len_4)
def all_lengths(l):
 lengths = 0
 for i in repeat_len_4:
 if i == l:
 lengths +=1
 return lengths
l = []
for i in range(4,21,1):
 l.append(all_lengths(i))
print(l)
total_sequences_containing_repeats_4 = unit_finder_4[1]
##print('TotalSeq 4')
##print(total_sequences_containing_repeats_4)
##units_found_4 = unit_finder_4[0]
##print(units_found_4)
##print(len(units_found_4))
print('Total 4')
print(len(repeat_len_4))
print('Results for Unit Size 5')
unit_finder_5 = unit_finder(5,1)
repeat_len_5 = []
36
repeat_5 = unit_finder_5[2]
##print(repeat_5)
for i in repeat_5:
 repeat_len_5.append(len(i))
##print('Length 5')
##print(repeat_len_5)
def all_lengths(l):
 lengths = 0
 for i in repeat_len_5:
 if i == l:
 lengths +=1
 return lengths
l = []
for i in range(4,21,1):
 l.append(all_lengths(i))
print(l)
total_sequences_containing_repeats_5 = unit_finder_5[1]
##print('TotalSeq 5')
##print(total_sequences_containing_repeats_5)
##units_found_5 = unit_finder_5[0]
##print(units_found_5)
##print(len(units_found_5))
print('Total 5')
print(len(repeat_len_5))
print('Results for Unit Size 6')
unit_finder_6 = unit_finder(6,1)
repeat_len_6 = []
repeat_6 = unit_finder_6[2]
##print(repeat_6)
for i in repeat_6:
 repeat_len_6.append(len(i))
##print('Length 6')
##print(repeat_len_6)
def all_lengths(l):
 lengths = 0
 for i in repeat_len_6:
 if i == l:
 lengths +=1
 return lengths
l = []
for i in range(4,21,1):
 l.append(all_lengths(i))
print(l)
total_sequences_containing_repeats_6 = unit_finder_6[1]
##print('TotalSeq 6')
##print(total_sequences_containing_repeats_6)
#units_found_6 = unit_finder_6[0]
#print(units_found_6)
##print(len(units_found_6))
print('Total 6')
print(len(repeat_len_6))
print('Results for Unit Size 7')
unit_finder_7 = unit_finder(7,1)
repeat_len_7 = []
repeat_7 = unit_finder_7[2]
##print(repeat_7)
for i in repeat_7:
 repeat_len_7.append(len(i))
##print('Length 7')
##print(repeat_len_7)
def all_lengths(l):
 lengths = 0
 for i in repeat_len_7:
 if i == l:
 lengths +=1
 return lengths
l = []
for i in range(4,21,1):
 l.append(all_lengths(i))
print(l)
total_sequences_containing_repeats_7 = unit_finder_7[1]
37
##print('TotalSeq 7')
##print(total_sequences_containing_repeats_7)
#units_found_7 = unit_finder_7[0]
#print(units_found_7)
##print(len(units_found_7))
print('Total 7')
print(len(repeat_len_7))
print('Results for Unit Size 8')
unit_finder_8 = unit_finder(8,1)
repeat_len_8 = []
repeat_8 = unit_finder_8[2]
##print(repeat_8)
for i in repeat_8:
 repeat_len_8.append(len(i))
##print('Length 8')
##print(repeat_len_8)
def all_lengths(l):
 lengths = 0
 for i in repeat_len_8:
 if i == l:
 lengths +=1
 return lengths
l = []
for i in range(4,21,1):
 l.append(all_lengths(i))
print(l)
total_sequences_containing_repeats_8 = unit_finder_8[1]
##print('TotalSeq 8')
##print(total_sequences_containing_repeats_8)
#units_found_8 = unit_finder_7[0]
#print(units_found_8)
##print(len(units_found_8))
print('Total 8')
print(len(repeat_len_8))