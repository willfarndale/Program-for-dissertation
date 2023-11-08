import re
import numpy

# Function to read FASTA file contents into lists
def FASTAFileSequences_into_List(file_contents):
    sequences = []
    seqFragments = []
    for line in file_contents:
        if line.startswith('>'):
            if seqFragments:
                sequence = ''.join(seqFragments)
                sequences.append(sequence)
            seqFragments = []
        else:
            seq = line.strip()
            seqFragments.append(seq)
    if seqFragments:
        sequence = ''.join(seqFragments)
        sequences.append(sequence)
    return sequences

def unit_finder(unit_size, length, sequences):
    full_repeat = []
    units = []
    false_unit = 0
    all_units = 0

    for seq in sequences:
        pattern = r"([atcg]{" + re.escape(str(unit_size)) + r"})(\1){" + re.escape(str(length)) + r",}(?!\1)?"

        if re.search(pattern, seq):
            repeat = re.finditer(pattern, seq)
            for r in repeat:
                unit = r.group(1)
                units.append(unit)
                all_units += 1
                full_repeat.append(r.group(0))

                for homo in full_repeat:
                    pattern_2 = r"([atgc])\1+"
                    if re.search(pattern_2, homo):
                        false_unit += 1
                        while homo in full_repeat:
                            full_repeat.remove(homo)

    total_sequences_containing_repeats = all_units - false_unit
    return units, total_sequences_containing_repeats, full_repeat

filename = "C:\\Users\\willf\\OneDrive\\Desktop\\Sequences.txt.txt"

with open(filename, 'r') as file_contents:
    headers = FASTAFileHeaders_into_List(file_contents)
    file_contents.seek(0)
    sequences = FASTAFileSequences_into_List(file_contents)

numbers_of_sequences = len(sequences)

# Your unit_finder calls for different unit sizes
print('Results for Unit Size 2')
unit_finder_2 = unit_finder(2, 1, sequences)
# Add more unit_finder calls as needed

# Rest of your code here

# Your previous code for calling unit_finder
# Example for unit size 2 and length 1
print('Results for Unit Size 2')
unit_finder_2 = unit_finder(2, 1, sequences)
print('Total Sequences Containing Repeats (Unit Size 2):', unit_finder_2[1])
print('Unique Units Found (Unit Size 2):', unit_finder_2[0])
print('Total Repeats Found (Unit Size 2):', len(unit_finder_2[2]))
print('Lengths of Repeats (Unit Size 2):', [len(repeat) for repeat in unit_finder_2[2]])

# Add similar print statements for other unit sizes and lengths as needed
