# Python for Genomic Data Science

## About
As part of the Python for Genomic Data Science specialization course, I was tasked with investigating a multi-FASTA formatted file, which contains DNA sequences. In my analysis, I will be using Python as my programming language.

I will answer the following questions in my analysis:
- How many records are in the file?
- What are the lengths of the sequences in the file? What is the longest sequence and what is the shortest sequence? Is there more than one longest or shortest sequence? What are their identifiers?
- Identify all ORFs present in each sequence of the FASTA file. What is the length of the longest ORF in the file? What is the identifier of the sequence containing the longest ORF? For a given sequence identifier, what is the longest ORF contained in the sequence represented by that identifier? What is the starting position of the longest ORF in the sequence that contains it?

## Analysis
### Importing File
First I will open the FASTA file, which is located on my computer's desktop. If the file is not found the system will return my specified error message.

```
try:
	f = open("/Users/khaitlinbernaldez/Desktop/dna.example.fasta")
except IOError:
	print("File dna.example.fasta does not exist!")
```

Each sequence in the FASTA file begins with a greater-than symbol (>) and a single-line description, followed by the DNA sequence. Knowing this, I will create a dicitionary called seqs, where the greater-than symbol in the file acts as the delimeter for each key in the dictionary.

```
seqs = {}
for line in f:
	line = line.rstrip()
	if line[0] == '>':
		words = line.split()
		name = words[0][1:]
		seqs[name] = ''
	else:
		seqs[name] = seqs[name] + line
```

After reading the FASTA file, I can now determine how many records are in the file using the len command.

```
print(len(seqs))
```

### DNA Sequences

Now, I would like to compare the lengths of each DNA sequence. First, I will be creating a dictionary of the DNA sequences, where the key is the identifier and the value is the DNA sequence. This dictionary will be organized by length.

```
sorted_values = sorted(seqs.values(), key = len)
sorted_seq = {}

for i in sorted_values:
		for k in seqs.keys():
				if seqs[k] == i:
					sorted_seq[k] = seqs[k]
					break
```

We will now generate an output of the identifier of the DNA sequences, followed by the length of the corresponding sequence. Since we have already sorted our sequences by length, we can determine the longest and shortest sequences and their identifiers from this output.

```
for name, seq in sorted_seq.items():
	 print('NAME:' + name, 'LENGTH:' + str(len(seq)))
```

### Open Reading Frames (ORFs)
ORFs can be found by identifying the locations of the start and stop codons in a sequence, and including all bases between and inlcuding these locations.

```
def find_orf(sequence, frame):
	start_indices = []
	stop_indices = []
	orf = []
	front = 0
	start_codon = 'ATG'
	stop_codon = ['TAA', 'TAG', 'TGA']
	for i in range(frame, len(sequence), 3):
		codon = sequence[i:i+3]
		if codon == start_codon:
			start_indices.append(i)
	## print("start:" + str(start_indices))
	for i in range(frame, len(sequence), 3):
		codon = sequence[i:i+3]
		if codon in stop_codon:
			stop_indices.append(i)
	## print ("stop:" + str(stop_indices))
	for i in range(0, len(start_indices)):
		for j in range(0, len(stop_indices)):
			if start_indices[i] < stop_indices[j] and start_indices[i] > front:
				orf.append(sequence[start_indices[i]:stop_indices[j] + 3])
				front = stop_indices[j] + 3
	return orf
```

From the idenitified ORFs, we can find the ORF with the longest length using the following lines of code:

```
def find_max_length(sequence, frame):
	orf_length = []
	for i in range(0, len(find_orf(sequence, frame))):
		orf_length.append(len(find_orf(sequence, frame)[i]))
	for i in orf_length:
		if i:
			return max(orf_length)
```

#### Identifier
If we wanted to find the longest ORF with a specific identifier, we would repeat the previous code with the sequence parameter equal to the DNA sequence associated with our target identifier.


#### Reading Frame
To find the length of the longest ORF in a specific reading frame, we use the following code:

```
max_length = []
frame = ** # assign target reading frame value (0, 1, 2)
for value in seqs.values():
	if find_max_length(value, frame) != None: # Removes nulls
		max_length.append(find_max_length(value, frame))
longest_length = max(max(length))
print(longest_length)
```

#### Starting Position
Now that we have identified the longest ORFs in different scenarios, it might be useful to know the location of the ORF. This can be done with the following code:

```
frame = ** # assign target reading frame value (0, 1, 2)
target_sequence = []
for value in seqs.values():
	for i in range(0, len(find_orf(value, frame))):
		if len(find_orf(value, frame)[i]) == longest_length:
			target_sequence.append(find_orf(value, frame)[i])

for value in seqs.values():
	for i in range(0, len(value)):
		seq = value[i:i+len(target_sequence)]
		if seq == target_sequence:
			print(i)
```
