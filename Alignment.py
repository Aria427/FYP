"""

"""
from Bio import SeqIO

#To efficiently read the genome:
def readGenome(filename):
    genome = ''   
    #rU - open file for reading in universal readline mode (works across platforms due to differing newline characters)
    filehandle = open(filename, 'rU', encoding="latin1") 
    for record in SeqIO.parse(filehandle, "fasta"): #SeqRecord iterator
        genome += record.seq #parses records one by one, without changing the file order
    filehandle.close()
    return genome

#To inefficiently read the genome:
def readGenomeInefficient(filename):
    genome = '' 
    with open(filename, 'r', encoding="latin1") as file: #opening a file for reading
        for line in file:
            if line[0] != '>': #ignore header line with genome information
                genome += line.rstrip() #add each line of bases to the string 
                #rstrip() removes any trailing whitespace from the ends of the string (trim off new line/tab/space)
    return genome
  
#To read and parse the reads:
def readSequence(filename):
    sequences = []
    with open(filename, encoding="latin1") as file:  
        while True: #loops every 4 lines (each read is a set of 4) indefinitely until EOF
            file.readline() #read tag
            seq = file.readline().rstrip() #string of DNA bases
            file.readline() #+
            file.readline() #quality sequence - each quality score is Phredd33 encoded (converts to ASCII character)
            if len(seq) == 0: #reached EOF
                break
            sequences.append(seq)
    return sequences

#The following is a naive algorithm for exact matching where all occurrences are recorded:
def naiveExact(pattern, text):
    matchOffsets = []
    #loop through every position from where P could start without running past the end of T 
    for i in range(len(text) - len(pattern) + 1): #loop over all possible alignments of P in T from left to right	 
        match = True
        for j in range(len(pattern)): #loop over characters in P from left to right
            #i'th alignment and j'th character            
            if text[i+j] != pattern[j]: #compare characters of P with T
                match = False #mismatch; reject alignment
                break
        if match:
            matchOffsets.append(i)  #all characters matched; record
    return matchOffsets

#Hamming distance = minimum no of substitutions required to change one string into another.
#The following is also a naive algorithm but for approximate matching using the Hamming distance:
def naiveApproxHamming(pattern, text, maxHammingDist=1):
    matchOffsets = []
    for i in range(len(text) - len(pattern) + 1): 
        mismatches = 0
        for j in range(len(pattern)):
            if text[i+j] != pattern[j]:
                mismatches += 1     #mismatch
                if mismatches > maxHammingDist:
                    break           #exceeded maximum distance
        if mismatches <= maxHammingDist: #approximate match
            matchOffsets.append(i)
    return matchOffsets      

#The genome is double stranded and so the reads can come from one strand or the other.    
#To match both the read and the reverse complement of the read to the genome: 
def reverseComplement(read):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'} #each base is associated with its complementary base
    reverseRead = ''
    for base in read:
        reverseRead = complement[base] + reverseRead #complement added to beginning in order to reverse the read from end to start
    return reverseRead

#To align the reads against the genome to see how many match:
def align(reads, genome):
    readsMatched = 0
    readsCount = 0
    readsOffsets = []
    for read in reads:
        read = read[:25] #prefix of read as all 100 bases have a smaller chance of matching
        matches = naiveApproxHamming(read, genome) #check if read matches in forward direction of genome
        matches.extend(naiveApproxHamming(reverseComplement(read), genome)) #add results of any matches in reverse complement of genome
        readsCount += 1
        if len(matches) > 0: #match - read aligned in at least one place
            readsMatched += 1
        readsOffsets.append(matches)
    return readsMatched, readsCount, readsOffsets
    