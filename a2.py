f = open('nc_045512.txt')
comment = f.readline()
print(comment)
seq = f.read()
seq = seq.replace('\n', '')
print('------------------------------------------------')
print('Length of SARS-CoV-2 viral sequence is', len(seq), 'nucleotides')
print('------------------------------------------------')
def firstORF(sequence, start_position) :
    n = len(sequence)
    # count will be the number of codons from ATG to a stop,
    # but not counting the stop
    count = 0
    index = ()
    for i in range(start_position, n, 3) :
        triplet = sequence[i:i + 3]
        # Check for a start codon
        if triplet == 'ATG' :
            index += (i,)
        # Look at triplets up to a stop codon
            for j in range(i, n, 3) :
                triplet = sequence[j:j + 3]
                if triplet in ('TAG', 'TGA', 'TAA') :
                    index += (j+3,)
                    return index
                else :
                    count += 1
            # Reaching this point, no stop codon was found, so return 0
            index += (None,)
            return index
    # Reaching this point, no start codon was found, so return 0
    index += (None,)
    return index
#seq = 'ATG222333444555TAG6'
for frame in range(1, 4) :
    firstorflength = firstORF(seq, frame-1)
    print(f'First ORF in frame {frame} is {firstorflength} of length {(firstorflength[1] - firstorflength[0]) / 3} codons.')
print('------------------------------------------------')
        
print('\nCDS ranges:\nCDS join(266..13468,13468..21555)\nCDS 266..13483')
print('CDS 21563..25384\nCDS 25393..26220\nCDS 26245..26472\nCDS 26523..27191')
print('CDS 27202..27387\nCDS 27394..27759\nCDS 27756..27887\nCDS 27894..28259\nCDS 28274..29533\nCDS 29558..29674\n')

cds = [(266,13468),(13468,21555),(266, 13483),(21563,25384),(25393,26220),(26245,26472),
        (26523,27191),(27202,27387),(27394,27759),(27756,27887),(27894,28259),(28274,29533),(29558,29674)]

f = open('output.txt', 'w+')

for frame in range(1, 4) :
    start = frame - 1
    while start != None:
        orflen = firstORF(seq, start)
        if orflen[0] == None or orflen[1] == None:
            break
        if((orflen[1] - orflen[0]) > 100):
            found = [item for item in cds if item[0] == (orflen[0]+1) or item[1] == (orflen[0]+1)
                    or item[0] == orflen[1] or item[1] == orflen[1]]
            if found:
                print(f'Matching ORF line in frame {frame}: {orflen[0] + 1}..{orflen[1]} with CDS: {found[0][0]}..{found[0][1]}\n')
                f.write(f'Matching ORF line in frame {frame}: {orflen[0] + 1}..{orflen[1]} with CDS: {found[0][0]}..{found[0][1]}\n')
        start = orflen[1]

f.close()