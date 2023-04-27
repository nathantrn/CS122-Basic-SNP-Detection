#project1a
# Nathan Tran

KMER_SIZE = 48
L_third = int(KMER_SIZE/3)

#read in reads file and puts it into a list
def read_in_reads():
    
    f = open('project1a_10000_with_error_paired_reads.fasta', 'r')
    raw = f.read().strip().split("\n")        
    f.close()
    all_reads = []
    for item in raw:
        if item[0] != ">":
            all_reads.append(item)
    return all_reads

#reads in genome reference and puts it into a string
def read_in_genome():
    f = open('project1a_10000_reference_genome.fasta', 'r')
    raw = f.read().splitlines()     
    f.close()
    ref_genome = ''
    for line in raw:
        if line[0] != ">":
            ref_genome = ref_genome + line
    return ref_genome

#divides genome into kmers and records indexes using a dict
def index_kmers_in_genome(reference_genome):
    index_of_kmers = {}
    for i in range(len(reference_genome) - L_third  + 1):
        kmer = reference_genome[i:(i + L_third )]
        if kmer in index_of_kmers:
            index_of_kmers[kmer].append(i)
        else:
            index_of_kmers[kmer] = [i]
    return index_of_kmers
    
# divides the reads into kmers
def convert_reads_to_kmers(reads):
    kmers = []
    for item in reads:
        for i in range(len(item) - KMER_SIZE + 1):
            kmers.append(item[i:(i + KMER_SIZE)])
    return kmers

# dictionary of kmer frequencies
def kmer_frequencies(kmers):
    kmer_frequenciesDICT = {}
    for kmer in kmers:
        if kmer in kmer_frequenciesDICT:
            kmer_frequenciesDICT[kmer] += 1
        else:
            kmer_frequenciesDICT[kmer] = 1
    return kmer_frequenciesDICT

# compare kmer and matching section, any mismatches can be considered a potential mutation
def mutDICT(kmer, matching_sections):
    mutDICT = {}
    for ref_section, start_index in matching_sections:
        for i in range(KMER_SIZE):
            if kmer[i] != ref_section[i]:
                ref_base = ref_section[i]
                m_base = kmer[i]
                ID = start_index + i
                mut = [ID, ref_base, m_base]
                if ID not in mutDICT:
                    mutDICT[ID] = mut
    return mutDICT

# compare kmer and matching section, record matches in a dict
# used for read error analysis
def normDICT(kmer, matching_sections):
    normDICT = {}
    for ref_section, start_index in matching_sections:
        for i in range(KMER_SIZE):
            if kmer[i] == ref_section[i]:
                ref_base = ref_section[i]
                r_base = kmer[i]
                ID = start_index + i
                norm = [ID, ref_base, r_base]
                if ID not in normDICT:
                    normDICT[ID] = norm

    return normDICT


def get_MUTS(kmers, reference_index, reference):
    MUTS = {}
    for kmer in kmers:
        potential_indexes = potential_IDS(kmer, reference_index)
        potential_indexes = list(dict.fromkeys(potential_indexes))
        matching_sections = get_matching_sections(kmer, reference, potential_indexes)
        potential_MUT = mutDICT(kmer, matching_sections)
        for index in potential_MUT:
            if index not in MUTS:
                MUTS[index] = potential_MUT[index]
    return MUTS

# mutation frequency
# used for read error analysis
# compare mutation frequency of index to normal frequency of index
def MUT_freq(kmers, reference_index, reference):
    freq_MUT = {}
    for kmer in kmers:
        potential_indexes = potential_IDS(kmer, reference_index)
        potential_indexes = list(dict.fromkeys(potential_indexes))
        matching_sections = get_matching_sections(kmer, reference, potential_indexes)
        potential_MUT = mutDICT(kmer, matching_sections)
        for index in potential_MUT:
            if index not in freq_MUT:
                freq_MUT[index] = 1
            else:
                freq_MUT[index] += 1
   
    return freq_MUT

# norm frequency
# used for read error analysis
# compare mutation frequency of index to normal frequency of index
def NORM_freq(kmers, reference_index, reference):
    freq_NORM = {}
    for kmer in kmers:
        potential_indexes = potential_IDS(kmer, reference_index)
        potential_indexes = list(dict.fromkeys(potential_indexes))
        matching_sections = get_matching_sections(kmer, reference, potential_indexes)
        norms = normDICT(kmer, matching_sections)
        for index in norms:
            if index not in freq_NORM:
                freq_NORM[index] = 1
            else:
                freq_NORM[index] += 1
    return freq_NORM

#matches kmers to potential starting positions
def potential_IDS(kmer, reference_index):
    potential_indexes = []
    first_third = kmer[0:L_third]
    second_third = kmer[L_third :(2 * L_third )]
    third_third = kmer[(2 * L_third):]
    if first_third in reference_index:
        potential_indexes += reference_index[first_third]
    if second_third in reference_index:
        second_align = [(index - L_third) for index in reference_index[second_third]]
        potential_indexes += second_align
    if third_third in reference_index:
        third_align = [(index - (2 * L_third )) for index in reference_index[third_third]]
        potential_indexes += third_align
    return potential_indexes

#returns matching sections of kmers and reference
def get_matching_sections(kmer, reference, potential_indexes):
    matching_sections = []
    for index in potential_indexes:
        matching_section = reference[index:(index + KMER_SIZE)]
        num_mismatches = sum(kmer[i] != matching_section[i] for i in range(KMER_SIZE))
        if num_mismatches > 0 and num_mismatches <= 2:
            matching_sections.append((matching_section, index))
    return matching_sections


if __name__ == "__main__":
   
    reads = read_in_reads()
    reference_genome = read_in_genome()
    reference_index = index_kmers_in_genome(reference_genome)
    kmers = convert_reads_to_kmers(reads)
    kmer_to_frequency = kmer_frequencies(kmers)
    kmers = list(kmer_to_frequency.keys())
    freqMUT = MUT_freq(kmers, reference_index, reference_genome)
    freqNORM = NORM_freq(kmers, reference_index, reference_genome)
    MUTS = get_MUTS(kmers, reference_index, reference_genome)
    possible_indel = {}
    for key in freqMUT:
        if freqMUT[key] <= 3:
            possible_indel[key] = MUTS[key]
            del MUTS[key]
            
    for key in freqMUT:
        if key in freqNORM:
            if freqNORM[key] > freqMUT[key]:
                if key in MUTS:
                    del MUTS[key]
    snps = list(MUTS.values())
    snps.sort(key=lambda x: x[0])
    for x in snps:
        print(">S" + str(x[0]) + " " + str(x[1]) + " " + str(x[2]))
 


