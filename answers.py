#Gene to test with from https://www.ncbi.nlm.nih.gov/gene/?term=Homo+sapiens+insulin
insulin_gene = """AGCCCTCCAGGACAGGCTGCATCAGAAGAGGCCATCAAGCAGGTCTGTTCCAAGGGCCTTTGCGTCAGGTGGGCTCAGGATTCCAGGGTGGCTGGACCCCAGGCCCCAGCTCTGCAGCAGGGAGGACGTGGCTGGGCTCGTGAAGCATGTGGGGGTGAGCCCAGGGGCCCCAAGGCAGGGCACCTGGCCTTCAGCCTGCCTCAGCCCTGCCTGTCTCCCAGATCACTGTCCTTCTGCCATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGGTGAGCCAACTGCCCATTGCTGCCCCTGGCCGCCCCCAGCCACCCCCTGCTCCTGGCGCTCCCACCCAGCATGGGCAGAAGGGGGCAGGAGGCTGCCACCCAGCAGGGGGTCAGGTGCACTTTTTTAAAAAGAAGTTCTCTTGGTCACGTCCTAAAAGTGACCAGCTCCCTGTGGCCCAGTCAGAATCTCAGCCTGAGGACGGTGTTGGCTTCGGCAGCCCCGAGATACATCAGAGGGTGGGCACGCTCCTCCCTCCACTCGCCCCTCAAACAAATGCCCCGCAGCCCATTTCTCCACCCTCATTTGATGACCGCAGATTCAAGTGTTTTGTTAAGTAAAGTCCTGGGTGACCTGGGGTCACAGGGTGCCCCACGCTGCCTGCCTCTGGGCGAACACCCCATCACGCCCGGAGGAGGGCGTGGCTGCCTGCCTGAGTGGGCCAGACCCCTGTCGCCAGGCCTCACGGCAGCTCCATAGTCAGGAGATGGGGAAGATGCTGGGGACAGGCCCTGGGGAGAAGTACTGGGATCACCTGTTCAGGCTCCCACTGTGACGCTGCCCCGGGGCGGGGGAAGGAGGTGGGACATGTGGGCGTTGGGGCCTGTAGGTCCACACCCAGTGTGGGTGACCCTCCCTCTAACCTGGGTCCAGCCCGGCTGGAGATGGGTGGGAGTGCGACCTAGGGCTGGCGGGCAGGCGGGCACTGTGTCTCCCTGACTGTGTCCTCCTGTGTCCCTCTGCCTCGCCGCTGTTCCGGAACCTGCTCTGCGCGGCACGTCCTGGCAGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCCGCCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC"""

def find_stop_codons(seq):
    #table of amino acids
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein_location = "start indices "
    # check length of DNA early on
    if len(seq)%3 == 0:
        # iterate over the sequence and take a step of 3 each time 
        for i in range(0, len(seq), 3):
            # find the current codon based on the value
            if seq[i:i + 3] in ['TAG', 'TAA', 'TGA']:
              protein_location += f'{i}, '
    print(protein_location)
    return protein_location

def translate_to_string(seq):
    #table of amino acids
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    # check length of DNA early on
    if len(seq)%3 == 0:
        # iterate over the sequence and take a step of 3 each time 
        for i in range(0, len(seq), 3):
            # find the current codon based on the dictionary key-value lookup
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein

  
def translate_to_list(seq):
    #table of amino acids
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = []
    # check length of DNA early on
    if len(seq)%3 == 0:
        # iterate over the sequence and take a step of 3 each time 
        for i in range(0, len(seq), 3):
            # find the current codon based on the dictionary key-value lookup
            codon = seq[i:i + 3]
            protein.append(table[codon])
    print(protein)
    return protein


find_stop_codons(insulin_gene)
# output: start indices 294, 426, 900, 1140, 1413, 1422, 
translate_to_string(insulin_gene)
# output: SPPGQAASEEAIKQVCSKGLCVRWAQDSRVAGPQAPALQQGGRGWAREACGGEPRGPKAGHLAFSLPQPCLSPRSLSFCHGPVDAPPAPAGAAGPLGT_PSRSLCEPTPVRLTPGGSSLPSVRGTRLLLHTQDPPGGRGPAG_ANCPLLPLAAPSHPLLLALPPSMGRRGQEAATQQGVRCTFLKRSSLGHVLKVTSSLWPSQNLSLRTVLASAAPRYIRGWARSSLHSPLKQMPRSPFLHPHLMTADSSVLLSKVLGDLGSQGAPRCLPLGEHPITPGGGRGCLPEWARPLSPGLTAAP_SGDGEDAGDRPWGEVLGSPVQAPTVTLPRGGGRRWDMWALGPVGPHPVWVTLPLTWVQPGWRWVGVRPRAGGQAGTVSP_LCPPVSLCLAAVPEPALRGTSWQWGRWSWAGALVQAACSPWPWRGPCRSVALWNNAVPASAPSTSWRTTATRRSPQAAPHPPPPAPREME_SP_TS
translate_to_list(insulin_gene)
# output: ['S', 'P', 'P', 'G', 'Q', 'A', 'A', 'S', 'E', 'E', 'A', 'I', 'K', 'Q', 'V', 'C', 'S', 'K', 'G', 'L', 'C', 'V', 'R', 'W', 'A', 'Q', 'D', 'S', 'R', 'V', 'A', 'G', 'P', 'Q', 'A', 'P', 'A', 'L', 'Q', 'Q', 'G', 'G', 'R', 'G', 'W', 'A', 'R', 'E', 'A', 'C', 'G', 'G', 'E', 'P', 'R', 'G', 'P', 'K', 'A', 'G', 'H', 'L', 'A', 'F', 'S', 'L', 'P', 'Q', 'P', 'C', 'L', 'S', 'P', 'R', 'S', 'L', 'S', 'F', 'C', 'H', 'G', 'P', 'V', 'D', 'A', 'P', 'P', 'A', 'P', 'A', 'G', 'A', 'A', 'G', 'P', 'L', 'G', 'T', '_', 'P', 'S', 'R', 'S', 'L', 'C', 'E', 'P', 'T', 'P', 'V', 'R', 'L', 'T', 'P', 'G', 'G', 'S', 'S', 'L', 'P', 'S', 'V', 'R', 'G', 'T', 'R', 'L', 'L', 'L', 'H', 'T', 'Q', 'D', 'P', 'P', 'G', 'G', 'R', 'G', 'P', 'A', 'G', '_', 'A', 'N', 'C', 'P', 'L', 'L', 'P', 'L', 'A', 'A', 'P', 'S', 'H', 'P', 'L', 'L', 'L', 'A', 'L', 'P', 'P', 'S', 'M', 'G', 'R', 'R', 'G', 'Q', 'E', 'A', 'A', 'T', 'Q', 'Q', 'G', 'V', 'R', 'C', 'T', 'F', 'L', 'K', 'R', 'S', 'S', 'L', 'G', 'H', 'V', 'L', 'K', 'V', 'T', 'S', 'S', 'L', 'W', 'P', 'S', 'Q', 'N', 'L', 'S', 'L', 'R', 'T', 'V', 'L', 'A', 'S', 'A', 'A', 'P', 'R', 'Y', 'I', 'R', 'G', 'W', 'A', 'R', 'S', 'S', 'L', 'H', 'S', 'P', 'L', 'K', 'Q', 'M', 'P', 'R', 'S', 'P', 'F', 'L', 'H', 'P', 'H', 'L', 'M', 'T', 'A', 'D', 'S', 'S', 'V', 'L', 'L', 'S', 'K', 'V', 'L', 'G', 'D', 'L', 'G', 'S', 'Q', 'G', 'A', 'P', 'R', 'C', 'L', 'P', 'L', 'G', 'E', 'H', 'P', 'I', 'T', 'P', 'G', 'G', 'G', 'R', 'G', 'C', 'L', 'P', 'E', 'W', 'A', 'R', 'P', 'L', 'S', 'P', 'G', 'L', 'T', 'A', 'A', 'P', '_', 'S', 'G', 'D', 'G', 'E', 'D', 'A', 'G', 'D', 'R', 'P', 'W', 'G', 'E', 'V', 'L', 'G', 'S', 'P', 'V', 'Q', 'A', 'P', 'T', 'V', 'T', 'L', 'P', 'R', 'G', 'G', 'G', 'R', 'R', 'W', 'D', 'M', 'W', 'A', 'L', 'G', 'P', 'V', 'G', 'P', 'H', 'P', 'V', 'W', 'V', 'T', 'L', 'P', 'L', 'T', 'W', 'V', 'Q', 'P', 'G', 'W', 'R', 'W', 'V', 'G', 'V', 'R', 'P', 'R', 'A', 'G', 'G', 'Q', 'A', 'G', 'T', 'V', 'S', 'P', '_', 'L', 'C', 'P', 'P', 'V', 'S', 'L', 'C', 'L', 'A', 'A', 'V', 'P', 'E', 'P', 'A', 'L', 'R', 'G', 'T', 'S', 'W', 'Q', 'W', 'G', 'R', 'W', 'S', 'W', 'A', 'G', 'A', 'L', 'V', 'Q', 'A', 'A', 'C', 'S', 'P', 'W', 'P', 'W', 'R', 'G', 'P', 'C', 'R', 'S', 'V', 'A', 'L', 'W', 'N', 'N', 'A', 'V', 'P', 'A', 'S', 'A', 'P', 'S', 'T', 'S', 'W', 'R', 'T', 'T', 'A', 'T', 'R', 'R', 'S', 'P', 'Q', 'A', 'A', 'P', 'H', 'P', 'P', 'P', 'P', 'A', 'P', 'R', 'E', 'M', 'E', '_', 'S', 'P', '_', 'T', 'S']
