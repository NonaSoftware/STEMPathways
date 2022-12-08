# STEMPathways
This is a learning assignment designed to help students learn to code in Python surrounding basic bioinformatics and the Central Dogma of biology

Below is a small take home assignment including some background reading for converting a DNA sequence into an amino acid sequence. Remember the central dogma of biology that DNA -> RNA -> Protein

[Article 1](https://newatlas.com/science/colossal-woolly-mammoth-de-extinction/).   
[Article 2](https://medlineplus.gov/genetics/understanding/genomicresearch/genomeediting/)

After finishing the articles above, we’ll start off with a dictionary called amino_acids (below) and iterate over a sequence to produce the desired protein. While a typical dictionary wil have you search for a keyword and find the definition(s), a dictionary in python has a unique key that produces a value. For example, if you wanted to access the dictionary amino_acids and get the value from ATA, you would refer to `amino_acids[‘ATA’]` and get the value `I` as below
``` 
amino_acids['ATA']
'I' 
``` 

```
def translate(seq):
    
   amino_acids = {
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
 
   return protein
 ```

1. Find the STOP codon in a given sequence
2. Return each codon in a given sequence
3. Given a sequence, return a list of amino acids
The answers to this assignment are available in this repository
