"""
Transcription and translation take the information in DNA and use it to produce proteins.
Transcription uses a strand of DNA as a template to build a molecule called RNA.
The RNA molecule is the link between DNA and the production of proteins.
During translation, the RNA molecule created in the transcription process delivers information
from the DNA to the protein-building machines.
                      DNA → RNA → Protein
"""

print("Welcome to Nucleic Acid Converter")
print("This machine will convert DNA to mRNA to Proteins from a fasta file \n")

def openfile(filename_fasta):
    """"
    Retrieve  a DNA from a fasta file

    Parameters:
    ___________
    fasta file : DNA sequence
        A fasta file contain the DNA sequence

    Returns:
    ________
    RNA : RNA sequence
        A RNA sequence contains AUCG

    """
    global initial_seq

    initial_seq= ""
    with open(filename_fasta) as file:
        for line in file:
            if line[0] not in ['>', ';']:
                initial_seq += line.strip()
        print("This is the DNA sequence : \n")
        print(initial_seq)
        length = len(initial_seq)
        print("The DNA sequence contains ", length, "nucleotids. ")
        print()

    return initial_seq

def transcription():
    """"
    Transcription is the process of copying a segment of DNA into RNA

    Parameters:
    ___________
    DNA : DNA sequence
        A DNA sequence contains ATCG

    Returns:
    ________
    RNA : RNA sequence
        A RNA sequence contains AUCG

    """
    global initial_seq
    global rna

    final_seq = ""
    #convert the DNA sequence to RNA
    trancribe_dna = initial_seq.maketrans("ATCG", "AUCG")
    rna = initial_seq.translate(trancribe_dna)
    #print the RNA sequence
    print("This is the RNA sequence: \n")
    print(rna)
    print("")

    return final_seq

def translation():
    """
    In translation the mRNA is "decoded" to build a protein

    Parameters:
    ___________
    RNA : RNA sequence
        A RNA sequence contains AUCG
    Returns:
    ________
    Proteins : Proteins  sequence
        A Proteins  sequence contains 'SDRAVAANPASMIIVIM' which denotes differen amino acid.


    """
    #codon table translate a genetic code into a sequence of amino acid
    rna_codons = {
        # 'M' - START CODON , '*' - STOP CODON
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "UGU": "C", "UGC": "C",
        "GAU": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "UUU": "F", "UUC": "F",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAU": "H", "CAC": "H",
        "AUA": "I", "AUU": "I", "AUC": "I",
        "AAA": "K", "AAG": "K",
        "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "AUG": "M", # first codon in the transcribed RNA and it codes for the amino acid methionine
        "AAU": "N", "AAC": "N",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "UGG": "W",
        "UAU": "Y", "UAC": "Y",
        "UAA": "*", "UAG": "*", "UGA": "*"
    }

    global rna

    protein = ""
    for i in range(0, len(rna), 3):
        codon = rna[i:i + 3]
        protein += rna_codons[codon]

    print("This is the protein: \n")
    print(protein)

    return protein

def main():
    #open a fasta or a txt file
    openfile("Agrobacterium.fasta")
    #transcription of DNA
    transcription()
    #translation of mRNA
    translation()

if __name__ == '__main__':
    main()