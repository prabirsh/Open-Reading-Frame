"""Convert a DNA sequence to RNA and translate it to protein
in a chosen forward reading frame (1, 2 or 3)."""

import sys

INPUT_FILE = "prabir_input.txt"
RNA_FILE = "prabir_output1_RNA.txt"
PROTEIN_FILE = "prabir_output2_protein.txt"

# codon dict
RNA_CODON_TABLE = {"UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
                   "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
                   "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
                   "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
                   "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
                   "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
                   "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
                   "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
                   "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
                   "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
                   "UAA": "STOP", "CAA": "Q", "AAA": "K", "GAA": "E",
                   "UAG": "STOP", "CAG": "Q", "AAG": "K", "GAG": "E",
                   "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
                   "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
                   "UGA": "STOP", "CGA": "R", "AGA": "R", "GGA": "G",
                   "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"}


def read_dna(fname):
    """Read a DNA sequence from a file.

    All whitespace (including newlines) is removed and the sequence is
    upper-cased. Raises ValueError if the file is empty or contains
    anything other than A, C, G, T.
    """
    with open(fname) as fh:
        dna = "".join(fh.read().split()).upper()
    if not dna:
        raise ValueError(f"{fname} contains no sequence")
    invalid = sorted(set(dna) - set("ACGT"))
    if invalid:
        raise ValueError(f"{fname} contains non-ACGT characters: {', '.join(invalid)}")
    return dna


def transcribe(dna):
    """Replace all occurrences of T with U."""
    return dna.replace("T", "U")


def translate(rna, frame):
    """Translate an RNA string in forward reading frame 1, 2 or 3.

    Stop codons are written as '*'. An incomplete trailing codon is ignored.
    """
    protein = []
    for i in range(frame - 1, len(rna) - 2, 3):
        amino_acid = RNA_CODON_TABLE[rna[i:i + 3]]
        protein.append("*" if amino_acid == "STOP" else amino_acid)
    return "".join(protein)


def ask_frame():
    """Prompt until the user enters a valid reading frame (1, 2 or 3)."""
    print("\n-::::::::::::::Protein:::::::::::::::-")
    while True:
        print("Enter 1 :-for Reading Frame 1 on the direct strand")
        print("Enter 2 :-for Reading Frame 2 on the direct strand")
        print("Enter 3 :-for Reading Frame 3 on the direct strand")
        print("Enter your Reading Frame:")
        try:
            frame = int(input())
        except ValueError:
            frame = 0
        if frame in (1, 2, 3):
            return frame
        print("INPUT//Invalid!! Please enter 1, 2 or 3.\n")


def main():
    fname = sys.argv[1] if len(sys.argv) > 1 else INPUT_FILE
    try:
        dna = read_dna(fname)
    except (OSError, ValueError) as err:
        sys.exit(f"Error: {err}")

    print("length of the DNA:-", len(dna))

    # generate RNA string and RNA file
    rna = transcribe(dna)
    with open(RNA_FILE, "w") as fout:
        fout.write(rna + "\n")

    print("\n -::::::::::::::::RNA:::::::::::::::-")
    print(">RNA Sequence from base 0 to ", len(dna))
    print(rna)
    print(f"\n \n :::-File:{RNA_FILE} created on your directory where RNA sequences are stored-::: \n \n")

    # generate protein string and protein file
    frame = ask_frame()
    protein = translate(rna, frame)
    with open(PROTEIN_FILE, "w") as fout:
        fout.write(protein + "\n")

    print(f">Protein Sequence from base {frame - 1} to ", len(dna), f"bp (reading frame {frame})")
    print(protein)
    print(f"\n \n :::-File:{PROTEIN_FILE} created on your directory where Protein sequences are stored-::: \n \n")


if __name__ == "__main__":
    main()
