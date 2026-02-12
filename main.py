# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # Align all species to humans using BLOSUM62, gap open = -10, gap extend = -1
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)

    species = [
        ("Gallus gallus", gg_seq),
        ("Mus musculus", mm_seq),
        ("Balaeniceps rex", br_seq),
        ("Tursiops truncatus", tt_seq),
    ]

    results = []
    for name, seq in species:
        score, _, _ = nw.align(hs_seq, seq)
        results.append((name, score))

    # Sort most similar to least similar
    results.sort(key=lambda x: x[1], reverse=True)

    # Print species ranked by similarity and their alignment scores
    print("Species BRD2 ranked by similarity to Human BRD2:")
    for name, score in results:
        print(f"  {name}: {score}")
    

if __name__ == "__main__":
    main()
