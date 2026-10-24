from Bio import SeqIO
import logomaker
import matplotlib.pyplot as plt
import pandas as pd


def generate_kozak_logo(fasta_path="kozak.fasta", output_plot="kozak_motif.png"):
    print(f"Reading sequences from {fasta_path}...")
    sequences = [str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")]

    if not sequences:
        print("Error: No sequences found in FASTA file.")
        return

    print(f"Loaded {len(sequences):,} sequences. Building frequency matrix...")

    # Define alphabet
    chars = ["A", "C", "G", "T"]
    seq_len = len(sequences[0])

    # Initialize a count matrix (rows = positions, columns = nucleotides)
    counts = pd.DataFrame(0.0, index=range(seq_len), columns=chars)

    for seq in sequences:
        if len(seq) == seq_len:
            for i, nt in enumerate(seq):
                if nt in chars:
                    counts.loc[i, nt] += 1.0

    # Convert counts to frequency matrix (probabilities)
    freq_matrix = counts.div(counts.sum(axis=1), axis=0)

    # Create the Sequence Logo plot
    fig, ax = plt.subplots(figsize=(12, 4))

    # Use logomaker to draw the logo
    logo = logomaker.Logo(
        freq_matrix,
        ax=ax,
        color_scheme="classic",  # Standard nucleotide colors (A:green, C:blue, G:orange, T:red)
        stack_order="big_on_top",
    )

    # Style the plot
    ax.set_title(
        f"Kozak Sequence Logo (n = {len(sequences):,})", fontsize=14, fontweight="bold"
    )
    ax.set_xlabel("Position Relative to Start Codon", fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)

    # Customize X-axis ticks so 'ATG' is clearly labeled in the center
    # Assuming upstream=10, ATG is at indices 10, 11, 12
    tick_positions = list(range(seq_len))
    # Label relative positions: -10 to +7 (or adjust based on your slicing)
    tick_labels = [str(i - 10) for i in tick_positions]
    # Highlight ATG positions
    tick_labels[10] = "A(0)"
    tick_labels[11] = "T(+1)"
    tick_labels[12] = "G(+2)"

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=0, fontsize=10)

    # Add grid lines and clean layout
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()

    # Save figure
    plt.savefig(output_plot, dpi=300)
    print(f"Sequence logo successfully saved to {output_plot}")
    plt.show()


if __name__ == "__main__":
    generate_kozak_logo()