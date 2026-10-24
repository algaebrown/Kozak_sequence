from Bio.Seq import Seq
import gzip
import yaml
import pysam


def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def extract_kozak_from_gff(
    config_path="config.yaml",
    output_fa="kozak.fasta",
    upstream=10,
    downstream=7,
):
    # Load paths from config.yaml
    config = load_config(config_path)
    gff_path = config["GFF"]
    genome_fa_path = config["GENOME"]

    print(f"Opening genome: {genome_fa_path}")
    fa = pysam.FastaFile(genome_fa_path)

    extracted_count = 0
    skipped_count = 0

    print(f"Processing GFF annotations from: {gff_path}")
    
    # Open gzipped GFF3 file
    opener = gzip.open if gff_path.endswith(".gz") else open

    with opener(gff_path, "rt") as gtf, open(output_fa, "w") as out_f:
        for line in gtf:
            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attributes = parts

            # Gencode GFF3 uses 'start_codon'
            if feature == "start_codon":
                g_start = int(start)
                g_end = int(end)

                # Validate start codon length is 3 bp
                if (g_end - g_start + 1) != 3:
                    skipped_count += 1
                    continue

                try:
                    if strand == "+":
                        # Upstream is to the left (lower coordinates)
                        fetch_start = g_start - 1 - upstream  # 0-based start
                        fetch_end = g_end + downstream  # 0-based end (exclusive)
                        seq = fa.fetch(chrom, fetch_start, fetch_end).upper()

                    elif strand == "-":
                        # On negative strand, upstream is towards higher coordinates
                        fetch_start = g_start - 1 - downstream
                        fetch_end = g_end + upstream
                        raw_seq = fa.fetch(chrom, fetch_start, fetch_end).upper()
                        # Reverse complement to maintain 5' -> 3' biological orientation
                        seq = str(Seq(raw_seq).reverse_complement())

                    # Quality check: verify correct length and core start codon triplet
                    expected_len = upstream + 3 + downstream
                    core_start = seq[upstream : upstream + 3]

                    if len(seq) == expected_len and core_start in {
                        "ATG",
                        "CTG",
                        "GTG",
                    }:
                        out_f.write(
                            f">chrom_{chrom}_pos_{g_start}_strand_{strand}\n{seq}\n"
                        )
                        extracted_count += 1
                    else:
                        skipped_count += 1

                except KeyError:
                    # Handles chromosomes present in GFF but missing from genome fasta index
                    skipped_count += 1
                    continue

    fa.close()
    print(f"\nDone!")
    print(f"Successfully extracted: {extracted_count} Kozak windows")
    print(f"Skipped: {skipped_count} entries (mismatched lengths or contigs)")
    print(f"Output saved to: {output_fa}")


if __name__ == "__main__":
    extract_kozak_from_gff()