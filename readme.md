# Human Kozak Sequence Dataset (21 bp) 🧬

A lightweight, high-signal genomic dataset of **21-nucleotide translation initiation windows** extracted from the human reference genome (`GRCh38`). Designed specifically as an elite toy benchmark for training and evaluating **Sequence Diffusion Models**, categorical language models, and biological sequence deep learning architectures.

---

## 📌 What is this Dataset?

This dataset consists of fixed-length genomic windows centered on annotated protein-coding start codons:
* **Window Composition:** `10 bp upstream` + `ATG` (core start codon) + `7 bp downstream` (Total length: **21 bp**).
* **Biological Relevance:** Mammalian translation initiation relies heavily on strict positional sequence constraints, famously known as the **Kozak Consensus Motif** (`gccRccATGG`, where `R` is a purine at position -3). 

### Why use this for Diffusion LMs?
1. **Perfect Fixed Shape:** Eliminates variable-length padding headaches, formatting neatly into clean tensors of shape `(B, 21, 4)` for continuous embedding diffusion or `(B, 21)` for discrete masking diffusion.
2. **Sharp Positional Invariants:** Features 100% invariant core codons (`ATG`), a hallmark purine preference (`A`/`G`) at position -3, and a downstream `G` at position +4. 
3. **Fast Convergence:** With a tiny alphabet (`A, C, G, T`) and low sequence length, models train in minutes on a single GPU, allowing you to instantly inspect whether your reverse diffusion process learns true biological grammar or just noise.

---

## 📂 Dataset Structure

Each record contains genomic coordinate metadata and the 21-bp nucleotide string:

* `id`: Genomic coordinate string (e.g., `chrom_chr1_pos_123456_strand_+`)
* `sequence`: 21-character string containing uppercase nucleotides (`A`, `C`, `G`, `T`)

---

## 🚀 How to Use

### 1. Load directly via Hugging Face `datasets`

```python
from datasets import load_dataset

# Load your uploaded dataset
dataset = load_dataset("your-username/human-kozak-21bp", split="train")

print(dataset[0])
# Output: {'id': 'chrom_chr1_pos_...', 'sequence': 'GCCGCCATGG...'}