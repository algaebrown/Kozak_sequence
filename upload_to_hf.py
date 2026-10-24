from datasets import Dataset, DatasetDict


def upload_kozak_to_hf(
    fasta_path="kozak.fasta", repo_id="algaebrown/human-kozak-21bp"
):
    from Bio import SeqIO

    # Parse headers and sequences into a dictionary format
    headers = []
    sequences = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq).upper())

    # Create Hugging Face Dataset
    raw_dict = {"id": headers, "sequence": sequences}
    hf_dataset = Dataset.from_dict(raw_dict)

    # Push to Hugging Face Hub
    print(f"Uploading dataset to Hugging Face: {repo_id}...")
    hf_dataset.push_to_hub(repo_id, private=False)
    print("Upload complete!")


# Usage (Replace with your Hugging Face username/repo)
upload_kozak_to_hf("kozak.fasta", repo_id="algaebrown/human-kozak-21bp")