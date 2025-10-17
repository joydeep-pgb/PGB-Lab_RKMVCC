from Bio import SeqIO
import os

def split_fasta(input_fasta, output_dir, chunk_size=100):
    """
    Split a FASTA file into smaller files each containing `chunk_size` sequences.

    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_dir (str): Directory to store the split FASTA files.
        chunk_size (int): Number of sequences per output file (default = 100).
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    records = SeqIO.parse(input_fasta, "fasta")
    batch = []
    file_count = 1

    for i, record in enumerate(records, start=1):
        batch.append(record)

        # Write batch when reaching the chunk size
        if i % chunk_size == 0:
            output_file = os.path.join(output_dir, f"subset_{file_count}.fasta")
            SeqIO.write(batch, output_file, "fasta")
            print(f"✅ Created {output_file} with {len(batch)} sequences.")
            batch = []
            file_count += 1

    # Write remaining sequences if total % chunk_size != 0
    if batch:
        output_file = os.path.join(output_dir, f"subset_{file_count}.fasta")
        SeqIO.write(batch, output_file, "fasta")
        print(f"✅ Created {output_file} with {len(batch)} sequences.")

if __name__ == "__main__":
    # Example usage
    input_fasta = "/media/pgb-lab/One_HDD/Mango_AS_rMATS/IsoformSwitch/Seq_Output/isoformSwitchAnalyzeR_isoform_AA_complete.fasta"
    output_dir = "/media/pgb-lab/One_HDD/Mango_AS_rMATS/IsoformSwitch/Seq_Output/Split_FASTA/"
    split_fasta(input_fasta, output_dir, chunk_size=100)
    print("All done!")