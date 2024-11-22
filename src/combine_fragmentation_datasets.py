from Bio import SeqIO
import random
import sys
import os

# Add the src directory to the system path
sys.path.append(os.path.dirname(__file__))

# Import the function
from concatenate_windows import parse_header


def main():
    random.seed(1)
    assembly_file1 = sys.argv[1]  # Path to your input assembly file1
    assembly_file2 = sys.argv[2]  # Path to your input assembly file2
    output_file = sys.argv[3]    # Path to your desired output file

    # Load the assembly file
    records1 = list(SeqIO.parse(assembly_file1, "fasta"))
    records2 = list(SeqIO.parse(assembly_file2, "fasta"))

    samples1 = set([parse_header(r.id)[0] for r in records1])
    samples2 = set([parse_header(r.id)[0] for r in records2])

    assert samples1 == samples2, "Not the same samples in assemblies"

    num_samples = len(samples1)
    print(num_samples)

    # CVM28_Exiguobacterium_aurantiacum was not concatenated to a long contig so I will force this into the framented assembly
    # {'CVM28_Exiguobacterium_aurantiacum', 'CVM64_Rhizobium_leguminosarum', 'DSMZ16553-Burkholderia_cenocepacia', 'DSMZ2661-Methanocaldococcus_jannaschii'}
    samples1.remove("CVM28_Exiguobacterium_aurantiacum")
    samples1.remove("CVM64_Rhizobium_leguminosarum")
    samples1.remove("DSMZ16553-Burkholderia_cenocepacia")
    samples1.remove("DSMZ2661-Methanocaldococcus_jannaschii")
    print(len(samples2))
    selected_samples = random.sample(list(samples1), num_samples // 2)
    

    sequences_to_write = []

    for record in records1:
        sample_id = parse_header(record.id)[0]

        if sample_id not in selected_samples:
            sequences_to_write.append(record)
    for record in records2:
        sample_id = parse_header(record.id)[0]

        if sample_id in selected_samples:
            sequences_to_write.append(record)

    samples_after = set([parse_header(r.id)[0] for r in sequences_to_write])
    assert samples_after == samples2

    SeqIO.write(sequences_to_write, output_file, "fasta")
    
if __name__ == '__main__':
    main()
    
