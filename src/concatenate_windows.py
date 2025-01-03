from Bio import SeqIO
from Bio.Seq import Seq
import sys

def parse_header(header):
    """
    Parses the header to extract sample name, contig ID, and sliding window range.
    """
    # Split header by '_sliding_'
    if '_sliding_' in header:
        sample_contig_part, sliding_window = header.rsplit('_sliding_', 1)
    else:
        sample_contig_part = header
        sliding_window = None

    # Split sample_contig_part by '_contig_'
    if '_contig_' in sample_contig_part:
        sample_name, contig_id = sample_contig_part.rsplit('_contig_', 1)
    else:
        sample_name = sample_contig_part
        contig_id = None

    return sample_name, contig_id, sliding_window

def parse_sliding_window(sliding_window):
    """
    Parses the sliding window string to get start and end positions as integers.
    """
    start_str, end_str = sliding_window.split('-')
    return int(start_str), int(end_str)

def main():
    num_windows = 80
    # Input and output file paths
    assembly_file = sys.argv[1]  # Path to your input assembly file
    output_file = sys.argv[2]    # Path to your desired output file

    # Load the assembly file
    records = list(SeqIO.parse(assembly_file, "fasta"))

    # Dictionary to store slides organized by sample and contig
    slides_by_sample_and_contig = {}

    # Group slides by sample and contig
    for record in records:
        header = record.id
        sample_name, contig_id, sliding_window = parse_header(header)

        if sample_name not in slides_by_sample_and_contig:
            slides_by_sample_and_contig[sample_name] = {}
        if contig_id not in slides_by_sample_and_contig[sample_name]:
            slides_by_sample_and_contig[sample_name][contig_id] = []

        slides_by_sample_and_contig[sample_name][contig_id].append((sliding_window, record))

    print(f"There are {len(slides_by_sample_and_contig.keys())} samples")
    all_sequences_to_write = []
    samples_w_concatenation = []

    
    # Process each sample and contig
    for sample_name in slides_by_sample_and_contig:
        concatenated = False
        for contig_id in slides_by_sample_and_contig[sample_name]:
            slides = slides_by_sample_and_contig[sample_name][contig_id]
            # Sort the slides by start position
            slides.sort(key=lambda x: parse_sliding_window(x[0])[0])

            num_slides = len(slides)

            used_slide_indices = set()

            if not concatenated and num_slides >= num_windows:
                i = 0
                while i <= num_slides - num_windows:
                    # Check if the next num_windows slides are in succession
                    consecutive = True
                    for j in range(i, i + num_windows - 1):
                        start1, end1 = parse_sliding_window(slides[j][0])
                        start2, end2 = parse_sliding_window(slides[j + 1][0])

                        # Check if the windows are consecutive
                        if end1 + 1 != start2:
                            consecutive = False
                            break

                    if consecutive:
                        # Concatenate the sequences of the num_windows consecutive slides
                        sequences = [slides[k][1].seq for k in range(i, i + num_windows)]
                        concatenated_seq = ''.join(map(str, sequences))

                        # Create a new header for the concatenated sequence
                        start_position = parse_sliding_window(slides[i][0])[0]
                        end_position = parse_sliding_window(slides[i + num_windows-1][0])[1]

                        new_header = f"{sample_name}_contig_{contig_id}_concatenated_{start_position}-{end_position}"
                        new_record = SeqIO.SeqRecord(
                            seq=Seq(concatenated_seq),
                            id=new_header,
                            description=''
                        )

                        # Add the concatenated sequence to the list
                        all_sequences_to_write.append(new_record)

                        # Mark these slides as used
                        for k in range(i, i + num_windows):
                            used_slide_indices.add(k)

                        concatenated = True
                        samples_w_concatenation.append(sample_name)
                        break
                    else:
                        # If not consecutive, move to the next slide
                        i += 1

            # Collect the unused slides
            for idx in range(num_slides):
                if idx not in used_slide_indices:
                    record = slides[idx][1]
                    all_sequences_to_write.append(record)


    samples_without_concatenation = set(slides_by_sample_and_contig.keys()) - set(samples_w_concatenation)
    print("Samples without a concatenation")
    print(samples_without_concatenation)

    
    # Write all sequences (concatenated and unused) to the output file
    SeqIO.write(all_sequences_to_write, output_file, "fasta")

if __name__ == "__main__":
    main()

