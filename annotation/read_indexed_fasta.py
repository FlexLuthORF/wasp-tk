# Read the specified sequence from an indexed fasta file


def read_sequence_from_fasta(fasta_file, index_file, sequence_id, start=1, end=None):
    # Read the index file to find the position of the sequence
    with open(index_file, 'r') as f_index:
        for line in f_index:
            line_parts = line.strip().split('\t')
            if line_parts[0] == sequence_id:
                seq_length = int(line_parts[1])     # length of the sequence
                seq_offset = int(line_parts[2])     # offset of the sequence in the file
                line_bases = int(line_parts[3])     # number of bases on each line
                line_width = int(line_parts[4])     # number of bytes on each line, including newline chars
                break
        else:
            raise ValueError("Sequence ID not found in index file")

    if end is None:
        end = start + seq_length - 1

    if end > seq_length:
        raise ValueError("End position is beyond the end of the sequence")

    # Calculate the position in the file
    line_width = line_bases + 1     # newline is always 1 character in Python
    start_line = (start - 1) // line_bases
    start_line_offset = (start - 1) % line_bases
    end_line = (end - 1) // line_bases
    end_line_offset = (end - 1) % line_bases

    pos = seq_offset + start_line * line_width + start_line_offset
    end_pos = seq_offset + end_line * line_width + end_line_offset

    # Read the sequence from the file
    sequence = ''
    with open(fasta_file, 'r') as f_fasta:
        f_fasta.seek(pos)
        sequence = f_fasta.read(end_pos - pos + 1).replace('\n', '')

    return sequence


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file')
    parser.add_argument('index_file')
    parser.add_argument('sequence_id')
    parser.add_argument('--start', type=int, default=1)
    parser.add_argument('--end', type=int, default=None)
    args = parser.parse_args()

    sequence = read_sequence_from_fasta(args.fasta_file, args.index_file, args.sequence_id, args.start, args.end)
    print(sequence)
