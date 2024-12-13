from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse
import pandas as pd


def read_lengths_distribution(fasta_file):
    # Read the fasta file and calculate the lengths
    lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    return lengths


def group_read_lengths(lengths):
    # Define bins for grouping read lengths
    bins = [0, 100, 500]  # Initial bins
    max_length = max(lengths)
    bins += list(range(1000, max_length + 500, 500))  # Extend bins every 500

    # Bin the lengths into the defined bins and create a DataFrame for tabulation
    df = pd.DataFrame({'Lengths': lengths})
    df['Groups'] = pd.cut(df['Lengths'], bins=bins, right=True)
    length_distribution = df['Groups'].value_counts().sort_index()
    return length_distribution


def plot_histogram(length_distribution, output_file):
    # Plot the histogram
    ax = length_distribution.plot(kind='bar', width=0.8)
    plt.xlabel('Read Length Groups (bp)')
    plt.ylabel('Frequency')
    plt.title('Read Length Distribution')

    # Customize x-axis labels
    x_labels = [str(int(interval.right)) for interval in length_distribution.index]
    ax.set_xticks(range(len(x_labels)))
    ax.set_xticklabels(x_labels, rotation=45, ha='right')

    # Show labels every 2000 bp
    ax.set_xticks([i for i, lbl in enumerate(x_labels) if int(lbl) % 2000 == 0])

    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Read length distribution and histogram plot from a fasta file.")
    parser.add_argument("fasta_file", help="Path to the input fasta file")
    parser.add_argument("output_file", help="Path to save the histogram plot")
    args = parser.parse_args()

    # Process read lengths
    lengths = read_lengths_distribution(args.fasta_file)
    length_distribution = group_read_lengths(lengths)

    # Print table
    print(length_distribution)

    # Plot histogram
    plot_histogram(length_distribution, args.output_file)


if __name__ == "__main__":
    main()