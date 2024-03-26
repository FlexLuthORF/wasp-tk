import csv
import sys
import os

def process_csv_to_fasta(csv_file, output_path):
    fasta_data = {'V': [], 'V_gapped': [], 'D': [], 'J': []}
    fasta_filenames = {
        'V': os.path.join(output_path, 'human_ig_V_alleles.fasta'),
        'V_gapped': os.path.join(output_path, 'human_ig_V_alleles_gapped.fasta'),
        'D': os.path.join(output_path, 'human_ig_D_alleles.fasta'),
        'J': os.path.join(output_path, 'human_ig_J_alleles.fasta')
    }
    processed_alleles = set()

    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            allele = row['vdjbase_allele']
            if len(allele) >= 4 and allele not in processed_alleles:
                if len(allele) > 40:  # Check for allele length
                    allele = allele[:40]  # Trim allele

                category = allele[3]
                if category in ['V', 'D', 'J']:
                    region = row.get(f'{category}-REGION', '')
                    if region:
                        fasta_data[category].append(f'>{allele}\n{region}')
                        processed_alleles.add(allele)

                    if category == 'V':
                        region_gapped = row.get('V-REGION-GAPPED', '')
                        if region_gapped:
                            fasta_data['V_gapped'].append(f'>{allele}\n{region_gapped}')
            else:
                print(f"Skipped allele: {allele}")

    for key, data in fasta_data.items():
        if data:
            with open(fasta_filenames[key], 'w') as fasta_file:
                fasta_file.write('\n'.join(data))
        else:
            print(f"No data for category {key}, no file created.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <csv_file> <output_path>")
    else:
        process_csv_to_fasta(sys.argv[1], sys.argv[2])
