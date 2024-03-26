import sys
import os
import shutil
from cyvcf2 import VCF, Writer

def read_data(filename):
    with open(filename, 'r') as file:
        data = [line.strip().split('\t') for line in file]
    return data

def read_bed_file(bed_filename):
    bed_entries = {}
    with open(bed_filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            bed_entries[parts[3]] = parts[:3]  # Store chromosome, start, and end
    return bed_entries

def create_sample_bed_file(hemi_regions, bed_entries, sample_id, outd):
    sample_bed_filename = os.path.join(outd, '{}_hemi.bed'.format(sample_id))
    with open(sample_bed_filename, 'w') as bed_file:
        for region in hemi_regions:
            bed_file.write('\t'.join(bed_entries[region]) + '\t' + region + '\n')
    return sample_bed_filename

def process_samples(data, bed_entries, input_vcf, outd, sample_id):
    hemi_regions = []
    for row in data:
        if row[0] == sample_id:
            sv_name = row[1]
            genotype = row[2]
            if genotype in ['0/1'] and sv_name in bed_entries:
                hemi_regions.append(sv_name)

    if hemi_regions:
        output_vcf_dir = os.path.join(outd, "change_to_hemi")
        os.makedirs(output_vcf_dir, exist_ok=True)
        output_vcf = os.path.join(output_vcf_dir, f'{sample_id}_hemi.vcf')
        sample_bed_filename = create_sample_bed_file(hemi_regions, bed_entries, sample_id, outd)
        change_genotypes(input_vcf, sample_bed_filename, output_vcf)
        return output_vcf
    else:
        return input_vcf

def read_bedfile(bedfile):
    coords = []
    with open(bedfile, 'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            coords.append([chrom, start, end])
    return coords

def change_genotypes(input_vcf, bedfile, output_path):
    vcf_reader = VCF(input_vcf, strict_gt=True)

    # Create a new VCF Writer using the input VCF as a template
    vcf_writer = Writer(output_path, vcf_reader)

    coords = read_bedfile(bedfile)

    for record in vcf_reader:
        for index, sample in enumerate(vcf_reader.samples):
            change_snp = False
            for chrom, start, end in coords:
                if record.CHROM != chrom:
                    continue
                if record.POS < start:
                    continue
                if record.POS > end:
                    continue
                change_snp = True
                break

            if change_snp:
                current_genotype = record.genotypes[index]
                if current_genotype is None or -1 in current_genotype:
                    record.genotypes[index] = [-1, -1, False]  # Set genotype to unknown
                elif 1 in current_genotype:
                    record.genotypes[index] = [-1, 1, False]  # Set genotype to ./1
                else:
                    record.genotypes[index] = [-1, 0, False]  # Set genotype to ./0
                
                # Reassign the genotypes field
                record.genotypes = record.genotypes

        vcf_writer.write_record(record)

    vcf_writer.close()
    vcf_reader.close()

def main(input_vcf, genotypes_filename, bed_filename, sample_id):
    outd = os.path.dirname(input_vcf)
    data = read_data(genotypes_filename)
    bed_entries = read_bed_file(bed_filename)
    output_vcf = process_samples(data, bed_entries, input_vcf, outd, sample_id)
    return output_vcf

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: python process_samples.py <input_vcf> <genotypes_filename> <bed_filename> <sample_id>")
        sys.exit(1)
    output_vcf = main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    print(output_vcf)