import os
from vcf_processing import change_genotypes, read_bedfile

def main():
    # Input VCF file path
    input_vcf = "/home/zmvanw01/projects/12-sample-comparison/hifiasm-path/geno_analysis/per_samp/2201410002/2201410002.vcf"

    # Create a fake bed file with both SVs
    fake_bed_file = "fake_sample.bed"
    with open(fake_bed_file, "w") as file:
        file.write("chr22\t22811963\t22821062\tIGLV5-39\n")
        file.write("chr2\t90225183\t90249908\tIGKV1-NL1\n")

    # Output VCF file path
    output_vcf = "test_double.vcf"

    # Call the change_genotypes function with the fake bed file
    change_genotypes(input_vcf, fake_bed_file, output_vcf)

    print(f"Modified VCF file saved as: {output_vcf}")

if __name__ == "__main__":
    main()