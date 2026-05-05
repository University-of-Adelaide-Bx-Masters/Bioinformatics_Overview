import os
import argparse 

def main(input_vcf):
    ####### 1. Checking if an input vcf file exists ###########
    vcf_file = input_vcf

    # check whether the input file exist
    if not os.path.exists(vcf_file):
        #Error is output below. Change it to make it more informative!
        raise Exception("Unknown Error")
    #else:
        #print(f'Successfully found input path for vcf!')


    ####### 2. Checking if the input vcf file has the right extension ###########

    if not vcf_file.lower().endswith(".vcf"):
        raise ValueError("Unknown Error number 2")
    #else:
        #print(f'Checked vcf ends with correct extension - passed')


    ####### 3. Checking if the format of the vcf file is correct ###########

    required_cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

    main_vcf_header = False
    with open(vcf_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            print(line)
        # Metadata lines
            if line.startswith("##"):
                if main_vcf_header:
                    # Should not be metadata after header
                    raise Exception("Header with columns should not be followed by metadata")
                continue

            # Column header line
            if line.startswith("#CHROM"):
                if main_vcf_header:
                    # Multiple headers → invalid
                    raise Exception("Too many headers are present in vcf file")
                main_vcf_header = True
                cols = line.split("\t")

            # Must have at least the 8 required columns
                if len(cols) < 8:
                    raise Exception("A column is missing from the vcf main header")

                # Check exact order of required columns
                if cols[:8] != required_cols:
                    raise Exception("The order of columns in the vcf is not correct")
                
                    # After this, we expect only data lines
                continue

                # Any other line before header is invalid
            if not main_vcf_header:
                raise Exception("Cannot have data lines preceding header")
            break
    #print("Checked header with columns comes after metadata - passed")
    #print("Checked correct number of headers - passed")
    #print("Checked all columns present - passed")
    #print("Checked order of columns meets vcfs specs - passed")
    #print("Checked data lines come after headers - passed")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_vcf")
    args = parser.parse_args()
    main(args.input_vcf)