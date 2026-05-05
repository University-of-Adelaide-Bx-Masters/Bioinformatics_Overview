!#bin/bash


python3 ./vcf_validator.py

awk -f extract_variants.awk ../4_refs/DPYD_variants_genome_location.csv ../1_vcfs/*.vcf > ../3_reports/variant_info.txt

python3 ./generate_DPYD_report.py

