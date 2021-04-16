# Normalize-RNA-Counts

This program has been designed to take a table of counts per gene (or feature) from RNA-seq and normalize the counts.
If a genome features file (-g option) containing information about transcript lengths is provided, the default
normalization is read counts per kilobase per million/fragments per kilobase per million (RPKM/FPKM), but can be
changed to transcripts per million (TPM) by using the -t or --tpm option. If no genome file is provided, counts will
be normalized to read counts per million (RCPM). Use the -r option to provide a list of features that you do not want
to be counted for normalization (excluding lines starting with either 'N_' or '__' added by HTSEq-count or STAR).

Example usage:
```
python normalize_rna_counts.py -c all_samples_counts.csv -g GENCODE_Longest_CDS.csv -t -o All_samples -r remove.txt
```

optional arguments:

  -h, --help
                        show this help message and exit

  -c COUNTS_FILE, --counts_file COUNTS_FILE

                        This is a file that contains counts per gene (or feature) for all samples. This must be in csv
                        format with the first row being sample names and the first column being feature names.
  -g GENOME_FEATURES, --genome_features GENOME_FEATURES

                        This is a file that contains information about the features contained in your dataset. The
                        column containing the feature name that matches the first column in the --counts_file must be
                        labeled as 'gene_id' and the column containing the length of the feature must be labeled as
                        'CDS Length'. Instructions on how to generate the genome_features file from a GTF file can be
                        found here: https://github.com/davidwsant/parse_gtf

  -t, --tpm

                        Specify this option if you would like to scale using TPM rather than RPKM/FPKM.

  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX

                        You can use this option to specify a prefix that you would like added to the normalized output
                        file. By default, the prefix will 'Normalized'.

  -r REMOVE, --remove REMOVE

                        Use this option to supply a file with specific feature names that you want to be excluded from
                        the total cell reads for normalization. Do not include in this file features added by HTSeq-
                        count that begin with '__' or features added by STAR that begin with 'N_' as they will already
                        be removed. Features that you may want to include are those that were identified as outliers
                        during differential analysis by a method such as Cook's cutoff.
