import sys
import pandas as pd

def read_samples(samples_file, build_sra_links = False):
    if '.tsv' in samples_file or '.csv' in samples_file:
        separator = '\t'
        if '.csv' in samples_file:
            separator = ','
        try:
            samples = pd.read_csv(samples_file, dtype=str, sep=separator)
            #samples['sample'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
            samples.drop(columns=["unit"], inplace=True)
            samples.set_index(["sample"], drop=False, inplace=True)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in samples_file:
        try:
            samples = pd.read_excel(samples_file, dtype=str, sep=separator)
            #samples['sample'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
            samples.drop(columns=["unit"], inplace=True)
            samples.set_index(["sample"], drop=False, inplace=True)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    # build sra links if necessary (reads)
    if build_sra_links:
        base_link = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
        if "SRR" in samples.columns and "LibraryLayout" in samples.columns:
            samples['fq1'] = samples['SRR'].apply(lambda x : base_link + x[0:6] + '/00' + x[-1] + '/' + x + '/' + x + '_1.fastq.gz')
            samples['fq2'] = df.apply(lambda row : build_fq2(row), axis=1)
        else:
            sys.stderr.write(f"\n\tError: To build SRA links, 'SRR' column must exist in samples file \n\n")
    return samples
