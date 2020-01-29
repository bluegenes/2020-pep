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
            #samples.drop(columns=["unit"], inplace=True)
            samples.set_index(["sample"], drop=False, inplace=True)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in samples_file:
        try:
            samples = pd.read_excel(samples_file, dtype=str, sep=separator)
            #samples['sample'] = samples["sample"].map(str) + '_' + samples["unit"].map(str)
            #samples.drop(columns=["unit"], inplace=True)
            samples.set_index(["sample"], drop=False, inplace=True)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    return samples

def write_yaml(yamlD, paramsfile):
    with open(paramsfile, 'w') as out:
        yaml.dump(yamlD, stream=out, indent=2, default_flow_style=False)
