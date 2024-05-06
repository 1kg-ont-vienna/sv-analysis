import sys
import pandas as pd

qc_table_file = sys.argv[1]
phased_qc_table_file = sys.argv[2]
output = sys.argv[3]

qc_table = pd.read_csv(qc_table_file, sep='\t', header=0)
qc_phased_table = pd.read_csv(phased_qc_table_file, sep='\t', header=0)

qc_phased_table=qc_phased_table.rename(columns={'AC': 'PHASED_AC',
                        'HWE': 'PHASED_HWE',
                        'EXCHET': 'PHASED_EXCHET',
                        'IRS': 'PHASED_IRS',
                        'MAXR2': 'PHASED_MAXR2'})

df = pd.merge(qc_table, qc_phased_table, on='VARIANT_ID', how ="left")

df.to_csv(output, sep='\t', index=False)