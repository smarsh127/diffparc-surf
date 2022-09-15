from glob import glob
import os

in_dir=snakemake.input.conn_dir
header_line=snakemake.params.header_line
out_csv=snakemake.output.conn_csv


csv_files = glob(os.path.join(in_dir,'*.csv'))
all_lines = list()

all_lines.append(header_line+'\n')

for csv in csv_files:
    with open(csv,'r') as f:
        skiphdr = f.readline()
        all_lines.append(f.readline())

with open(out_csv,'w') as f:
    f.writelines(all_lines)


