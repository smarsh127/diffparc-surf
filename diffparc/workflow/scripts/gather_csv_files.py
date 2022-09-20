from glob import glob
import os
import numpy as np

in_dir = snakemake.input.conn_dir
header_line = snakemake.params.header_line
out_csv = snakemake.output.conn_csv


csv_files = sorted(glob(os.path.join(in_dir, "*.csv")))
all_lines = list()


#make this function backwards-compatible -- newer versions of mrtrix use csv and have command_history line
# older version is space-delimited, no header line
with open(csv_files[0],"r") as f:
    first_line = f.readline()
    print(first_line)
    if 'command_history' in first_line:
        load_opts={'delimiter': ',', 'skiprows':1}
    else:
        load_opts={}

all_lines.append(header_line + "\n")

for csv in csv_files:
    row_data = np.loadtxt(csv,**load_opts,dtype='int16')
    all_lines.append(','.join([f'{r}' for r in row_data])+'\n' )


with open(out_csv, "w") as f:
    f.writelines(all_lines)
