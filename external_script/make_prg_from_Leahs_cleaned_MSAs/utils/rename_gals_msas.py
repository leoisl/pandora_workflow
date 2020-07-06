# configs
input_dir = "/hps/nobackup2/iqbal/leandro/gal_msa/10perc"
output_dir = "/hps/nobackup2/iqbal/leandro/gal_msa/10perc_renamed"

import os
from shutil import copyfile

os.makedirs(output_dir, exist_ok=True)
for file_index, file in enumerate(os.listdir(input_dir)):
    source_file_full_path = os.path.join(input_dir, file)
    dest_file_full_path = os.path.join(output_dir, f"group_{file_index}.fa")
    copyfile(source_file_full_path, dest_file_full_path)