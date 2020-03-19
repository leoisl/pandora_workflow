import glob
import re
import sys
import pickle

out_files = glob.glob("*.out")

gene_to_RAM_usage = {}
for file in out_files:
    print(f"Procesing {file} ...", file=sys.stderr)
    with open(file) as fh:
        lines = fh.readlines()
        lines = "".join(lines)
    
    if "add_denovo_paths" in lines and "Successfully completed." in lines:
         memory_as_str = re.search(r"Max Memory :\s+(.+?)\s+", lines).group(1)
         memory_is_dash_or_int = memory_as_str=="-" or str.isdigit(memory_as_str)
         assert memory_is_dash_or_int, f"Could not get memory argument: {memory_as_str}."

         if memory_as_str != "-":
              memory = int(memory_as_str)
              gene = re.search(r"gene=(.*?)>", lines).group(1)
              gene_to_RAM_usage[gene]=memory


with open('gene_to_RAM_usage.pickle', 'wb') as handle:
    pickle.dump(gene_to_RAM_usage, handle, protocol=pickle.HIGHEST_PROTOCOL)
