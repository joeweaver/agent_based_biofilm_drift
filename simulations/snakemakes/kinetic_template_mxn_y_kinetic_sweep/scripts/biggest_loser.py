import pandas as pd
import sys
import re
import os

results_dir = sys.argv[1] 
results_file = sys.argv[2]
lineages = int(sys.argv[3])
print(lineages)
results = pd.read_csv(os.path.join(results_dir, results_file))
final_result = results.iloc[-1:]

het_cols = [f'het{n}_rv' for n in range(1, lineages+1)]
min_het_idx = final_result[het_cols].idxmin(axis=1)

min_het_colname = min_het_idx.iloc[0]
min_het_number = re.search("het([0-9]+)_rv", min_het_colname)[1]


with open(os.path.join(results_dir, 'biggest_loser.txt'), 'w') as f:
    f.write(min_het_number)
