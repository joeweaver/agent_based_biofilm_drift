import os
import sys
import pandas as pd

sweep_dir = sys.argv[1]
csv_files = sys.argv[2:]

frames = []

for csv_file in csv_files:
    frames.append(pd.read_csv(csv_file.strip()))

final = pd.concat(frames, ignore_index = True)
final.to_csv(os.path.join(sweep_dir,'sweep_colony_outcomes.csv'), index = False)
