import pandas as pd
import numpy as np
from pathlib import Path

per_sample = pd.read_csv("/exports/sasc/project-247-fetus_methylation/analysis"
                "/04_methAnalysisHisat2/PQ5/methylator/perSample/charlotte_chr11_1999844.Total.tsv", sep="\t")

r = pd.pivot_table(per_sample, values="counts", index="methStatesCount",
                 aggfunc=
np.sum)
print(r)