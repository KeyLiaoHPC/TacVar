import pandas as pd
import numpy as np

met_file = '~/03-Project/_Research_Perf_Var/10-Accuracy/experiments/8-FilTStencil/amd9654-2.4G/TL_CALC_W-C2048-PAPI.csv'
df = pd.read_csv(met_file, header=None)[1]
df.to_csv(path_or_buf='met.csv', index=False, header=False)
