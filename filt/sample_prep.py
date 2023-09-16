import pandas as pd
import numpy as np

tf_file = '~/03-Project/_Research_Perf_Var/10-Accuracy/experiments/8-FilTStencil/amd9654-2.4G/TL_CALC_W-C2048-PAPI-sample.csv'
df = pd.read_csv(tf_file, header=None)
ovh = df[2] - df[1] * 2 / 2.4
ovh.astype('int').to_csv(path_or_buf='tf.csv', index=False, header=False)
