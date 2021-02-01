import sys
sys.path.insert(1, "/".join(sys.path[0].split("/")[:-1]))

from support.progress_info import progress_bar

import pandas as pd
import numpy as np

from numpy import transpose as t

data = pd.read_csv("energy_data.csv")

n = len(data.index)

# i + 1 protoze index jde v Pythonu od nuly
year = [i + 1 for i in data.index]

q = [val for val in data.log_quantity]
p = [val for val in data.log_price]
i = [val for val in data.income]
t = [val for val in data.date]

y = np.matrix(q)
x = np.matrix(([[1]*n], year, p, np.log(i)))

print(x)
