import sys
sys.path.insert(1, "/".join(sys.path[0].split("/")[:-1]))
import support

from pandas import read_csv
data = read_csv 
