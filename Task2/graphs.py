from cmath import log
from email.quoprimime import header_decode
import matplotlib.pyplot as plt
import pandas as pd

def get_var_name(var, var_names = locals()):
    return [var_name for var_name in var_names if id(var) == id(var_names[var_name])]

data = pd.read_csv('appr.csv')

headers = data.columns
print(headers)
base = [data.get(headers[i]) for i in range(3)]

for i in range(3, len(headers)):
    plt.figure(headers[i])
    plt.xlabel('x')
    plt.ylabel('y')
    for j in range(1, 3):
        plt.plot(base[0], base[j], label=headers[j])
    plt.plot(base[0], data.get(headers[i]), label = headers[i])

    plt.legend()
    plt.show()

