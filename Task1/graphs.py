from cmath import log
import matplotlib.pyplot as plt
import pandas as pd

def get_var_name(var, var_names = locals()):
    return [var_name for var_name in var_names if id(var) == id(var_names[var_name])]

data = [ pd.read_csv(f"uint{2**i}.csv") for i in range(4, 7)]

linear = lambda n : n
logarithmic = lambda n : log(n)
nlogn = lambda n : n * log(n)
quadratic = lambda n : n ** 2
cubic = lambda n : n ** 2

func = [linear, linear, linear, linear, linear, linear, quadratic, nlogn, nlogn, cubic] 


headers = data[0].columns
print(headers)

for head, f in zip(headers, func):
    plt.figure(head)
    plt.xlabel("n")
    plt.ylabel("ms")
    for i in range(3):
        plt.plot(data[i].get(head), label=(f"uint{2**(4 + i)}"))
    plt.legend()
    plt.show()

    plt.plot([f(i) for i in range(1, 100)], label=get_var_name(f)[0]) 
    plt.legend()
    plt.show()

