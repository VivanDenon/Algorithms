from PIL import Image
import numpy as np
import random as rd
import string as str

img = Image.open("1.png")
m = np.asarray(img)
r = [[False for i in row] for row in m]

for row in enumerate(m):
    for p in enumerate(row[1]):
        if(p[1][3] > 100):
            r [row[0]][p[0]] = True

for i in r:
    s = ''
    for j in i:
        if not j:
            s += rd.choice(str.ascii_letters) + rd.choice(str.digits)
        else:
            s += '  '
    print(s)
