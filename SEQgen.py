import random

str = "ACTG"
out = ''

for i in range(0,100):
    out += random.choice(str)

print(out)

