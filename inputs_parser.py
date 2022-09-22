#Parse inputs from inputs file
import sys
import os

file = open("inputs", "r")
lines = file.readlines()
file.close()

#Parse inputs
v = float(lines[0])
a = float(lines[1])
b = float(lines[2])
e1 = float(lines[3])
e2 = float(lines[4])
eps11 = float(lines[5])
eps22 = float(lines[6])
eps12 = float(lines[7])

c11m = e1*(1-v)/((1+v)*(1-2*v))
c12m = e1*v/((1+v)*(1-2*v))
c44m = e1/(2*(1+v))

c11p = e2*(1-v)/((1+v)*(1-2*v))
c12p = e2*v/((1+v)*(1-2*v))
c44p = e2/(2*(1+v))

#Write elastic inputs to input/elastic_constants by first reading the file
file = open("input/elastic_constants", "r")
lines = file.readlines()
file.close()

print(lines)
#Write elastic inputs to input/elastic_constants
lines[0] = str(eps11)+'\n'
lines[1] = str(eps22)+'\n'
lines[2] = str(eps12)+'\n'

lines[8] = str(c11m)+'\n'
lines[9] = str(c12m)+'\n'
lines[10] =str(c44m)+'\n'
lines[11] =str(c11p)+'\n'
lines[12] =str(c12p)+'\n'
lines[13] =str(c44p)+'\n'

print(lines)

#Write elastic inputs to input/elastic_constants
file = open("input/elastic_constants", "w")
file.writelines(lines)
file.close()