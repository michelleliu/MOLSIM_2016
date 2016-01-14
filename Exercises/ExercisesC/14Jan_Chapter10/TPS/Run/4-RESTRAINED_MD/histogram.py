#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt("committor.dat",comments="#")

bins = np.linspace(0.35, 0.65, 10)
n,bins,patches = plt.hist(data[:,2], bins, normed = True, facecolor='green')

plt.xlabel('Commitor')
plt.ylabel('Fraction of frames')
plt.xlim((0.35,0.65))

plt.show()
plt.savefig("histogram.png")
