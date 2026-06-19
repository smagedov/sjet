import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

efp_list = []
efpdist_list = []

for i in range(50):
    efp_list.append(np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/efp/" + str(i) + "_efp.npy"))
    efpdist_list.append(np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/efpdist/" + str(i) + "_efpdist.npy"))

efp = np.concatenate(efp_list, axis=0)
efpdist = np.concatenate(efpdist_list, axis=0)

plt.hist2d(efpdist, efp, cmin=1, bins=50)
plt.colorbar(label='Count in bin')
plt.xlabel("Distance")
plt.ylabel("EFP")
plt.title("EFP as a function of distance")
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/EFP.png")
plt.show()

