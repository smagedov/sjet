import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

stability_list = []
stabdists_list = []

for i in range(200):
    stability_list.append(np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/stability/" + str(i) + "_stability.npy"))
    stabdists_list.append(np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/stabdists/" + str(i) + "_stabdists.npy"))

stability = np.concatenate(stability_list, axis=0)
stabdists = np.concatenate(stabdists_list, axis=0)

print(stabdists)
print(stability)

mask = (stability != 0.0) & (stability > -3000)
plt.hist2d(stabdists[mask], stability[mask], cmin=1, bins=50)
plt.colorbar(label='Count in bin')
plt.xlabel("Distance")
plt.ylabel("Stability")
plt.title("Stability as a function of distance")
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/instability_01.png")
plt.show()
