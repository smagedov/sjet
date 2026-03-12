import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ratios = np.load("/home/smagedov/local/src/sjet-master/GaussianMultijet/npyarrays/ratios.npy")
dists = np.load("/home/smagedov/local/src/sjet-master/GaussianMultijet/npyarrays/dists.npy")

plt.scatter(dists, ratios)

plt.xlabel("Distance")
plt.ylabel("Ratio")
plt.title("Ratio vs Distance")

plt.show()
