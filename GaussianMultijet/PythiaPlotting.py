import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

vecsumratios = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/vecsumratios.npy")
ratios = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/ratios.npy")
dists = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/dists.npy")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

ax1.scatter(dists, ratios, color='blue', label='Scalar Sum PT Ratio')
ax2.scatter(dists, vecsumratios, color='red', label='Vec Sum PT Ratio')

ax1.set_xlabel("Distance")
ax1.set_ylabel("PT Ratio")
ax1.set_title("PT Sum Ratio vs Distance")

ax2.set_xlabel("Distance")
ax2.set_ylabel("PT Ratio")
ax2.set_title("PT Sum Ratio vs Distance")

plt.tight_layout()
plt.show()
