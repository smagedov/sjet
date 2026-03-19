import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

vecsumratios = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/vecsumratios.npy")
ratios = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/ratios.npy")
dists = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/dists.npy")
deltar = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/deltar.npy")
clusdis = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/clusdis.npy")

plt.hist2d(dists, ratios, cmin=1, bins=30)
plt.colorbar(label='Count in bin')
plt.xlim(0,5)
plt.xlabel("Distance")
plt.ylabel("Scalar Sum PT Ratio")
plt.title("Scalar Sum PT Ratio vs Distance")
plt.show()

plt.hist2d(dists, vecsumratios, cmin=1, bins=30)
plt.colorbar(label='Count in bin')
plt.xlim(0,5)
plt.xlabel("Distance")
plt.ylabel("Vec Sum PT Ratio")
plt.title("Vec Sum PT Ratio vs Distance")
plt.show()

plt.hist2d(dists, deltar, cmin=1, bins=30)
plt.colorbar(label='Count in bin')
plt.xlim(0,5)
plt.xlabel("Distance")
plt.ylabel("Delta R")
plt.title("Delta R vs Distance")
plt.show()

plt.hist2d(dists, clusdis, cmin=1, bins=30)
plt.colorbar(label='Count in bin')
plt.xlim(0,5)
plt.xlabel("Distance")
plt.ylabel("Delta R")
plt.title("Delta R (PT weight) vs Distance")
plt.show()

