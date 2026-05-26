import numpy as np

dists = np.load("conc_dists.npy")
ratios = np.load("conc_ratios.npy")
masses = np.load("conc_masses.npy")
deltar = np.load("conc_deltar.npy")

features = np.vstack((dists, ratios, masses, deltar)).T

print(features.shape)

np.save("conc_features.npy", features)
