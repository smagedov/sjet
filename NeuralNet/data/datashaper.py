import numpy as np

dists = np.load("conc_dists.npy")
deltar = np.load("conc_deltar.npy")
massd = np.load("conc_massd.npy")
massp1 = np.load("conc_massp1.npy")
massp2 = np.load("conc_massp2.npy")
ptd = np.load("conc_ptd.npy")
pt1 = np.load("conc_pt1.npy")
pt2 = np.load("conc_pt2.npy")
eta1 = np.load("conc_eta1.npy")
eta2 = np.load("conc_eta2.npy")
phi1 = np.load("conc_phi1.npy")
phi2 = np.load("conc_phi2.npy")

features = np.vstack((ptd, massd, massp1, massp2, pt1, pt2, eta1, eta2)).T

print(features.shape)

np.save("conc_features.npy", features)
