import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

stability_list = []
stabdists_list = []
for i in range(27):
    stability_list.append(np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/" + str(i) + "_stability.npy"))
    stabdists_list.append(np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/" + str(i) + "_stabdists.npy"))
stability = np.concatenate(stability_list, axis=0)
stabdists = np.concatenate(stabdists_list, axis=0)

print(stability)

plt.hist2d(stabdists, stability, cmin=1, bins=50)
plt.colorbar(label='Count in bin')
plt.xlabel("Distance")
plt.ylabel("Stability")
plt.title("Stability as a function of distance")
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/stability.png")
plt.show()

vecsumratios_list = []
ratios_list = []
dists_list = []
deltar_list = []
clusdis_list = []
masses_list = []
for i in range(1000):
    tmp_vec = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/" + str(i) + "_vecsumratios.npy")
    tmp_rat = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/" + str(i) + "_ratios.npy")
    tmp_dis = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/" + str(i) + "_dists.npy")
    tmp_del = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/" + str(i) + "_deltar.npy")
    tmp_clu = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/" + str(i) + "_clusdis.npy")
    tmp_mas = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/" + str(i) + "_masses.npy")
    vecsumratios_list.append(tmp_vec)
    ratios_list.append(tmp_rat)
    dists_list.append(tmp_dis)
    deltar_list.append(tmp_del)
    clusdis_list.append(tmp_clu)
    masses_list.append(tmp_mas)
vecsumratios = np.concatenate(vecsumratios_list, axis=0)
ratios = np.concatenate(ratios_list, axis=0)
dists = np.concatenate(dists_list, axis=0)
deltar = np.concatenate(deltar_list, axis=0)
clusdis = np.concatenate(clusdis_list, axis=0)
masses = np.concatenate(masses_list, axis=0)

plt.hist2d(dists, ratios, cmin=1, bins=100)
plt.colorbar(label='Count in bin')
plt.xlabel("Distance")
plt.ylabel("Scalar Sum PT Ratio")
plt.title("Scalar Sum PT Ratio vs Distance")
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/clusqual_scalarsumptratio.png")
plt.show()

mask = (vecsumratios >= 0.3)
plt.hist2d(dists[mask], vecsumratios[mask], cmin=1, bins=100)
plt.colorbar(label='Count in bin')
plt.xlabel("Distance")
plt.ylabel("Vec Sum PT Ratio")
plt.title("Vec Sum PT Ratio vs Distance")
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/clusqual_vecsumptratio.png")
plt.show()

mask = (deltar < 0.5)
plt.hist2d(dists[mask], deltar[mask], cmin=1, bins=100)
plt.colorbar(label='Count in bin')
plt.xlabel("Distance")
plt.ylabel("Delta R")
plt.title("Delta R vs Distance")
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/clusqual_deltar.png")
plt.show()

plt.hist2d(dists, clusdis, cmin=1, bins=100)
plt.colorbar(label='Count in bin')
plt.xlabel("Distance")
plt.ylabel("Delta R")
plt.title("Delta R (PT weight) vs Distance")
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/clusqual_deltarptratio.png")
plt.show()

mask = (masses >= 0.0) & (masses <= 300.0)
plt.hist2d(dists[mask], masses[mask], cmin=1, bins=100)
plt.colorbar(label='Count in bin')
plt.xlabel("Distance")
plt.ylabel("Mass")
plt.title("Mass vs Distance")
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/clusqual_mass.png")
plt.show()
