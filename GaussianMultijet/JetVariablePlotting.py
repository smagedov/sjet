import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pts_list = []
ms_list = []
phis_list = []
etas_list = []
dists_list = []

for i in range(100):
    flatPts = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/flatPts/" + str(i) + "_flatPts.npy")
    flatMs = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/flatMs/" + str(i) + "_flatMs.npy")
    flatEtas = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/flatEtas/" + str(i) + "_flatEtas.npy")
    flatPhis = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/flatPhis/" + str(i) + "_flatPhis.npy")
    flatDists = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/flatDists/" + str(i) + "_flatDists.npy")
    offsets = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/offsets/" + str(i) + "_offsets.npy")
    pts = []
    ms = []
    etas = []
    phis = []
    dists = []
    for i in range(len(offsets)-1):
        pts.append(flatPts[offsets[i]:offsets[i+1]])
        ms.append(flatMs[offsets[i]:offsets[i+1]])
        etas.append(flatEtas[offsets[i]:offsets[i+1]])
        phis.append(flatPhis[offsets[i]:offsets[i+1]])
        dists.append(flatDists[offsets[i]:offsets[i+1]])
    pts_list.append(pts)
    ms_list.append(ms)
    etas_list.append(etas)
    phis_list.append(phis)
    dists_list.append(dists)

# Loop over every jet in the first event
for jet in range(len(dists_list[10])):

    d = np.asarray(dists_list[0][jet])
    eta = np.asarray(etas_list[0][jet])
    phi = np.asarray(phis_list[0][jet])

    mask = d > 0.0

    d = d[mask]
    eta = eta[mask]
    phi = phi[mask]

    # Sort by distance
    order = np.argsort(d)
    d = d[order]
    eta = eta[order]
    phi = phi[order]

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))

    # eta-phi trajectory
    ax[0].scatter(eta, phi)
    ax[0].set_xlabel(r"$\eta$")
    ax[0].set_ylabel(r"$\phi$")
    ax[0].set_title(f"Event 0 Jet {jet}: Trajectory")
    ax[0].grid(True)

    # coordinate evolution
    ax[1].plot(d, eta, "-o", label=r"$\eta$")
    ax[1].plot(d, phi, "-s", label=r"$\phi$")
    ax[1].set_xlabel("Distance")
    ax[1].set_ylabel("Coordinate")
    ax[1].set_title(f"Event 0 Jet {jet}: History")
    ax[1].legend()
    ax[1].grid(True)

    plt.tight_layout()
    plt.show()

    plt.close()

# Pt Plotting
plt.figure(figsize=(8, 6))

all_d = []
all_p = []

for event_dists, event_pts in zip(dists_list, pts_list):
    for d, p in zip(event_dists, event_pts):
        d = np.asarray(d)
        p = np.asarray(p)

        mask = (d > 0) & (d < 3)
        all_d.extend(d[mask])
        all_p.extend(p[mask])

all_d = np.asarray(all_d)
all_p = np.asarray(all_p)

nbins = 50
bins = np.linspace(all_d.min(), all_d.max(), nbins + 1)
centers = 0.5 * (bins[:-1] + bins[1:])

mean_pt = np.full(nbins, np.nan)

for i in range(nbins):
    mask = (all_d >= bins[i]) & (all_d < bins[i+1])
    if np.any(mask):
        mean_pt[i] = all_p[mask].mean()

plt.plot(centers, mean_pt, '-o')
plt.xlabel("Distance")
plt.ylabel("Average $p_T$")
plt.title("Average Jet History")
plt.grid(True)
plt.show()

#for event_dists, event_pts in zip(dists_list, pts_list):
#    for d, p in zip(event_dists, event_pts):
#        d = np.asarray(d)
#        p = np.asarray(p)
#
#        mask = d > 0.0
#        d = d[mask]
#        p = p[mask]
#
#        order = np.argsort(d)
#        d = d[order]
#        p = p[order]
#
#        plt.plot(d, p, '-o')
#
#plt.title("Jet History Pt")
#plt.xlabel("dist")
#plt.ylabel("pt")
#plt.grid(True)
#plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/jetHistory/jetHistory_pt.png")
#plt.show()


# Mass Plotting

all_d = []
all_m = []

for event_dists, event_ms in zip(dists_list, ms_list):
    for d, m in zip(event_dists, event_ms):
        d = np.asarray(d)
        m = np.asarray(m)

        mask = (d > 0) & (d < 3)
        all_d.extend(d[mask])
        all_m.extend(m[mask])

all_d = np.asarray(all_d)
all_m = np.asarray(all_m)

nbins = 50
bins = np.linspace(all_d.min(), all_d.max(), nbins + 1)
centers = 0.5 * (bins[:-1] + bins[1:])

mean_m = np.full(nbins, np.nan)

for i in range(nbins):
    mask = (all_d >= bins[i]) & (all_d < bins[i+1])
    if np.any(mask):
        mean_m[i] = all_m[mask].mean()

plt.plot(centers, mean_m, '-o')
plt.xlabel("Distance")
plt.ylabel("Average $m$")
plt.title("Average Jet History")
plt.grid(True)
plt.show()

#for event_dists, event_ms in zip(dists_list, ms_list):
#    for d, m in zip(event_dists, event_ms):
#        d = np.asarray(d)
#        m = np.asarray(m)
#
#        mask = d > 0.0
#        d = d[mask]
#        m = m[mask]
#
#        order = np.argsort(d)
#        d = d[order]
#        m = m[order]
#
#        plt.plot(d, m, '-o')
#
#plt.title("Jet History Mass")
#plt.xlabel("dist")
#plt.ylabel("mass")
#plt.grid(True)
#plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/jetHistory/jetHistory_m.png")
#plt.show()

# Eta-Phi Plotting

all_etas = []
all_phis = []
all_d = []

for event_dists, event_etas, event_phis in zip(dists_list, etas_list, phis_list):
    for d, etas, phis in zip(event_dists, event_etas, event_phis):
        etas = np.asarray(etas)
        phis = np.asarray(phis)

        mask = (d > 0) & (d < 3)
        all_etas.extend(etas[mask])
        all_phis.extend(phis[mask])

all_etas = np.asarray(all_etas)
all_phis = np.asarray(all_phis)

nbins = 50
bins = np.linspace(all_etas.min(), all_etas.max(), nbins + 1)
centers = 0.5 * (bins[:-1] + bins[1:])

mean_phis = np.full(nbins, np.nan)

for i in range(nbins):
    mask = (all_etas >= bins[i]) & (all_etas < bins[i+1])
    if np.any(mask):
        mean_phis[i] = all_phis[mask].mean()

plt.plot(centers, mean_phis, '-o')
plt.xlabel("Eta")
plt.ylabel("Average Phi")
plt.title("Jet Location History")
plt.grid(True)
plt.show()

# Eta-Phi Plotting

#for event_dists, event_etas, event_phis in zip(dists_list, etas_list, phis_list):
#    for d, eta, phi in zip(event_dists, event_etas, event_phis):
#        d = np.asarray(d)
#        eta = np.asarray(eta)
#        phi = np.asarray(phi)
#        mask = d > 0.0
#        plt.scatter(eta[mask], phi[mask], marker='o')
#
#plt.title("Jet Location History")
#plt.xlabel("eta")
#plt.ylabel("Phi")
#plt.grid(True)
#plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/jetHistory/jetHistory_location.png")
#plt.show()
