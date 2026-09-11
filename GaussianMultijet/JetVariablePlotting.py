import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
import re

pts_list = []
ms_list = []
phis_list = []
etas_list = []
dists_list = []
gamma_list = []

for i in range(5):
    flatPts = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/oraclejets/0/flatPts/" + str(i) + "_flatPts.npy")
    flatMs = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/oraclejets/0/flatMs/" + str(i) + "_flatMs.npy")
    flatEtas = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/oraclejets/0/flatEtas/" + str(i) + "_flatEtas.npy")
    flatPhis = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/oraclejets/0/flatPhis/" + str(i) + "_flatPhis.npy")
    flatDists = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/oraclejets/0/flatDists/" + str(i) + "_flatDists.npy")
    flatGamma = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/oraclejets/0/flatGamma/" + str(i) + "_flatGamma.npy")
    offsets = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/jethistory/oraclejets/0/offsets/" + str(i) + "_offsets.npy")
    pts = []
    ms = []
    etas = []
    phis = []
    dists = []
    gamma = []
    for i in range(len(offsets)-1):
        pts.append(flatPts[offsets[i]:offsets[i+1]])
        ms.append(flatMs[offsets[i]:offsets[i+1]])
        etas.append(flatEtas[offsets[i]:offsets[i+1]])
        phis.append(flatPhis[offsets[i]:offsets[i+1]])
        dists.append(flatDists[offsets[i]:offsets[i+1]])
        gamma.append(flatGamma[offsets[i]:offsets[i+1]])
    pts_list.append(pts)
    ms_list.append(ms)
    etas_list.append(etas)
    phis_list.append(phis)
    dists_list.append(dists)
    gamma_list.append(gamma)

# Loop over every jet in the first event
for jet in range(len(dists_list[0])):

    d = np.asarray(dists_list[0][jet])
    g = np.asarray(gamma_list[0][jet])
    eta = np.asarray(etas_list[0][jet])
    phi = np.asarray(phis_list[0][jet])
    m = np.asarray(ms_list[0][jet])
    pt = np.asarray(pts_list[0][jet])

    mask = (d > 0.0) & (d < 4.0)

    d = d[mask]
    g = g[mask]
    eta = eta[mask]
    phi = phi[mask]
    m = m[mask]
    pt = pt[mask]

    # Sort by distance
    order = np.argsort(d)
    d = d[order]
    g = g[order]
    eta = eta[order]
    phi = phi[order]
    m = m[order]
    pt = pt[order]

    fig, ax = plt.subplots(2, 3, figsize=(12, 8))

    # eta-phi trajectory
    scatter = ax[0][0].scatter(
    eta,
    phi,
    c=d,
    cmap="viridis",
    s=40)
    fig.colorbar(scatter, ax=ax[0][0], label="Distance")
    ax[0][0].set_xlabel(r"$\eta$")
    ax[0][0].set_ylabel(r"$\phi$")
    ax[0][0].set_title(f"Event 0 Jet {jet}: Trajectory")
    ax[0][0].grid(True)

    # coordinate evolution
    ax[0][1].scatter(g, eta, label=r"$\eta$")
    ax[0][1].set_xlabel("Gamma")
    ax[0][1].set_ylabel("Eta")
    ax[0][1].set_title(f"Event 0 Jet {jet}: Eta History")
    ax[0][1].legend()
    ax[0][1].grid(True)

    ax[0][2].scatter(g, phi, label=r"$\phi$")
    ax[0][2].set_xlabel("Gamma")
    ax[0][2].set_ylabel("Phi")
    ax[0][2].set_title(f"Event 0 Jet {jet}: Phi History")
    ax[0][2].legend()
    ax[0][2].grid(True)

    ax[1][0].scatter(d, g, label=r"gamma")
    ax[1][0].set_xlabel("Distance")
    ax[1][0].set_ylabel("Gamma")
    ax[1][0].set_title(f"Event 0 Jet {jet}: Gamma Evolution")
    ax[1][0].legend()
    ax[1][0].grid(True)

    ax[1][1].scatter(g, m, label=r"$m$")
    ax[1][1].set_xlabel("Gamma")
    ax[1][1].set_ylabel("Mass")
    ax[1][1].set_title(f"Event 0 Jet {jet}: Mass History")
    ax[1][1].legend()
    ax[1][1].grid(True)

    ax[1][2].scatter(g, pt, label=r"$pT$")
    ax[1][2].set_xlabel("Gamma")
    ax[1][2].set_ylabel("pT")
    ax[1][2].set_title(f"Event 0 Jet {jet}: pT History")
    ax[1][2].legend()
    ax[1][2].grid(True)

    plt.savefig(
        "/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/"
        f"variableHist/variableHist_oracle_gamma_Event0_{jet}.png"
    )

    plt.tight_layout()
    plt.show()

    plt.close()    


# Distance-Jet Size Plotting:

jetdist_list = []
jetsize_list = []

for i in range(100):
    try:
        tmp_jetdist = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/distcomp/0/imdist/" + str(i) + "_jet0_imdist.npy")
        tmp_jetsize = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/distcomp/0/radii/" + str(i) + "_jet0_radii.npy")
        jetdist_list.append(tmp_jetdist)
        jetsize_list.append(tmp_jetsize)
    except:
        continue

jetdist = np.concatenate(jetdist_list)
jetsize = np.concatenate(jetsize_list)


plt.hist2d(jetdist, jetsize, cmin=1, bins=(50, 50))
plt.colorbar()
plt.xlabel("Distance")
plt.ylabel("Characteristic Jet Size")
plt.title("Jet Size Evolution")
plt.grid(True)
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/jetSize_evolution.png")
#plt.show()
plt.close()

# Oracle and Matched Invariant Moment printing

# Paths to the two folders
matched_path = "/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/im_matched/"
oracle_path = "/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/im_oracle/"

# Get all .npy files
matched_files = glob.glob(os.path.join(matched_path, "*.npy"))
oracle_files = glob.glob(os.path.join(oracle_path, "*.npy"))


def get_jet_number(filename):
    match = re.search(r"_jet(\d+)_im\.npy$", os.path.basename(filename))
    return int(match.group(1))


def get_event_number(filename):
    match = re.search(r"^(\d+)_jet", os.path.basename(filename))
    return int(match.group(1))


# Keep only jets 0 through 5
matched_files = [
    f for f in matched_files
    if get_jet_number(f) <= 5
]

oracle_files = [
    f for f in oracle_files
    if get_jet_number(f) <= 5
]

# Sort numerically by event number, then jet number
matched_files.sort(key=lambda f: (get_event_number(f), get_jet_number(f)))
oracle_files.sort(key=lambda f: (get_event_number(f), get_jet_number(f)))


# Load arrays into lists
im_matched_list = []
im_oracle_list = []

for filename in matched_files:
    tmp = np.load(filename)
    im_matched_list.append(tmp[0])

for filename in oracle_files:
    tmp = np.load(filename)
    im_oracle_list.append(tmp[0])

im_matched = np.asarray(im_matched_list)
im_oracle = np.asarray(im_oracle_list)

plt.hist2d(im_oracle, im_matched, cmin=1, bins=(50, 50))
plt.colorbar()
plt.xlabel("Oracle Invariant Moment")
plt.ylabel("Matched Invariant Moment")
plt.title("Invariant Moment 1 - Oracle vs Matched")
plt.grid(True)
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/oracle_vs_matched_im_1.png")
plt.show()

# =======================================================
# Debugging Histograms
# =======================================================

deltar_list = []
ptratio_list = []
matchedpt_list = []
oraclept_list = []
score_list = []
dist_list = []
wmass_list = []
wminusmass_list = []
topmass_list = []
topbarmass_list = []
rwmass_list = []
rwminusmass_list = []
rtopmass_list = []
rtopbarmass_list = []
cwmass_list = []
cwminusmass_list = []
ctopmass_list = []
ctopbarmass_list = []


for i in range(1000):
    try:
        tmp_deltar= np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/deltar/" + str(i) + "_deltar.npy")
        tmp_ptratio = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/ptratio/" + str(i) + "_ptratio.npy")
        tmp_matchedpt = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/matchedpt/" + str(i) + "_matchedpt.npy")
        tmp_oraclept = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/oraclept/" + str(i) + "_oraclept.npy")
        tmp_score = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/score/" + str(i) + "_score.npy")
        tmp_dist = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/dist/" + str(i) + "_dist.npy")
        tmp_invmass = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/invmass/" + str(i) + "_invmass.npy")
        tmp_recmass = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/recmass/" + str(i) + "_recmass.npy")
        tmp_clusmass = np.load("/home/smagedov/local/src/sjet-local/GaussianMultijet/npyarrays/debugging/clusmass/" + str(i) + "_clusmass.npy")
        deltar_list.append(tmp_deltar)
        ptratio_list.append(tmp_ptratio)
        matchedpt_list.append(tmp_matchedpt)
        oraclept_list.append(tmp_oraclept)
        score_list.append(tmp_score)
        dist_list.append(tmp_dist)
        wmass_list.append(tmp_invmass[0])
        wminusmass_list.append(tmp_invmass[1])
        topmass_list.append(tmp_invmass[2])
        topbarmass_list.append(tmp_invmass[3])
        rwmass_list.append(tmp_recmass[0])
        rwminusmass_list.append(tmp_recmass[1])
        rtopmass_list.append(tmp_recmass[2])
        rtopbarmass_list.append(tmp_recmass[3])
        cwmass_list.append(tmp_clusmass[0])
        cwminusmass_list.append(tmp_clusmass[1])
        ctopmass_list.append(tmp_clusmass[2])
        ctopbarmass_list.append(tmp_clusmass[3])
    except:
        continue


deltar = np.concatenate(deltar_list)
ptratio = np.log(np.concatenate(ptratio_list))
matchedpt = np.concatenate(matchedpt_list)
oraclept = np.concatenate(oraclept_list)
score = np.concatenate(score_list)
dist = np.concatenate(dist_list)
wmass = np.asarray(wmass_list)
wminusmass = np.asarray(wminusmass_list)
topmass = np.asarray(topmass_list)
topbarmass = np.asarray(topbarmass_list)
rwmass = np.asarray(rwmass_list)
rwminusmass = np.asarray(rwminusmass_list)
rtopmass = np.asarray(rtopmass_list)
rtopbarmass = np.asarray(rtopbarmass_list)
cwmass = np.asarray(cwmass_list)
cwminusmass = np.asarray(cwminusmass_list)
ctopmass = np.asarray(ctopmass_list)
ctopbarmass = np.asarray(ctopbarmass_list)

# -------------------------
# Invariant Mass Histograms
# -------------------------
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
ax[0][0].hist(wmass, bins=50, alpha=0.5, label="Pythia Mass", color='blue')
ax[0][0].hist(rwmass, bins=50, alpha=0.5, label="Reco Mass", color='orange')
ax[0][0].set_xlabel(r"W Mass")
ax[0][0].set_ylabel("Number of matches")
ax[0][0].set_title(f"W Mass Pythia Distribution")
ax[0][0].legend()
ax[0][0].grid(True)

ax[0][1].hist(wminusmass, bins=50, alpha=0.5, label="Pythia Mass", color='blue')
ax[0][1].hist(rwminusmass, bins=50, alpha=0.5, label="Reco Mass", color='orange')
ax[0][1].set_xlabel(r"W Minus Mass")
ax[0][1].set_ylabel("Number of matches")
ax[0][1].set_title(f"W Minus Mass Distribution")
ax[0][1].legend()
ax[0][1].grid(True)

ax[1][0].hist(topmass, bins=50, alpha=0.5, label="Pythia Mass", color='blue')
ax[1][0].hist(rtopmass, bins=50, alpha=0.5, label="Reco Mass", color='orange')
ax[1][0].set_xlabel(r"Top Mass")
ax[1][0].set_ylabel("Number of matches")
ax[1][0].set_title(f"Top Mass Distribution")
ax[1][0].legend()
ax[1][0].grid(True)

ax[1][1].hist(topbarmass, bins=50, alpha=0.5, label="Pythia Mass", color='blue')
ax[1][1].hist(rtopbarmass, bins=50, alpha=0.5, label="Reco Mass", color='orange')
ax[1][1].set_xlabel(r"Top Bar Mass")
ax[1][1].set_ylabel("Number of matches")
ax[1][1].set_title(f"Top Bar Mass Distribution")
ax[1][1].legend()
ax[1][1].grid(True)

plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/debugging/debugging_hist_invmass.png")
plt.show()
plt.close()


# -------------------------
# Clustered Mass Histograms
# -------------------------
fig, ax = plt.subplots(2, 2, figsize=(12, 8))
ax[0][0].hist(cwmass, bins=50, alpha=0.5, label="Clustered Mass", color='red')
ax[0][0].hist(rwmass, bins=50, alpha=0.5, label="Reco Mass", color='orange')
ax[0][0].set_xlabel(r"W Mass")
ax[0][0].set_ylabel("Number of matches")
ax[0][0].set_title(f"W Mass Pythia Distribution")
ax[0][0].legend()
ax[0][0].grid(True)

ax[0][1].hist(cwminusmass, bins=50, alpha=0.5, label="Clustered Mass", color='red')
ax[0][1].hist(rwminusmass, bins=50, alpha=0.5, label="Reco Mass", color='orange')
ax[0][1].set_xlabel(r"W Minus Mass")
ax[0][1].set_ylabel("Number of matches")
ax[0][1].set_title(f"W Minus Mass Distribution")
ax[0][1].legend()
ax[0][1].grid(True)

ax[1][0].hist(ctopmass, bins=50, alpha=0.5, label="Clustered Mass", color='red')
ax[1][0].hist(rtopmass, bins=50, alpha=0.5, label="Reco Mass", color='orange')
ax[1][0].set_xlabel(r"Top Mass")
ax[1][0].set_ylabel("Number of matches")
ax[1][0].set_title(f"Top Mass Distribution")
ax[1][0].legend()
ax[1][0].grid(True)

ax[1][1].hist(ctopbarmass, bins=50, alpha=0.5, label="Clustered Mass", color='red')
ax[1][1].hist(rtopbarmass, bins=50, alpha=0.5, label="Reco Mass", color='orange')
ax[1][1].set_xlabel(r"Top Bar Mass")
ax[1][1].set_ylabel("Number of matches")
ax[1][1].set_title(f"Top Bar Mass Distribution")
ax[1][1].legend()
ax[1][1].grid(True)

plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/debugging/debugging_hist_clusinvmass.png")
plt.show()
plt.close()



# -------------------------
# Delta R histogram
# -------------------------
plt.figure(figsize=(8, 6))
plt.hist(deltar, bins=50)
plt.xlabel(r"$\Delta R$")
plt.ylabel("Number of matches")
plt.title(r"Distribution of $\Delta R$")
plt.grid()
plt.tight_layout()
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/debugging/debugging_hist_deltar.png")
#plt.show()


# -------------------------
# pT ratio histogram
# -------------------------
plt.figure(figsize=(8, 6))
plt.hist(ptratio, bins=50)
plt.xlabel(r"$|log(p_T^{\mathrm{oracle}} / p_T^{\mathrm{cluster}})|$")
plt.ylabel("Number of matches")
plt.title(r"Distribution of $p_T$ Ratio")
plt.grid()
plt.tight_layout()
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/debugging/debugging_hist_ptratio.png")
#plt.show()

# -------------------------
# Matched Jet pT histogram
# -------------------------
plt.figure(figsize=(8, 6))
plt.hist(matchedpt, bins=50)
plt.xlabel(r"$p_T^{\mathrm{cluster}}$")
plt.ylabel("Number of matches")
plt.title(r"Distribution of Matched Jet $p_T$")
plt.grid()
plt.tight_layout()
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/debugging/debugging_hist_matchedpt.png")
#plt.show()

# -------------------------
# Oracle Jet pT histogram
# -------------------------
plt.figure(figsize=(8, 6))
plt.hist(oraclept, bins=50)
plt.xlabel(r"$p_T^{\mathrm{oracle}}$")
plt.ylabel("Number of matches")
plt.title(r"Distribution of Oracle Jet $p_T$")
plt.grid()
plt.tight_layout()
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/debugging/debugging_hist_oraclept.png")
#plt.show()

# -------------------------
# Score histogram
# -------------------------
plt.figure(figsize=(8, 6))
plt.hist(score, bins=50)
plt.xlabel(r"$\Delta R + p_T^{\mathrm{oracle}} / p_T^{\mathrm{cluster}}$")
plt.ylabel("Number of matches")
plt.title(r"Distribution of the Score")
plt.grid()
plt.tight_layout()
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/debugging/debugging_hist_score.png")
#plt.show()

# -------------------------
# Distance histogram
# -------------------------
plt.figure(figsize=(8, 6))
plt.hist(dist, bins=50)
plt.xlabel(r"$Dist$")
plt.ylabel("Number of matches")
plt.title(r"Distribution of Clustering Distance")
plt.grid()
plt.tight_layout()
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/debugging/debugging_hist_dist.png")
#plt.show()

# -------------------------
# Log of Pt vs Delta R
# -------------------------
plt.figure(figsize=(8, 6))
plt.scatter(ptratio, deltar)
plt.xlim(-1, 1)
plt.ylim(0, 0.5)
plt.ylabel(r"$Delta R$")
plt.xlabel(r"$log($p_T^{\mathrm{oracle}} / p_T^{\mathrm{cluster}})$")
plt.title(r"pT Ratio Distribution")
plt.grid()
plt.tight_layout()
plt.savefig("/home/smagedov/local/src/sjet-local/GaussianMultijet/plots/debugging/debugging_deltar_vs_ratio.png")
#plt.show()

