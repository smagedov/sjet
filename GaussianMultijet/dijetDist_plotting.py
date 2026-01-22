import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

distcutoff = np.arange(0.1, 3.1, 0.1).round(1).tolist()
for d in distcutoff:
    fileprefix = "/home/smagedov/local/src/sjet-master/GaussianDijet/dijetDisStudy/" + str(d) + "/dijetStudyExampleOutput_"
    deltar = np.arange(0.1, 3.1, 0.1).round(1).tolist()
    filelist = []
    for r in deltar:
        file = fileprefix + str(r) + ".txt"
        filelist.append(file)

    arrays = []
    for filename in filelist:
        arr = np.loadtxt(filename, comments="#")
        arrays.append(arr)

    data = np.stack(arrays)

    # Obtain vectors of distances from file arrays
    ranvec = []
    jacvec = []
    mmdvec = []
    for i in range(len(data)):
        tempran = []
        tempjac = []
        tempmmd = []
        for j in range(len(data[i])):
            tempran.append(data[i][j][3])
            tempjac.append(data[i][j][4])
            tempmmd.append(data[i][j][5])
        ranvec.append(tempran)
        jacvec.append(tempjac)
        mmdvec.append(tempmmd)

    ranvec = np.array(ranvec)
    jacvec = np.array(jacvec)
    mmdvec = np.array(mmdvec)

    # Compute mean and standard deviation
    ranmean = np.mean(ranvec, axis=1)
    ranstd = np.std(ranvec, axis=1)
    jacmean = np.mean(jacvec, axis=1)
    jacstd = np.std(jacvec, axis=1)
    mmdmean = np.mean(mmdvec, axis=1)
    mmdstd = np.std(mmdvec, axis=1)

    # X positions
    r = np.arange(0.1, 3.1, 0.1).round(1)

    # Plot with error bars
    plt.ylim([-0.1, 1])
    plt.errorbar(r, ranmean, yerr=ranstd, fmt='o-', capsize=5, label='Rand Distance', color = 'red')
    plt.errorbar(r, jacmean, yerr=jacstd, fmt='o-', capsize=5, label='Jaccard Distance', color = 'yellow')
    plt.errorbar(r, mmdmean, yerr=mmdstd, fmt='o-', capsize=5, label='MMD Distance', color = 'blue')
    plt.xlabel('Delta R')
    plt.ylabel('Index')
    plt.title('Dist Cutoff ' + str(d) + ': Indecies vs. Delta R')
    plt.legend()
    plt.grid(True)
    plt.savefig("plots/distComp/dijetDistStudy_Cutoff_" + str(d) + ".png")
    plt.show()
