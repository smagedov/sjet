import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def read_events_from_file(filename):
    with open(filename, 'r') as file:
        events = []
        current_event = []

        for line in file:
            line = line.strip()

            if line.startswith('#'):
                if current_event:  # Save the previous event if it exists
                    events.append(current_event)
                    current_event = []
            elif line:  # Skip empty lines
                # Convert the line into a list of numbers
                data = list(map(float, line.split()))
                current_event.append(data)

        # Append the last event if file doesn't end with '#Event'
        if current_event:
            events.append(current_event)

    return events

filename = 'dijetStudyExampleOutput.txt'
event_arrays = read_events_from_file(filename)

for i in range(len(event_arrays)):
    eta_00 = []
    phi_00 = []
    pt_00 = []
    eta_11 = []
    phi_11 = []
    pt_11 = []
    eta_01 = []
    phi_01 = []
    pt_01 = []
    eta_10 = []
    phi_10 = []
    pt_10 = []
    for j in range(len(event_arrays[i])):
        if (event_arrays[i][j][3] == 0):
            if (event_arrays[i][j][4] == 0):
                pt_00.append(event_arrays[i][j][0])
                eta_00.append(event_arrays[i][j][1])
                phi_00.append(event_arrays[i][j][2])
            else:
                pt_01.append(event_arrays[i][j][0])
                eta_01.append(event_arrays[i][j][1])
                phi_01.append(event_arrays[i][j][2])
        else:
            if (event_arrays[i][j][4] == 1):
                pt_11.append(event_arrays[i][j][0])
                eta_11.append(event_arrays[i][j][1])
                phi_11.append(event_arrays[i][j][2])
            else:
                pt_10.append(event_arrays[i][j][0])
                eta_10.append(event_arrays[i][j][1])
                phi_10.append(event_arrays[i][j][2])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(eta_00, phi_00, np.zeros(len(pt_00)), 0.1, 0.1, pt_00, color='c', label="Reference Jet")
    ax.bar3d(eta_11, phi_11, np.zeros(len(pt_11)), 0.1, 0.1, pt_11, color='b', label="Probe Jet")
    if len(pt_01) != 0:
        ax.bar3d(eta_01, phi_01, np.zeros(len(pt_01)), 0.1, 0.1, pt_01, color='y', label="Reference Misclass")
    if len(pt_10) != 0:
        ax.bar3d(eta_10, phi_10, np.zeros(len(pt_10)), 0.1, 0.1, pt_10, color='r', label="Probe Misclass")
    ax.set_title("Event " + str(i) + " Jet Clustering")
    ax.set_xlabel("eta")
    ax.set_ylabel("phi")
    ax.set_zlabel("pt")
    ax.legend()
    plt.savefig("plots/dijetStudyTrueClustering_" + str(i) + ".png")



#df = pd.read_csv('dijetStudyExampleOutput_cleaned.txt', sep=' ')

#eta = df['eta'].tolist()
#phi = df['phi'].tolist()
#pt = df['pt'].tolist()
#true = df['true_cluster_assignment'].tolist()
#clustering = df['clustering_cluster_assignment'].tolist()

#reference_eta = []
#reference_phi = []
#reference_pt = []

#probe_eta = []
#probe_phi = []
#probe_pt = []

#background_eta = []
#background_phi = []
#background_pt = []

#for i in range(len(true)):
#    if (true[i] == 0):
#        reference_eta.append(eta[i])
#        reference_phi.append(phi[i])
#        reference_pt.append(pt[i])
#    elif (true[i] == 1):
#        probe_eta.append(eta[i])
#        probe_phi.append(phi[i])
#        probe_pt.append(pt[i])
#    else:
#        background_eta.append(eta[i])
#        background_phi.append(phi[i])
#        background_pt.append(pt[i])

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.bar3d(reference_eta, reference_phi, np.zeros(len(reference_pt)), 0.1, 0.1, reference_pt, color='y', label="Reference Jet")
#ax.bar3d(probe_eta, probe_phi, np.zeros(len(probe_pt)), 0.1, 0.1, probe_pt, color='r', label="Probe Jet")
#if len(background_eta) != 0:
#    ax.bar3d(background_eta, background_phi, np.zeros(len(background_pt)), 0.1, 0.1, background_pt, color='c', label="Background")
#ax.set_title("True Clustering")
#ax.set_xlabel("eta")
#ax.set_ylabel("phi")
#ax.set_zlabel("pt")
#plt.savefig("dijetStudyTrueClustering.png")
#plt.show()

#reference_eta.clear()
#reference_phi.clear()
#reference_pt.clear()

#probe_eta.clear()
#probe_phi.clear()
#probe_pt.clear()

#background_eta.clear()
#background_phi.clear()
#background_pt.clear()

#for i in range(len(clustering)):
#    if (true[i] == 0):
#        reference_eta.append(eta[i])
#        reference_phi.append(phi[i])
#        reference_pt.append(pt[i])
#    elif (true[i] == 1):
#        probe_eta.append(eta[i])
#        probe_phi.append(phi[i])
#        probe_pt.append(pt[i])
#    else:
#        background_eta.append(eta[i])
#        background_phi.append(phi[i])
#        background_pt.append(pt[i])

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.bar3d(reference_eta, reference_phi, np.zeros(len(reference_pt)), 0.1, 0.1, reference_pt, color='y', label="Reference Jet")
#ax.bar3d(probe_eta, probe_phi, np.zeros(len(probe_pt)), 0.1, 0.1, probe_pt, color='r', label="Probe Jet")
#if len(background_eta) != 0:
#    ax.bar3d(background_eta, background_phi, np.zeros(len(background_pt)), 0.1, 0.1, background_pt, color='c', label="Background")
#ax.set_title("Algorithm Clustering")
#ax.set_xlabel("eta")
#ax.set_ylabel("phi")
#ax.set_zlabel("pt")
#plt.savefig("dijetStudyAlgoClustering.png")
#plt.show()
