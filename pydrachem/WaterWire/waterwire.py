import numpy as np
import pytraj as pt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
from scipy.ndimage import gaussian_filter

def normalize_array(array):
    a_min = array.min()
    a_max = array.max()
    norm_array = (array - a_min)/(a_max - a_min)
    return norm_array

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def proton_list(coord_dict,ox_string):
    protons=[]
    for key in coord_dict.keys():
        if ox_string+" " in key:
            protons.append(key.split()[1])
    return protons

def oxygen_list(coord_dict,proton_string):
    oxygens=[]
    for key in coord_dict.keys():
        if key[-len(proton_string):] == proton_string:
            oxygens.append(key.split()[0])
    return oxygens

def Get_Collective_Distances(coord_dict,crunched_S,start_ox,traj,threshold=2.85):
    if type(start_ox) == int:
        start_ox = "@"+str(start_ox+1)
    Oxygen_wire = crunched_S[start_ox].copy()
    oxy_list = [start_ox]
    for key in crunched_S.keys():
        if crunched_S[key].max() > threshold:
            Oxygen_wire += crunched_S[key]
            oxy_list.append(key.split()[0])
    new_oxy_list=[start_ox]
    for oxy in oxy_list:
        protons = proton_list(coord_dict,oxy)
        for prot in protons:
            oxygens = oxygen_list(coord_dict,prot)
            if oxy in oxygens:
                oxygens.remove(oxy)
            for ox in oxygens:
                if ox in oxy_list and ox not in new_oxy_list:
                    new_oxy_list.append(ox)
    print("Water Wire Follows: ",new_oxy_list)
    oxy_list = new_oxy_list.copy()
    coord_list=[]
    for oxy in oxy_list[1:]:
        coord_list.append(crunched_S[oxy].copy()/Oxygen_wire.max())
    dist_list = []
    for i in range(0,len(oxy_list)-1):
        dist_string=oxy_list[i]+" "+oxy_list[i+1]
        new_distances = pt.distance(traj,dist_string)*coord_list[i]
        dist_list.append(new_distances)
    coll_dist = np.zeros(len(dist_list[0]))
    for dist in dist_list:
        coll_dist += dist
    return coll_dist

def Get_Coord_Data(infile,topfile,min_OH_dist=1.0,max_OH_dist=1.8) :
    traj = pt.load(infile,top=topfile) #Load parameter and coordinates into trajectory object
    traj.top.set_reference(traj[0]) #Set reference for rest of data (e.g. if you want to align), this aligns to first frame
    closest_atoms = traj.top.select(":1<:6.") #Getting everything within 6A of residue 1, the HPTS molecule.
    OW_indices=[]
    HW_indices=[]
    for atom in traj.top.atoms:    
        if (atom.type == "OW" or atom.type == "oh" or atom.type == "O1") and atom.index in closest_atoms:
            OW_indices.append(atom.index)
        elif (atom.type == "HW" or atom.type == "ho") and atom.index in closest_atoms:
            HW_indices.append(atom.index)
    coord_dict={}
    for oxygen in OW_indices:
        for hydrogen in HW_indices:
            dist_string = "@"+str(oxygen+1)+" @"+str(hydrogen+1)
            distances = pt.distance(traj,dist_string)
            closest_dist = distances.min()
            if closest_dist < max_OH_dist:
                S_val = np.zeros(len(distances),dtype=float)
                for i in range(len(distances)):
                    if distances[i] <= min_OH_dist:
                        S_val[i] = 1
                    elif distances[i] >= max_OH_dist:
                        S_val[i] = 0
                    else:
                        S_val[i] = 1.0 - ((distances[i]-min_OH_dist)/(max_OH_dist - min_OH_dist))
                coord_dict[dist_string] = S_val
    crunched_S = {}
    for key in coord_dict.keys():
        crunched_key = key.split()[0]
        if crunched_key in crunched_S.keys():
            crunched_S[crunched_key] += coord_dict[key].copy()
        else:
            crunched_S[crunched_key] = coord_dict[key].copy()
    return traj,coord_dict,crunched_S

def Plot_Collective_Dist_v_Energy(collective_distance,energy,trash_threshold=0.25):
    trash =int(trash_threshold*len(collective_distance))
    test_x = np.linspace(collective_distance[trash:].min(),collective_distance[trash:].max(),len(collective_distance[trash:]))
    test_y = np.zeros(len(test_x))
    for i in range(len(test_x)):
        index = find_nearest(collective_distance[trash:],test_x[i])
        test_y[i] = energy[trash:][index]#/total_coord[i]
    fig = plt.figure(figsize=(5,4))
    ax = fig.add_subplot(111)
    ax.plot(test_x,test_y,color="grey",alpha=0.5)
    smoothed = gaussian_filter(test_y,sigma=50)
    ax.plot(test_x,smoothed,color="blue")
    ax.set_xlabel(r"Collective Distance ($\AA$)")
    ax.set_ylabel("Emission Energy (eV)")
    plt.show()

def Water_Wire_Processing(infile,topfile,energy_file):
    traj,coord_dict,crunched_S = Get_Coord_Data(infile,topfile)
    energy = np.genfromtxt(energy_file,usecols=2)
    collective_distance = Get_Collective_Distances(coord_dict,crunched_S,8,traj,threshold=2.8)
    Plot_Collective_Dist_v_Energy(collective_distance,energy)
