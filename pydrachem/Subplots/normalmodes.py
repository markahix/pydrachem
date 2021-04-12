import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import prody as prd
import glob

def plot_normal_modes(filename,ax=None,title="",num_of_modes=4,colormap="viridis"):
    """
    Plots a line plot of the first N largest modes (default=4) from a given NMD file.
    Parameters
    ----------
    filename : str
        location of file (.nmd format)
    ax : matplotlib.Axes object
        ax object to plot data
    title : str
        plot title
    num_of_modes : int
        number of modes to plot, starting with the largest and moving down.
    colormap : str
        matplotlib.colors color name
    Returns
    -------
    None
    """
    if ax == None:
        ax = plt.gca()
    if not glob.glob(filename):
        print("File not found: " + filename)
        return
    NMA_data,Atom_Group = prd.parseNMD(filename)
    del Atom_Group
    eigens = NMA_data.getEigvals()
    scales = []
    temp = open(filename)
    lines = temp.readlines()
    temp.close()
    for line in lines:
        if 'mode' in line[:5]:
            scales.append(float(line.split()[:3][-1]))
    cmap = cm.get_cmap(colormap)
    color_range=np.linspace(0,1,num_of_modes)
    for i in range(num_of_modes):
        dataset = [np.linalg.norm(NMA_data.getEigvecs()[:,i][n:n+3])*scales[i]*eigens[i] for n in range(0, len(NMA_data.getEigvecs()[:,i]), 3)]
        ax.plot(dataset,label="Mode "+str(i+1),color=cmap(color_range[i]))
    ax.set_title(title)
    ax.legend()
    ax.set_xlim(0,len(dataset))
    plt.xticks(rotation=90)
    return

def plot_PCA_from_NMD(filename,ax=None,title="",plot_range=None,colormap="viridis"):
    """
    Plots the first two normal modes against each other as a scatter plot.
    Parameters
    ----------
    filename : str
        location of file (.nmd format)
    ax : matplotlib.Axes object
        ax object to plot data
    title : str
        plot title
    plot_range : float
        square range of the plot, in case of specific outliers to be discarded.
    colormap : str
        matplotlib.colors color name
    Returns
    -------
    None
    """
    if ax == None:
        ax = plt.gca()
    if not glob.glob(filename):
        print("File not found: " + filename)
        return
    data = np.genfromtxt(filename,delimiter=None,skip_header=9)
    area = np.pi * (2)**2
    x = data[0,2] * data[0,3:]
    y = data[1,2] * data[1,3:]
    z = -(-x**2 - y**2)
    ax.set_title(title)
    ax.set_xlabel("Mode 1")
    ax.set_ylabel("Mode 2")
    if plot_range != None:
        ax.set_xlim( -1*plot_range,plot_range)
        ax.set_ylim( -1*plot_range,plot_range)
    ax.scatter(x,y,marker='o', s=area, zorder=10, alpha=0.4, c=z, edgecolors = 'black', cmap=colormap)
    plt.xticks(rotation=90)
    return
