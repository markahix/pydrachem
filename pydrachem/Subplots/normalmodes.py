import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import prody as prd
import glob
def parse_normal_modes(filename):
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
    dataset = []
    for i in range(len(scales)):
        dataset.append([np.linalg.norm(NMA_data.getEigvecs()[:,i][n:n+3])*scales[i]*eigens[i] for n in range(0, len(NMA_data.getEigvecs()[:,i]), 3)])
    return scales, dataset, eigens
  
def Plot_Total_Contributions_Normal_Modes(nmdfile,ax=None,title=""):
    if ax == None:
        ax = plt.gca()
    scales,dataset,eigens = parse_normal_modes(nmdfile)
    sums=[]
    for i in range(len(eigens)):
        sums.append(sum(eigens[:i]))
    x = np.arange(1,100,1)
    ax.plot(x[:10],eigens[:10]/sum(eigens),
           marker=".",
           color="blue")
    ax.set_xlabel("Mode")
    ax.set_ylabel(title)
    ax.set_xlim(1,10)
    ax.set_ylim(0,1)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.set_yticks([0,.5,1.0])
    ax.set_yticklabels(["0%","50%","100%"])
    ax.set_xticks(x[:10])

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

def parsePCA(filename):
    data = np.genfromtxt(filename,delimiter=None,skip_header=9)
    area = np.pi * (2)**2
    x = data[0,2] * data[0,3:]
    y = data[1,2] * data[1,3:]
    z = np.linspace(0,1,len(data[1,3:]))
    return x,y,z

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
    z = np.linspace(0,1,len(data[1,3:]))
    ax.set_title(title)
    ax.set_xlabel("Mode 1")
    ax.set_ylabel("Mode 2")
    if plot_range != None:
        ax.set_xlim( -1*plot_range,plot_range)
        ax.set_ylim( -1*plot_range,plot_range)
    ax.scatter(x,y,marker='o', s=area, zorder=10, alpha=0.4, c=z, edgecolors = 'black', cmap=colormap)
    plt.xticks(rotation=90)
    return

class NormalModeData():
    def __init__(self,filename:str):
        self.scales,self.dataset,self.eigens = parse_normal_modes(filename)
        self.filename = filename
        self.x, self.y, self.z = parsePCA(filename)
    def PCA(self,ax=None,title="",plot_range=None,colormap="viridis"):
        if ax == None:
            ax = plt.gca()
        area = np.pi * (2)**2
        ax.set_title(title)
        ax.set_xlabel("Mode 1")
        ax.set_ylabel("Mode 2")
        if plot_range != None:
            ax.set_xlim( -1*plot_range,plot_range)
            ax.set_ylim( -1*plot_range,plot_range)
        ax.scatter(self.x,self.y,marker='o', s=area, zorder=10, alpha=0.4, c=self.z, edgecolors = 'black', cmap=colormap)
        plt.xticks(rotation=90)
        return
    def plot_modes(self,ax=None,title="",n_modes=4,colormap="viridis"):
        if ax == None:
            ax = plt.gca()
        cmap = cm.get_cmap(colormap)
        color_range=np.linspace(0,1,n_modes)
        for i in range(n_modes):
            ax.plot(self.dataset[i],label="Mode "+str(i+1),color=cmap(color_range[i]))
        ax.set_title(title)
        ax.legend()
        ax.set_xlim(0,len(self.dataset[0]))
        plt.xticks(rotation=90)
        return
    def contributions(self,ax=None,title=""):
        if ax == None:
            ax = plt.gca()
        x = np.arange(1,len(self.eigens)+1,1)
        ax.plot(x,self.eigens/sum(self.eigens),
            marker=".",
            color="blue")
        ax.set_xlabel("Mode")
        ax.set_ylabel(title)
        ax.set_xlim(1,len(self.eigens))
        ax.set_ylim(0,1)
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.tick_params(axis='both', which='minor', labelsize=8)
        ax.set_yticks([0,.5,1.0])
        ax.set_yticklabels(["0%","50%","100%"])
        ax.set_xticks(x[:10])
        return