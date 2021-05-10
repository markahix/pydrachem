import matplotlib.pyplot as plt
import glob
import numpy as np
def parse_RMSF_file(rmsf_file: str):
    if glob.glob(rmsf_file):
        return np.genfromtxt(rmsf_file,skip_header=1,usecols=1,dtype=float)
    else:
        print("File not found: ",rmsf_file)
        return None        

def rmsf(temp_rmsf_array,ax=None,title="",plotrange=0,num_res=0,barcolor="blue",threshold=0,res_offset=0):
    """
    description
    Parameters
    ----------
    rmsf_array : NumPy array or str filename of array.
        either the array itself or a filename containing the array
    ax : matplotlib ax
        matplotlib.Axes object to plot data on.
    title : str
        Title of the plot.
    plotrange : int
        x-axis maximum limit, usually the number of residues in the array.
    num_res : int
        number of residues in the array
    barcolor : str
        matplotlib.colors color name
    Returns
    -------
    None
    """
    
    if type(temp_rmsf_array) == str:
            if glob.glob(temp_rmsf_array):
                data = np.genfromtxt(temp_rmsf_array,skip_header=1,usecols=1,dtype=float)
            else:
                print("File not found: " + temp_rmsf_array)
                return
    elif type(temp_rmsf_array) != str:
        data = temp_rmsf_array.copy()
    
    if ax == None:
        ax = plt.gca()
    if plotrange != 0:
        ax.set_ylim(0,plotrange)
    if num_res != 0:
        data = data[:num_res]
    ax.set_title(title)
    ax.set_xlim(0+res_offset,len(data)+1+res_offset)
    ax.set_xlabel("Residue")
    ax.set_ylabel(r"RMSF ($\AA$)")
    if threshold != 0:
        ax.bar(np.arange(1+res_offset,len(data)+1+res_offset,1),data,align="center",color="grey",alpha=0.5)
        for i in range(len(data)):
            if (data[i] < threshold) and (data[i] > -threshold):
                data[i] = 0
    ax.bar(np.arange(1,len(data)+1,1),data,align="center",color=barcolor)
