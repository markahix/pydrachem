import matplotlib.pyplot as plt
import glob
import numpy as np

def parse_corr_file(corrfile):
    return np.genfromtxt(corrfile)

def correlation(corrfile,ax=None,colormap="RdBu",title="Correlated Movement",fontsize="10"):
    """
    matplotlib-based function to put correlation matrix plots into subplots.
    Parameters
    ----------
    corrfile : str, directed to filename or 2D NumPy array
        datafile formatted as cpptraj matrix output
    ax : matplotlib ax, default: None
        directs output to the given (default: currently active) matplotlib ax.
    colormap : str, matplotlib.cm colormap, default: "viridis"
        sets the colormap to use in the matrix plot.
    Returns
    -------
    None
    """
    if ax == None:
        ax = plt.gca()
    ax.set_xlabel("Residue Number",fontsize=fontsize)
    ax.set_ylabel("Residue Number",fontsize=fontsize)
    if type(corrfile) == str:
        if glob.glob(corrfile):
            data = np.genfromtxt(corrfile)
        else:
            print("File not found: " + corrfile)
            return
    else:
        data = corrfile
    dims = data.shape[0]
    X, Y = np.mgrid[0:dims:complex(0, dims), 0:dims:complex(0, dims)]
    ax.set_title(title)
    im = ax.pcolormesh(X,Y,data,cmap=colormap,vmin=-1.,vmax=1.)
    plt.colorbar(im,ticks=[-1,0,1],ax=ax)
    plt.gca().xaxis.tick_bottom()
    plt.xticks(rotation=90)


