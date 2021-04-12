import matplotlib.pyplot as plt
import glob
import numpy as np
from scipy.ndimage import gaussian_filter
def hbonds_v_time(filename,ax=None,title="",steps_per_ns=100,linecolor="red",minmax=(0,20)):
    """
    Plots hydrogen bonds over time with smoothing function overlaid on raw data.
    Raw data is in grey with 50% transparency and thin lineweight.
    Smoothed data is in colored overlay at 100% opacity.
    Parameters
    ----------
    filename : str or NumPy array of rmsd values
        either an array of data or a filename containing the data
    ax : matplotlib.Axes object, default: None
        axis to direct plot.
    title : str, default=""
        Title for the subplot
    steps_per_ns : int, default: 100
        number of steps per time unit in the array.
    linecolor : str, default: "red"
        matplotlib.colors color name

    Returns
    -------
    None
    """
    if type(filename) == str:
        data = np.genfromtxt(filename,skip_header=1,usecols=1)
    else:
        data = filename
    if ax == None:
        ax = plt.gca()
    x = np.arange(0,len(data),dtype=float)/steps_per_ns
    yavg = np.empty(len(data))
    yavg.fill(np.average(data))

    smoothed = gaussian_filter(data,sigma=100)
    ax.set_title(title)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Hydrogen Bonds")
    ax.set_xlim(0,x.max())
    (ymin,ymax) = minmax
    ax.set_ylim(ymin,ymax)
    ax.plot(x,data,color="grey",alpha=0.5)
    ax.plot(x,smoothed,color=linecolor,alpha=1.0)
    ax.plot(x,yavg,color="black",linestyle="--",label=str(round(np.average(data),3)))
    ax.legend()
    return
