import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import math

def OpenMM_statefile_parser(filename,outname=""):
    test = np.genfromtxt(filename,names=True,delimiter=",")
    # test = dataset
    colnames = test.dtype.names[2:]
    x = test['Time_ps']/1000
    num_subplots = len(colnames)
    fig = plt.figure(figsize=(15,num_subplots*1.5),dpi=300)
    for i in range(num_subplots):
        colname = colnames[i]
        ax = fig.add_subplot(math.ceil(num_subplots/2),2,i+1)
        ax.plot(x[20:],test[colname][20:],color="grey",alpha=0.5,label=colname)
        avg = np.sum(test[colname][20:])/len(test[colname][20:])
        data=gaussian_filter(test[colname][20:],sigma=100)
        ax.plot(x[20:],data,color="blue")
        yavg = np.empty(len(data))
        yavg.fill(avg)
        ax.plot(x[20:],yavg,color="red",linestyle="--",label=str(round(avg,3)))
        ax.set_xlim(0,x.max())
        ax.set_xlabel("Time (ns)")
        ax.legend(loc=2)
    fig.subplots_adjust(hspace=0.3)
    if outname != "":
        fig.savefig(outname,dpi=300)
    else:
        plt.show()
    return

