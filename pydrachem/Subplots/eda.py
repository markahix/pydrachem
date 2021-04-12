import numpy as np
def plot_EDA(courow,vdwrow,ax=None):
    if ax == None:
        ax = plt.gca()
    totalrow = courow + vdwrow
    BarWidth = 0.33
    cou_x = np.arange(1,len(vdwrow)+1,1)
    vdw_x = [x+BarWidth for x in cou_x]
    tot_x = [x+BarWidth for x in vdw_x]
    ax.set_xlim(0,len(vdwrow)+1)
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_xlabel("Residue Number")
    vdw_label = "Van der Waals Energy - Total: "+str(round(sum(vdwrow),2))+" kcal/mol"
    cou_label = "Coulomb Energy - Total: "+str(round(sum(courow),2))+" kcal/mol"
    tot_label = "Combined Nonbonded Energy - Total: "+str(round(sum(totalrow),2))+" kcal/mol"

    ax.bar(vdw_x,vdwrow,align="center",color="green",label=vdw_label)
    ax.bar(cou_x,courow,align="center",color="blue",label=cou_label)
    ax.bar(tot_x,totalrow,align="center",color="red",label=tot_label)

    ax.legend()

    return ax
