import numpy as np

def cou_file_to_EDA_grid(filename):
    data=np.genfromtxt(filename, skip_header=2,delimiter=None,usecols=(1,2,3))
    rescount = int(data.max(axis=0)[1])
    cou_grid = np.zeros(shape=(rescount,rescount),dtype=float)
    for line in data:
        x = int(line[0]-1)
        y = int(line[1]-1)
        z = float(line[2])
        cou_grid[x][y] = z
        cou_grid[y][x] = z
    return cou_grid

def vdw_file_to_EDA_grid(filename):
    data=np.genfromtxt(filename,delimiter=None,skip_header=1,usecols=(1,2,3))
    rescount = int(data.max(axis=0)[1])
    vdw_grid = np.zeros(shape=(rescount,rescount),dtype=float)
    for line in data:
        x = int(line[0]-1)
        y = int(line[1]-1)
        z = float(line[2])
        vdw_grid[x][y] = z
        vdw_grid[y][x] = z
    return vdw_grid

def single_residue_vdw(vdw,residue_number,buffer_space=2):
    rescount = len(vdw)
    resnum = residue_number -1
    front = resnum + buffer_space
    back = resnum - buffer_space
    # if back < 2:
        # back = 2
    vdwrow = np.zeros(rescount,dtype=float)
    if back > 0:
        vdwrow[:back] = vdw[resnum,:back]## + vdw[:back,resnum]
    if front < len(vdw[0]):
        vdwrow[front:] = vdw[resnum,front:]## + vdw[front:,resnum]
    return vdwrow

def single_residue_cou(cou,residue_number,buffer_space=2):
    rescount = len(cou)
    resnum = residue_number - 1
    front = resnum + buffer_space
    back = resnum - buffer_space
    # if back < 2:
        # back = 2
    courow = np.zeros(rescount,dtype=float)
    if back > 0:
        courow[:back] = cou[resnum,:back]## + cou[:back,resnum]
    if front < len(cou[0]):
        courow[front:] = cou[resnum,front:]## + cou[front:,resnum]
    return courow


