from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from simtk.openmm import app
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
import math
from scipy.ndimage import gaussian_filter
import parmed


def Output_Cleanup(filename):
    lines_seen = set() # holds lines already seen
    with open("temp_clean.txt", "w") as output_file:
        for each_line in open(filename, "r"):
            if each_line not in lines_seen: # check if line is not duplicate
                output_file.write(each_line)
                lines_seen.add(each_line)
    os.system("mv temp_clean.txt "+filename)
    return

def Initialize_System(prmtopfile,inpcrdfile,temperature):
    prmtop = app.AmberPrmtopFile(prmtopfile)
    inpcrd = app.AmberInpcrdFile(inpcrdfile)
    system = prmtop.createSystem(nonbondedMethod = app.PME,
        nonbondedCutoff = 1.0*nanometers, # pylint: disable=undefined-variable
        constraints = app.HBonds,
        ewaldErrorTolerance = 0.001,
        rigidWater = False,
        removeCMMotion = True)
    thermostat = LangevinIntegrator(float(temperature)*kelvin, 1/picoseconds, 0.001*picoseconds)# pylint: disable=undefined-variable
    barostat = MonteCarloBarostat(1.0*atmosphere, float(temperature)*kelvin, 25)
    system.addForce(barostat)
    return prmtop,inpcrd,system,thermostat

def Test_State_File(statefile,outname):
    Output_Cleanup(statefile)
    test = np.genfromtxt(statefile, names = True, delimiter=',')
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
    fig.savefig(outname,dpi=300)
    return

def amber_selection_to_atomidx(structure, selection):
    """
    Converts AmberMask selection [amber-syntax]_ to list of atom indices.

    Parameters
    ----------
    structure : parmed.Structure()
        Structure of the system, used for atom selection.
    selection : str
        AmberMask selection that gets converted to a list of atom indices.

    Returns
    -------
    mask_idx : list of int
        List of atom indices.

    References
    ----------
    .. [amber-syntax] J. Swails, ParmEd Documentation (2015). http://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax

    """
    mask = parmed.amber.AmberMask(structure, str(selection))
    mask_idx = [i for i in mask.Selected()]
    return mask_idx

def Production(topology,system,temperature,checkpoint,bead,step_write_interval=10000,slurm=False):
    thermostat = LangevinIntegrator(temperature*kelvin,1/picoseconds,0.001*picoseconds)# pylint: disable=undefined-variable
    if slurm == False:
        platform = Platform.getPlatformByName("CUDA")
        properties = {"CudaPrecision":"mixed"}
        sim = Simulation(prmtop.topology, system, thermostat,platform,properties)
    elif slurm == True:
        sim = Simulation(prmtop.topology, system, thermostat)
    sim.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
    sim.reporters.append(StateDataReporter("temp.out",step_write_interval,step=True,time=True,potentialEnergy=True,
                         kineticEnergy=True,totalEnergy=True,temperature=True,volume=True,density=True))
    sim.reporters.append(CheckpointReporter("prod_"+str(bead)+".chk",step_write_interval))
    sim.reporters.append(DCDReporter("prod_"+str(bead)+".dcd",step_write_interval))
    with open(checkpoint, 'rb') as f:
        sim.context.loadCheckpoint(f.read())
    sim.step(step_write_interval*100)
    os.system("cat temp.out >> output.log")
    if slurm == True:
        os.system("cp /scratch/$USER/$SLURM_JOBID/* $SLURM_SUBMIT_DIR")
        os.system("rm /scratch/$USER/$SLURM_JOBID/*.dcd")
    os.system("rm temp.out")

def Get_At_Idx(prmtop,residue,atom):
    return prmtop.residues[residue].atoms[atom].idx

def Add_Restraints(system,prminfo,mask_A,mask_B,distance,force):
    restraint = HarmonicBondForce()
    Aindex = amber_selection_to_atomidx(prminfo, mask_A)[0]#.index(1) 
    Bindex = amber_selection_to_atomidx(prminfo, mask_B)[0]#.index(1)
    restraint.addBond(int(Aindex),int(Bindex),float(distance)*angstroms,float(force)*kilocalories_per_mole/angstroms**2)
    system.addForce(restraint)

def Specific_Restraint_File(system,prm_info,filename):
    f=open(filename,"r")
    restraints=f.readlines()
    f.close()
    rest_weight=0.0
    for line in restraints:
        if "weight" in line:
            rest_weight = line.split()[1]
        if "bond" in line:
            row = line.split()[1:]
            Add_Restraints(system,prm_info,row[0],row[1],row[2],rest_weight)
            print("Added distance restraint ("+str(row[2])+" Angstroms) of "+str(rest_weight)+"kcal between "+str(row[0])+" and "+str(row[1]))


def Run_Production(prmtopfile,inpcrdfile,restart,temperature,Start_Bead,Number_of_Steps,restraint_file=None,slurm = False):
    prmtop, inpcrd, system, thermostat = Initialize_System(prmtopfile, inpcrdfile, temperature)
    prm_info = parmed.amber.AmberParm(prmtopfile)
    if restraint_file != None:
        Specific_Restraint_File(system,prm_info,restraint_file)
    Production(prmtop.topology,system,temperature,restart,Start_Bead,10000,slurm)
    bead = Start_Bead + 1
    while bead < Number_of_Steps:
        checkpoint = "prod_"+str(bead-1)+".chk"
        Production(prmtop.topology,system,temperature,checkpoint,bead,10000,slurm)
        bead = bead+1
        print("Completed Production Step "+str(bead))
