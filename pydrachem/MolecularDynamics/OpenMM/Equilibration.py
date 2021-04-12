#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from simtk.openmm import app
import os
import numpy as np
import sys
#import matplotlib.pyplot as plt
import math
#from scipy.ndimage import gaussian_filter
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
        nonbondedCutoff = 1.0*nanometers,
        constraints = app.HBonds,
        ewaldErrorTolerance = 0.001,
        rigidWater = False,
        removeCMMotion = True)
    thermostat = LangevinIntegrator(float(temperature)*kelvin, 1/picoseconds, 0.001*picoseconds)
    barostat = MonteCarloBarostat(1.0*atmosphere, float(temperature)*kelvin, 25)
    system.addForce(barostat)
    return prmtop,inpcrd,system,thermostat

def Check_Convergence(logfile,tolerance=0.01,history_to_check=100):
    Output_Cleanup(logfile)
    output = np.genfromtxt(logfile,names=True,delimiter=",")
    
    potential = output['Potential_Energy_kJmole'][-history_to_check:]
    kinetic = output['Kinetic_Energy_kJmole'][-history_to_check:]
    totalEnergy = output['Total_Energy_kJmole'][-history_to_check:]
    temperature = output['Temperature_K'][-history_to_check:]
    volume = output['Box_Volume_nm3'][-history_to_check:]
    density = output['Density_gmL'][-history_to_check:]
    
    dec_potential = abs(np.std(potential)/np.average(potential))
    dec_kinetic = abs(np.std(kinetic)/np.average(kinetic))
    dec_total = abs(np.std(totalEnergy)/np.average(totalEnergy))
    dec_temperature = abs(np.std(temperature)/np.average(temperature))
    dec_volume = abs(np.std(volume)/np.average(volume))
    dec_density = abs(np.std(density)/np.average(density))
    print("Convergence Criteria:  All St. Dev. within "+str(float(tolerance*100))+"% of Average")

    bool_potential = bool(dec_potential < tolerance)
    if bool_potential:
        print("Potential Energy Converged")
    else:
        print("Potential Energy Not Converged")
    print(dec_potential*100)

    bool_kinetic = bool(dec_kinetic < tolerance)
    if bool_kinetic:
        print("Kinetic Energy Converged")
    else:
        print("Kinetic Energy Not Converged")
    print(dec_kinetic*100)

    bool_total = bool(dec_total < tolerance)
    if bool_total:
        print("Total Energy Converged")
    else:
        print("Total Energy Not Converged")
    print(dec_total*100)

    bool_temperature = bool(dec_temperature < tolerance)
    if dec_temperature < tolerance:
        print("Temperature Converged")
    else:
        print("Temperature Not Converged")
    print(dec_temperature*100)

    bool_volume = bool(dec_volume < tolerance)
    if bool_volume:
        print("Box Volume Converged")
    else:
        print("Box Volume Not Converged")
    print(dec_volume*100)

    bool_density = bool(dec_density < tolerance)
    if bool_density:
        print("Density Converged")
    else:
        print("Density Not Converged")
    print(dec_density*100)
    
    converged = [bool_potential, bool_kinetic, bool_total, bool_temperature, bool_volume, bool_density]
    if all(converged):
        return True
    else:
        return False


'''def Test_State_File(statefile,outname):
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
'''
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

def Get_Original_Atom_Masses(system,index_array):
    original_atom_masses=[]
    for i in range(len(index_array)):
        original_atom_masses.append(system.getParticleMass(i)._value)
    return original_atom_masses

def Freeze_Atoms(system,index_array,positions,force_constant=100):
    rstrn_force=CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")
    rstrn_force.addGlobalParameter("k",force_constant*kilocalories_per_mole/angstroms**2)
    rstrn_force.addPerParticleParameter("x0")
    rstrn_force.addPerParticleParameter("y0")
    rstrn_force.addPerParticleParameter("z0")
    for i in range(len(index_array)):
        if index_array[i] == 1:
            rstrn_force.addParticle(i,positions[i])
    system.addForce(rstrn_force)
    return len(system.getForces())-1

def Minimize_Simulation(prmtop, inpcrd, system, thermostat, step_write_interval=10000, slurm=False):
    if slurm == False:
        platform = Platform.getPlatformByName("CUDA")
        properties = {"CudaPrecision":"mixed"}
        sim = Simulation(prmtop.topology, system, thermostat,platform,properties)
    elif slurm == True:
        sim = Simulation(prmtop.topology, system, thermostat)
    sim.context.setPositions(inpcrd.getPositions(asNumpy=True))
    sim.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    sim.reporters.append(StateDataReporter("temp.out", step_write_interval, step = True, time = True, 
        potentialEnergy = True, kineticEnergy = True, totalEnergy = True, temperature = True,
        volume = True, density = True))
    sim.reporters.append(CheckpointReporter("temp_0.chk",step_write_interval))
    sim.reporters.append(DCDReporter("temp_0.dcd",step_write_interval))
    sim.minimizeEnergy(maxIterations=2000)
    keep_going=False
    iteration_count=1
    os.system("rm output.log")
    while not keep_going:
        sim.step(step_write_interval*100)
        os.system("cat temp.out >> output.log")
        keep_going = Check_Convergence("output.log")
        print("Completed Iteration "+str(iteration_count)+" of Frozen System")
        iteration_count = iteration_count+1
        if slurm == True:
            os.system("cp /scratch/$USER/$SLURM_JOBID/* $SLURM_SUBMIT_DIR")

    
def Relax_System(topology,system,index_array,temperature,checkpoint,bead,force_constant,force_index, step_write_interval=10000,slurm=False):
    system.removeForce(force_index)
    Freeze_Atoms(system,index_array,positions,force_constant)
    
    thermostat = LangevinIntegrator(temperature*kelvin,1/picoseconds,0.001*picoseconds)
    if slurm == False:
        platform = Platform.getPlatformByName("CUDA")
        properties = {"CudaPrecision":"mixed"}
        sim = Simulation(prmtop.topology, system, thermostat,platform,properties)
    elif slurm == True:
        sim = Simulation(prmtop.topology, system, thermostat)
    sim.context.getState(getPositions=True,enforcePeriodicBox=True).getPositions()
    sim.reporters.append(StateDataReporter("temp.out",step_write_interval,step=True,time=True,potentialEnergy=True,
                         kineticEnergy=True,totalEnergy=True,temperature=True,volume=True,density=True))
    sim.reporters.append(CheckpointReporter("temp_"+str(bead)+".chk",step_write_interval))
    sim.reporters.append(DCDReporter("temp_"+str(bead)+".dcd",step_write_interval))
    with open(checkpoint, 'rb') as f:
        sim.context.loadCheckpoint(f.read())
    converged=False
    while not converged:
        sim.step(step_write_interval*100)
        os.system("cat temp.out >> output.log")
        converged = Check_Convergence("output.log")
        if slurm == True:
            os.system("cp /scratch/$USER/$SLURM_JOBID/*.chk $SLURM_SUBMIT_DIR")
            os.system("cp /scratch/$USER/$SLURM_JOBID/*.dcd $SLURM_SUBMIT_DIR")

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
            print("First Atom in Bond: "+row[0])
            print("Second Atom in Bond: "+row[1])
            print("Distance between atoms: "+row[2]+"\n------\n")
            Add_Restraints(system,prm_info,row[0],row[1],row[2],rest_weight)
            print("Added distance restraint ("+str(row[2])+" Angstroms) of "+str(rest_weight)+"kcal between "+str(row[0])+" and "+str(row[1]))


if __name__ == "__main__":
    if len(sys.argv) > 6:
        slurm = False
        prmtopfile = sys.argv[1]
        inpcrdfile = sys.argv[2]
        frozen_mask = sys.argv[3]
        temperature = float(sys.argv[4])
        Restraint_Force = float(sys.argv[5])
        step_write_interval = int(sys.argv[6])
        if len(sys.argv) > 7:
            restraint_file=sys.argv[7]
        if len(sys.argv) > 8:
            slurm = True

        prmtop, inpcrd, system, thermostat = Initialize_System(prmtopfile, inpcrdfile, temperature)
        index_array = parmed.amber.AmberMask(parmed.load_file(prmtopfile),frozen_mask).Selection()
        original_atom_masses = Get_Original_Atom_Masses(system,index_array)
        positions = inpcrd.getPositions(asNumpy=True)
        prm_info = parmed.amber.AmberParm(prmtopfile)
        if len(sys.argv) > 7:
            Specific_Restraint_File(system,prm_info,sys.argv[7])

        force_index = Freeze_Atoms(system,index_array,positions,Restraint_Force)
        Minimize_Simulation(prmtop,inpcrd,system,thermostat,step_write_interval,slurm)
        bead = 1
        while Restraint_Force > 0.1:
            Restraint_Force = Restraint_Force*.5
            checkpoint = "temp_"+str(bead-1)+".chk"
            Relax_System(prmtop.topology,system,index_array,temperature,
                                 checkpoint,bead,Restraint_Force,force_index,step_write_interval,slurm)
            bead = bead+1
            print("Completed Relaxation at " + str(Restraint_Force)+".")

   #     Test_State_File("output.log","Relax_Restraints_Test.png")
    else:
        print("Syntax Error:  Expected 6 arguments\n\nEquilibration.py <.prmtop> <.inpcrd> <frozen mask> <temperature> <starting restraint> <timesteps to save>[<system_specific_restraints_file> <SLURM>]\n")
