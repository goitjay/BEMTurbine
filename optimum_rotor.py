# 09 Dec 2019
# Author: Jay Prakash Goit
# Compute the Betz Optimum Rotor
# The script uses Wind Energy Explained 2nd Edition Section 3.7 as references

# Assumptions made in thsi analysis
# 1. Wake rotation included as a Flag
# 2. There is no drag, Cd=0
# 3. There is no losses due to a finite number of blades (i.e. no tip loss)
# 4. Axial induction factor (a) =1/3 in each annular stream tube

import numpy as np
import matplotlib.pyplot as plt # plots
import matplotlib.cm as cm
from os import path # e.g. to join one or more path components
import os
import warnings # to suppress warning at specific place
import sys
import re # to rearrange number


# Developed function
from functions_general import cm2inch
from functions_general import numericSort
from functions_general import str2bool
from functions_general import current_function 

#------------------------------
# The main function
#------------------------------
def optimum_rotor():
    
    #====INPUT starts here================
    #setup_file=input("Enter the setup file: ")
    setup_file='setup/optimum_rotor.inp'
    
    #---Read input setup file
    with open(setup_file,'r') as f:
        lines = f.readlines()
        #---save file related
        save_folder= lines[3].split(' ')[0]
        save_file=   lines[4].split(' ')[0]
        save_file=save_folder+'/'+save_file
        
        #---Computation related
        wake_rotation=str2bool(lines[6].split(' ')[0]) # Flag to include wake rotation in optimum rotor calculation. True or False (logic)
        
        #---Rotor parameters
        Drotor =float(lines[11].split(' ')[0]) # Rotor diameter
        Rad_rotor=Drotor/2.0 # Radius of rotor
        Nblade =int(lines[12].split(' ')[0])   # Number of blades
        tsr    =float(lines[13].split(' ')[0]) # Tip speed ratio
        
        #---Airfoil properties
        Nairfoil  =int(lines[17].split(' ')[0]) # Number of airfoils
        Airfoilfile=np.empty((Nairfoil),dtype=object)
        Airfoil_no =np.empty((Nairfoil),dtype=int)
        AoA_col  =np.empty((Nairfoil),dtype=int) # Column with angle of attack in Cl Cd file
        Cl_col     =np.empty((Nairfoil),dtype=int)
        Cd_col     =np.empty((Nairfoil),dtype=int)
        header     =np.empty((Nairfoil),dtype=int)
        footer     =np.empty((Nairfoil),dtype=int)
        
        #---Airoifl file information
        for n in range(Nairfoil):
            Airfoilfile[n]=lines[19+n].split(',')[0]
            Airfoil_no[n] =int(lines[19+n].split(',')[1])-1
            AoA_col[n]    =int(lines[19+n].split(',')[2])-1
            Cl_col[n]     =int(lines[19+n].split(',')[3])-1
            Cd_col[n]     =int(lines[19+n].split(',')[4])-1
            header[n]     =int(lines[19+n].split(',')[5]) #No. lines as header in airfoil dat file
            footer[n]     =int(lines[19+n].split(',')[6]) 
                
        #---Blade properties
        Nnode =int(lines[30].split(' ')[0])
        Rnode =np.empty((Nnode),dtype=float) # Center of node 
        dRnode =np.empty((Nnode),dtype=float) # Span of node
        Airfoil_id=np.empty((Nnode),dtype=int) # Airfoil number to use at the node
        for n in range(Nnode):
            Rnode[n]=float(lines[32+n].split(',')[0])
            dRnode[n]=float(lines[32+n].split(',')[1])
            Airfoil_id[n]=int(lines[32+n].split(',')[2])-1
    #====INPUT ends here==================

    #---Some checks
    # Check if Rnode values in setup files are correct.
    if np.abs(Rnode[-1]-Rad_rotor)/Rad_rotor*100>5:
        print("The last Rnode value ",str(Rnode[-1]), "m is different from Radius of the rotor", str(Rad_rotor),"m by more than",str(5),"%.\nThe calculation will terminate here.")
        sys.exit()

    # Check if dRnode values in setup files are correct.
    if np.abs(np.sum(dRnode)-Rad_rotor)/Rad_rotor*100>5:
        print("Sum of dRnode values of",str(np.sum(dRnode)), "m is different from Radius of the rotor", str(Rad_rotor),"m by more than",str(5),"%.\nThe calculation will terminate here.")
        sys.exit()
    
    
    Chord =np.empty((Nnode),dtype=float) # Chord length of the section
    BlAngleRelWind=np.empty((Nnode),dtype=float) # Angle of relative wind at node position
    BlSectionPitch=np.empty((Nnode),dtype=float) # Pitch/twist angle at node position
    BlAoA =np.empty((Nnode),dtype=float) # Angle of attack with min(Cd/Cl)
    
    #---Run the iteration for each blade node
    for n in range(Nnode):
        
        # read airfoil file
        data_clcd=np.genfromtxt(Airfoilfile[Airfoil_id[n]],skip_header=header[Airfoil_id[n]],skip_footer=footer[Airfoil_id[n]],delimiter="",dtype=float,encoding="ISO-8859-1") # read data file
        AoA=data_clcd[:,AoA_col[Airfoil_id[n]]]
        Cl =data_clcd[:,Cl_col[Airfoil_id[n]]]
        Cd =data_clcd[:,Cd_col[Airfoil_id[n]]]
        
        imin=np.argmin(np.abs(Cd/Cl))
        imax=np.argmax(Cl/Cd)
        BlAoA[n]=AoA[imax] # AoA with maximum Cd/Cl value
        Clopt = Cl[imax] # Cl with maximum Cd/Cl value
        Cdopt = Cd[imax] # Cd with maximum Cd/Cl value
        # BlAoA[n]=7.0;Clopt=1.0; # just for checking
        print('AoA for max Cl/Cd',BlAoA[n],'   Clopt',Clopt,'   Cdopt',Cdopt)
        
        tsr_node= tsr*Rnode[n]/Rad_rotor

        #---Calculate angle of relative wind and chord lenght
        if wake_rotation:
            #---Wake rotation included
            BlAngleRelWind[n]=2.0/3.0*np.arctan(1.0/tsr_node) # Eq. (3.105) of reference
            Chord[n] = 8*np.pi*Rnode[n]*(1.0-np.cos(BlAngleRelWind[n]))/(Nblade*Clopt) # Eq. (3.106)
        else:
            # Wake rotation ignored
            BlAngleRelWind[n]=np.arctan(2.0/(3.0*tsr_node)) # Eq. (3.78) of reference
            Chord[n] = 8*np.pi*Rnode[n]*np.sin(BlAngleRelWind[n])/(3.0*Nblade*Clopt*tsr_node) # Eq. (3.79) 
        
        BlAngleRelWind[n] =BlAngleRelWind[n]*180.0/np.pi
        BlSectionPitch[n]=BlAngleRelWind[n]-BlAoA[n]
        
        del AoA,Cl,Cd,imin,Clopt,Cdopt
    
    save_dat = np.array([Rnode[:],dRnode[:],Chord[:],BlSectionPitch[:],Airfoil_id[:]+1]).T
    
    with open(save_file,'w') as f:
        f.write('Parameters for Betz Optimum blade'+'\n')
        f.write('Rotor diameter:'+str(Drotor)+' Tip speed ratio: '+str(tsr)+'\n')
        f.write('Rnode,dRnode, Chord,Twist, Afoil no.'+'\n')
        f.write(' [m] , [m]  ,  [m] ,[deg],  [-]     '+'\n')
    
    with open(save_file,'ab') as f:
        np.savetxt(f,save_dat[:,:],fmt='%s,%s,%s,%s,%s')
            
    return
#------------------------------



#--To call the main function everytime this file is run                      
if __name__=='__main__':
    optimum_rotor()
