# 11 Dec 2019
# Author: Jay Prakash Goit
# Compute the wind turbine rotor performance based on Blade Element Momentum Theory
# The script uses Wind Energy Explained 2nd Edition Section 3.6
# and Wind Energy Handbook 2nd Edition Section 3.5 for reference

# Assumptions made in this analysis
# 1. Divides blades into elements of length dr at radial positions from the root r
# 2. The aerodynamic lift and drag force on these elements are responsible for the axial
#    and angular momentum of all of the air which passes through the annulus swept by the blade
# 3. The force on a blade element can by be calculated from aerofoil characteristics using
#    angle of attack and correspodning lift and drag coefficients
# 4. There is no radial interatcion between the flows through the contiguous annuli

#Steps
# 1. Guess the values of axial and angular induction factor
# 2. Calculate the angle of  the relative wind
# 3. Calculate the angle of attack, i.e.
#  AoA= Angle relative wind - twist angle
# 4. Update axial and angular induction factor
# 5. Repeat 1 to 5 until the specified tolerance is reached

import numpy as np
import matplotlib.pyplot as plt # plots
import matplotlib.cm as cm
from os import path # e.g. to join one or more path components
import os
import warnings # to suppress warning at specific place
import sys
import re # to rearrange number
from scipy.optimize import fsolve


# Developed function
from functions_general import cm2inch
from functions_general import numericSort
from functions_general import str2bool
from functions_general import current_function 

#------------------------------
# The main function
#------------------------------
def BEM_analysis():
    
    #====INPUT starts here================
    #setup_file=input("Enter the setup file: ")
    #setup_file='setup/BEM_analysis.inp'
    setup_file='setup/BEM_analysis_KDWT300kW.inp'
    
    #---Read input setup file
    global save_element_param,file_element_param
    
    with open(setup_file,'r') as f:
        lines = f.readlines()

        computation_to_perform=int(lines[2].split(' ')[0]) # 1: tsrCpCt, 2: WindSpeed_Character Pick the number to perform desired computation
        
        #---COMPUTE Cp, Cq AND Ct AS A FUNCTION OF TSR
        save_folder_tsr_CpCt= lines[4].split(' ')[0]
        file_tsr_CpCt=   lines[5].split(' ')[0]
        file_tsr_CpCt=save_folder_tsr_CpCt+'/'+file_tsr_CpCt
        save_element_param=str2bool(lines[6].split(' ')[0]) # Flag for saving element parameters. If true it will save parameters for each element and for tsr
        file_element_param=lines[7].split(' ')[0]
        file_element_param=save_folder_tsr_CpCt+'/'+file_element_param
        tsrmin    =float(lines[8].split(' ')[0]) # Tip speed ratio used for computation
        tsrmax    =float(lines[8].split(' ')[1])
        tsrintrv  =float(lines[8].split(' ')[2])
        tsr=np.arange(tsrmin,tsrmax+tsrintrv,tsrintrv)
        
        Uin       =float(lines[9].split(' ')[0]) # Inflow wind speed
        opt_BladePitch=float(lines[10].split(' ')[0])*np.pi/180.0 # optimum blade pitch angle. Usually 0 degree, but use the design value 
        
        #---COMPUTE POWER AND OTHER PROPERTIES AS A FUNCTION OF WIND SPEED
        file_WindSpeed_input=lines[12].split(' ')[0]
        WindSpeed_col, RotorSpeed_col,BladePitch_col=int(lines[13].split(' ')[0])-1,int(lines[13].split(' ')[1])-1,int(lines[13].split(' ')[2])-1
        header_WindSpeed_inp,footer_WindSpeed_inp=int(lines[14].split(' ')[0]),int(lines[14].split(' ')[1])
        save_folder_WindSpeed_Character= lines[15].split(' ')[0]
        save_file_WindSpeed_Character=   lines[16].split(' ')[0]
        save_file_WindSpeed_Character=save_folder_WindSpeed_Character+'/'+save_file_WindSpeed_Character
        
        #---Computation related
        global a_ax_ini, a_ang_ini, a_tol, a_iter_max, tip_loss, AirDensity
        a_ax_ini  =float(lines[19].split(' ')[0]) # Initial axial induction factor
        a_ang_ini =float(lines[19].split(' ')[1]) # Initial angular induction factor
        a_tol     =float(lines[20].split(' ')[0]) # Tolerance for induction factor calculation
        a_iter_max=int(lines[21].split(' ')[0]) # Maximum iteration for induction factor calculation
        tip_loss=str2bool(lines[22].split(' ')[0]) # Flag for including tip loss correction factor. Only Prandtl's tip loss is implemented. True or False (logic)
        AirDensity=float(lines[23].split(' ')[0]) # Density of Air [kg/m^3]
        
        #---WIND TURBINE PARAMETERS
        global Drotor, Rad_rotor, Nblade, Efficiency_GBox,Efficiency_Generator
        Drotor =float(lines[26].split(' ')[0]) # Rotor diameter
        Rad_rotor=Drotor/2.0 # Radius of rotor
        Nblade =int(lines[27].split(' ')[0])   # Number of blades
        Efficiency_GBox=float(lines[28].split(' ')[0]) # Gearbox efficiency [%]
        Efficiency_Generator=float(lines[29].split(' ')[0]) # Generator efficiency [%]
        
        #---Airfoil properties
        global Nairfoil, Airfoilfile, Airfoil_no, AoA_col, Cl_col, Cd_col, header_Airfoilfile, footer_Airfoilfile
        Nairfoil  =int(lines[31].split(' ')[0]) # Number of airfoils
        Airfoilfile=np.empty((Nairfoil),dtype=object)
        Airfoil_no =np.empty((Nairfoil),dtype=int)
        AoA_col  =np.empty((Nairfoil),dtype=int) # Column with angle of attack in Cl Cd file
        Cl_col     =np.empty((Nairfoil),dtype=int)
        Cd_col     =np.empty((Nairfoil),dtype=int)
        header_Airfoilfile     =np.empty((Nairfoil),dtype=int)
        footer_Airfoilfile     =np.empty((Nairfoil),dtype=int)
        
        #---Airoifl file information
        for n in range(Nairfoil):
            Airfoilfile[n]=lines[33+n].split(',')[0]
            Airfoil_no[n] =int(lines[33+n].split(',')[1])-1
            AoA_col[n]    =int(lines[33+n].split(',')[2])-1
            Cl_col[n]     =int(lines[33+n].split(',')[3])-1
            Cd_col[n]     =int(lines[33+n].split(',')[4])-1
            header_Airfoilfile[n]     =int(lines[33+n].split(',')[5]) #No. lines as header in airfoil dat file
            footer_Airfoilfile[n]     =int(lines[33+n].split(',')[6]) 
        
        #---Blade properties
        global Nnode, Rnode, dRnode, Chord, Twist, Airfoil_id
        Nnode =int(lines[44].split(' ')[0]) # No. nodes for analysis
        Rnode =np.empty((Nnode),dtype=float) # Center of node 
        dRnode=np.empty((Nnode),dtype=float) # Span of node
        Chord =np.empty((Nnode),dtype=float) # Chord length
        Twist =np.empty((Nnode),dtype=float) # Twist angle/blade section pitch angle
        Airfoil_id=np.empty((Nnode),dtype=int) # Airfoil number to use at the node
        
        for n in range(Nnode):
            Rnode[n] =float(lines[46+n].split(',')[0])
            dRnode[n]=float(lines[46+n].split(',')[1])
            Chord[n] =float(lines[46+n].split(',')[2])
            Twist[n] =float(lines[46+n].split(',')[3])*np.pi/180.0
            Airfoil_id[n]=int(lines[46+n].split(',')[4])-1
        
    global myu
    myu = Rnode/Rad_rotor # r/R
    #====INPUT ends here==================
    
    #---Some checks
    # Check if Rnode values in setup files are correct.
    if np.abs(Rnode[-1]-Rad_rotor)/Rad_rotor*100>5:
        print("WARNING: The last Rnode value ",str(Rnode[-1]), "m is different from Radius of the rotor", str(Rad_rotor),"m by more than",str(5),"%.\nThe calculation will terminate here.")
        sys.exit()

    # Check if dRnode values in setup files are correct.
    if np.abs(np.sum(dRnode)-Rad_rotor)/Rad_rotor*100>15:
        print("WARNING: Sum of dRnode values of",str(np.sum(dRnode)), "m is different from Radius of the rotor", str(Rad_rotor),"m by more than",str(15),"%.\nThe calculation will terminate here.")
        sys.exit()
    #---
    
    
    #---Compute and save Ct, Cq, Cp as a function of tsr
    #---Used for producing Cp-tsr(lambda) curve
    if computation_to_perform==1:
        compute_save_tsr_vs_Ct_Cq_Cp(Uin,opt_BladePitch,tsr,file_tsr_CpCt)
    #---

    #---Used for producing power curve as a function of wind speed
    elif computation_to_perform==2:
        compute_save_WindSpeed_vs_Thrust_Torque_Power(file_WindSpeed_input,WindSpeed_col,RotorSpeed_col,BladePitch_col,header_WindSpeed_inp,footer_WindSpeed_inp,save_file_WindSpeed_Character)
    #---

    else:
        print("WARNING: computation_to_perform parameter is not correct. Choose the correct number correcponding to the desired computation.\nThe calculation will terminate here.")
        sys.exit()
    
    return
#------------------------------



#-----------------------------
#Compute and save Ct, Cq, Cp as a function of tsr
#-----------------------------
def compute_save_tsr_vs_Ct_Cq_Cp(Uin,opt_BladePitch,tsr,file_tsr_CpCt):
    Ct = np.zeros(tsr.shape[0]) # Thrust coefficient
    Cq = np.zeros(tsr.shape[0]) # Thrust coefficient
    Cp = np.zeros(tsr.shape[0]) # power coefficient    
    
    #---Iteration for each tsr
    for j in range(tsr.shape[0]):
        print('==============================')
        print('Computing for tsr ',tsr[j])
        Ang_vel=tsr[j]*Uin/Rad_rotor # Angular velocity
        
        Ct[j],Cq[j],Cp[j]=compute_Ct_Cq_Cp(Uin,opt_BladePitch,Ang_vel,tsr[j])
        
    dat_CpCt = np.array([tsr[:],Ct[:],Cq[:],Cp[:]]).T
    
    with open(file_tsr_CpCt,'w') as f:
        f.write('Ct, Cq and Cp as a function of tip speed ratio'+'\n')
        f.write('Rotor diameter: '+str(Drotor)+'\n')
        f.write('tsr,     Ct,           Cq,                 Cp'+'\n')
    
    with open(file_tsr_CpCt,'ab') as f:
        np.savetxt(f,dat_CpCt[:,:],fmt='%s,%s,%s,%s')


    return
#-----------------------------



#-----------------------------
#Compute and save Thrust, Torque and Power as a function of wind speed
# Wind speed characteristics
#-----------------------------
def compute_save_WindSpeed_vs_Thrust_Torque_Power(file_WindSpeed_input,WindSpeed_col,RotorSpeed_col,BladePitch_col,header_WindSpeed_inp,footer_WindSpeed_inp,save_file):

    data_WindSpeed_input=np.genfromtxt(file_WindSpeed_input,skip_header=header_WindSpeed_inp,skip_footer=footer_WindSpeed_inp,delimiter=",",dtype=float,encoding="ISO-8859-1") # read data file
    
    WindSpeed=data_WindSpeed_input[:,WindSpeed_col] #[m/s]
    RotorSpeed=data_WindSpeed_input[:,RotorSpeed_col] #[rpm]
    BladePitch=data_WindSpeed_input[:,BladePitch_col]*np.pi/180 # Changed to radian
    
    Thrust = np.zeros(WindSpeed.shape[0]) # Thrust 
    Torque = np.zeros(WindSpeed.shape[0]) # Torque
    Power = np.zeros(WindSpeed.shape[0]) # Power 
    Power_generator = np.zeros(WindSpeed.shape[0]) # Generator Power after considering gearbox and generator efficiency
    
    #---Iteration for each wind speed
    for j in range(WindSpeed.shape[0]):
        print('==============================')
        print('Computing for Wind Speed ',WindSpeed[j],"m/s")
        Ang_vel=2*np.pi*RotorSpeed[j]/60.0
        tsr=Ang_vel*Rad_rotor/WindSpeed[j]
        
        Ct,Cq,Cp=compute_Ct_Cq_Cp(WindSpeed[j],BladePitch[j],Ang_vel,tsr)
        Thrust[j] = Ct*0.5*AirDensity*WindSpeed[j]**2*np.pi*Rad_rotor**2 # Thrust [N] 
        Torque[j] = Cq*0.5*AirDensity*WindSpeed[j]**2*np.pi*Rad_rotor**3 # Torque [Nm]
        Power[j]  = Cp*0.5*AirDensity*WindSpeed[j]**3*np.pi*Rad_rotor**2 # Power [W]
        
    Power_generator=Efficiency_GBox/100*Efficiency_Generator/100*Power
    save_dat = np.array([WindSpeed[:],RotorSpeed[:],BladePitch[:]*180/np.pi,Thrust[:],Torque[:],Power[:],Power_generator[:]]).T
    
    with open(save_file ,'w') as f:
        f.write('Wind turbine output as a function of wind speed'+'\n')
        f.write('Rotor diameter: '+str(Drotor)+'\n')
        f.write('Wind speed (m/s), Rotor speed (rpm), Blade Pitch Angle (deg), Thrust (N), Torque (Nm), Power Rotor (W), Power generator (W)'+'\n')
    
    with open(save_file,'ab') as f:
        np.savetxt(f,save_dat[:,:],fmt='%s,%s,%s,%s,%s,%s,%s')


    return
#-----------------------------



#-----------------------------
# BEM analysis: Compute thrust, torque and power coefficients using BEM analysis
# Uin: Inflow streamwise velocity
#-----------------------------
def compute_Ct_Cq_Cp(Uin,BladePitch,Ang_vel,tsr):
    
    Ct_const=8/Rad_rotor**2 
    Cq_const=8*Ang_vel/(Uin*Rad_rotor**3)
    Cp_const=8*Ang_vel**2/(Uin**2*Rad_rotor**2)

    Ct,Cq,Cp=0.0,0.0,0.0 # Used for integrtion
    
    #---Run the iteration for each blade node
    for n in range(Nnode):
        
        # read airfoil file
        data_clcd=np.genfromtxt(Airfoilfile[Airfoil_id[n]],skip_header=header_Airfoilfile[Airfoil_id[n]],skip_footer=footer_Airfoilfile[Airfoil_id[n]],delimiter="",dtype=float,encoding="ISO-8859-1") # read data file
        AoA=data_clcd[:,AoA_col[Airfoil_id[n]]]*np.pi/180
        Cl =data_clcd[:,Cl_col[Airfoil_id[n]]]
        Cd =data_clcd[:,Cd_col[ Airfoil_id[n]]]
        
        tsr_node=tsr*myu[n]
        
        a_ax=a_ax_ini; a_ang=a_ang_ini;
        a_ax_bef=0.0;a_ang_bef=0.0;
        
        for ia in range(a_iter_max):
            #print('ia',ia)
            if np.abs(a_ax-a_ax_bef)<=a_tol and np.abs(a_ang-a_ang_bef)<=a_tol and \
               a_ax>0.0 and a_ang>0.0:
                break
            
            if ia==a_iter_max-1:
                print('Axial induction calculation iteration reached the limit of ', str(a_iter_max),' iterations.')
                print('element number=',n+1,' Rnode=',Rnode[n],' a_ax=',a_ax,' a_ang=',a_ang)
                break
            
            a_ax_bef=a_ax; a_ang_bef=a_ang;
            AngleRelWind=np.arctan((1-a_ax)/((1+a_ang)*tsr_node)) # Eq. (3.63)
            
            AoA_calc=AngleRelWind-Twist[n]-BladePitch # calculated angle of attack using Eq. (3.62)
            i1 = (np.abs(AoA[:]-AoA_calc)).argmin()
            if AoA[i1]>AoA_calc:
                i1=i1-1
                
            Cl_calc=(Cl[i1+1]-Cl[i1])/(AoA[i1+1]-AoA[i1])*(AoA_calc-AoA[i1])+Cl[i1]
            Cd_calc=(Cd[i1+1]-Cd[i1])/(AoA[i1+1]-AoA[i1])*(AoA_calc-AoA[i1])+Cd[i1]
            
            if tip_loss:
                F_tl=2.0/np.pi*np.arccos(np.exp(-Nblade/2*(1-myu[n])/(myu[n]*np.sin(AngleRelWind)))) #Eq. (3.98) Prandtl's tip loss
            else:
                F_tl=1.0
                
            c1= Nblade*Chord[n]/Rad_rotor*(Cl_calc*np.cos(AngleRelWind)+Cd_calc*np.sin(AngleRelWind))/\
                ((np.sin(AngleRelWind))**2*8.0*np.pi*myu[n]*F_tl) #Equating 3.69 & 3.58 and substituting for Urel using 3.64
            
            a_ax=fsolve(func_a_ax,a_ax_bef,args=c1)[0]
            #a_ax=c1/(1+c1) # this can also be used alternative to fsolve function above
            a_ang= (1.0-a_ax)*Nblade*Chord[n]/Rad_rotor*(Cl_calc*np.sin(AngleRelWind)-Cd_calc*np.cos(AngleRelWind))/\
                   ((np.sin(AngleRelWind))**2*8.0*np.pi*myu[n]**2*tsr*F_tl) # Equating Eq. (3.71) with (3.59)
            
        Ct=Ct+Ct_const*a_ax*(1-a_ax)*Rnode[n]*dRnode[n]
        Cq=Cq+Cq_const*a_ang*(1-a_ax)*Rnode[n]**3*dRnode[n]        
        Cp=Cp+Cp_const*a_ang*(1-a_ax)*Rnode[n]**3*dRnode[n]
        
        #---Save a_ax,a_ang etc for each element to file
        if save_element_param:
            dat_element_param=np.array([[Rnode[n],AoA_calc*180/np.pi,Cl_calc,Cd_calc,a_ax,a_ang]])
            
            if n == 0:
                file_element_param1=file_element_param+'_tsr'+str(tsr)+'.dat'
                with open(file_element_param1,'w') as f:
                    f.write('Saves computed parameters for each element along the blade'+'\n')
                    f.write('Rotor diameter: '+str(Drotor)+' tsr: '+str(tsr)+'\n')
                    f.write('Rnode [m], AoA [deg],    Cl,                Cd,                   a_ax,               a_ang'+'\n')
            with open(file_element_param1,'ab') as f:
                np.savetxt(f,dat_element_param[:,:],fmt='%s,%s,%s,%s,%s,%s')
                
            del dat_element_param # for safety
            
    return Ct,Cq,Cp
#-----------------------------



#-----------------------------
# Define the function for calculating axial induction factor
#-----------------------------
def func_a_ax(a_ax,c1):
    return a_ax/(1-a_ax)-c1
#-----------------------------



#--To call the main function everytime this file is run                      
if __name__=='__main__':
    BEM_analysis()
