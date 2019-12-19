# 19 Nov 2018
# Author: Jay Prakash Goit
# General functions for processing Lidar data


import numpy as np
import matplotlib.pyplot as plt
from os import path
import datetime
from matplotlib import dates
import os
from matplotlib.ticker import FuncFormatter
import re # to rearrange number
import warnings # to suppress warning at specific place
import sys
import matplotlib.cm as cm



#---------------------------
# Convert from cm2inch because Python only understand inch
#----------------------------
def cm2inch(value): # define rule for converting from cm to inch
    return value/2.54 #value is what you need to enter...see below
    hfont = {'fontname':'Helvetica'}
    return
#-----------------------------



#-----------------------------------
# Convert axis value into percentage
#-----------------------------------             
def to_percent(y,position):
    return round(100*y,1)
#----------------------------------- 



#-----------------------------
# Sort file according to 1,2,3.. instead of 1,10,11,..
#-----------------------------
def numericSort(value): # value is name of file
    numbers = re.compile(r'(\d+)') # decimal digits [0-9]
    parts = numbers.split(value) # split string into numbers and non-numbers
    parts[1::2] = map(int,parts[1::2]) # 1::2 means map only numbers
    return parts
#-----------------------------



#------------------------------
# Convert from string to boolean i.e. True and False
#------------------------------
def str2bool(string):
    if string=='True':
        return True
    elif string=='False':
        return False
#------------------------------



#------------------------------
# Current function name
#------------------------------
def current_function():
    return sys._getframe(1).f_code.co_name
#------------------------------
