#from scipy import special
import cell
import enviornment
import utils
import simulation
import MICalc
import numpy
import random
import randomClass
import matplotlib.pyplot as plt
import math
import pickle
#import numpy as np
import gc
import run
import time
import runSuite
plt.rc('pdf', fonttype=42)
locationStep = 0.5
# size scaled from 100
scale = 1/locationStep
SimParams = [100,-1]
# [diffCoeff, absorb/reflect/periodic, repeatFrequency, magnitude, locationA, locationB]
ConcParams = [5, "periodic", 0, 100000.0, 35, 70, "static", 10 * ((2*math.pi)/(SimParams[0]*scale)), 500, -25.5 * scale, 25.5 * scale, 0, 50]

cell_stress = .01

CellMetaStats = [1.0/cell_stress, 0.0, 250, 1.0, "non", "non", 0, "simulate", 1, 1, 1, 0]
CellStats = [200,200,400,0,0,0,5,1,0,0]
Celllocations = []

# Enviornemnt paramter array: [simulationLengthTime, full divide, full Die, sim Time Step, sim location step]
EnviornmentParams = [20, True, True, 0.1, locationStep]

CellMetaStats[5] = "measured"
simLength = EnviornmentParams[0]
simRepeat = 1
parameters = [SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, simRepeat]

noise_array = [0.0]
cell_stress_array = [0.2, 0.4, 0.6, 0.8, 1]
preset_bool = False

presets = []
presetA = []
presetB = []
presetA_loc = []
presetB_loc = []
presets.append(presetA)
presets.append(presetB)
presets.append(presetA_loc)
presets.append(presetB_loc)

totalCells, MI_entropy, MI_return, presets, MI_sums = run.testRun(SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, False, preset_bool, presets)

print(totalCells)

growth_arr = []
for i in range(len(totalCells)-1):
    growth_arr.append(math.log2(totalCells[i+1]/totalCells[i]))

ave_growth = 0.0
for i in range(len(growth_arr)):
    ave_growth = growth_arr[i]+ave_growth
ave_growth = ave_growth/len(growth_arr)
print(ave_growth)

print(MI_return)
print(MI_sums)

