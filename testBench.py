from scipy.stats import poisson
import math
import random
import matplotlib.pyplot as plt
import numpy as np
import statistics

Anaught = 400
"""""
#Run the main Simulation here
#Tyler Barker
#3-24-2021
#Demo.py

NEEDED DIRECTORIES under src:
Cells
Concentrations
Output
Output/Data
Output/Plots
Receptors
Receptors/A
Receptors/B
"""""

#from scipy import special
import cell
import enviornment
import utils
import simulation
import configuration
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
ConcParams = [5, "periodic", 0, 100000.0, 35, 70, "static", 10 * ((2*math.pi)/SimParams[0]), 500, -25.5 * scale, 25.5 * scale, 0, 50]

#cell meta stats: [ Absorbtion rate, receptorConsumptionRate, survivalCost, velocities, mutate bool, decisiontype ("naive", "measured", "prediction", "non"-no change), noise(as a protion of the binomial distribution), receptor_mode(simulate, gaussian), scale, combined_portion, ABdivide_equal]
    #Stats: [Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distTrav, StartID]
    #cell location Stats: [array of Locations (empty if everywhere)(zero indexed)

cell_stress = .01

CellMetaStats = [1.0/cell_stress, 0.0, 250, 1.0, "non", "non", 0, "simulate", 1, 1, 1]
CellStats = [200,200,400,0,0,0,5,1,0,0]
Celllocations = []

# Enviornemnt paramter array: [simulationLengthTime, full divide, full Die, sim Time Step, sim location step]
EnviornmentParams = [20, True, True, 0.1, locationStep]

CellMetaStats[5] = "prediction"
simLength = EnviornmentParams[0]
simRepeat = 1
parameters = [SimParams, ConcParams, CellMetaStats, CellStats, Celllocations, EnviornmentParams, simRepeat]

lam = 0.3
gamma = 0.2

T = 10*lam/gamma
random.randint(0, T)


non_10_growth = [0.007550252967483332,
0.008356608104132808,
0.007597482605685589,
0.006802195751832659,
0.006749470255261589,
0.009528463261384242]
measured_10_growth = [0.057359123275275155,
0.05665024431260903,
0.05803877789914299,
0.05782256065515469,
0.0579338186519314,
0.057744743166179564]

adjusted_10_growth = [0.0036223413507684844,
0.003783452954950896,
0.0028494502146929226]

non_20_growth = [0.05939902181160824,
0.05857036110087121,
0.06039139793412079]

measured_20_growth = [0.12536553238074455,
0.12607505383818338,
0.1259488714791388]

adjusted_20_growth = [0.025445670948148508,
0.01630958073684221,
0.027062962986337067]

non_y = []
measured_y = []
adjusted_y = []
non10sum = 0.0
measured10sum = 0.0
adjusted10sum = 0.0

non20sum = 0.0
measured20sum = 0.0
adjusted20sum = 0.0
for i in range(3):
    non10sum += non_10_growth[i]
    measured10sum += measured_10_growth[i]
    adjusted10sum += adjusted_10_growth[i]

    non20sum += non_20_growth[i]
    measured20sum += measured_20_growth[i]
    adjusted20sum += adjusted_20_growth[i]

non_y.append(non10sum/3)
non_y.append(non20sum/3)

measured_y.append(measured10sum/3)
measured_y.append(measured20sum/3)

adjusted_y.append(adjusted10sum/3)
adjusted_y.append(adjusted20sum/3)


test_list = non_10_growth
mean = sum(test_list) / len(test_list)
res_non10 = sum((i - mean) ** 2 for i in test_list) / len(test_list)
test_list = measured_10_growth
mean = sum(test_list) / len(test_list)
res_measured10 = sum((i - mean) ** 2 for i in test_list) / len(test_list)
test_list = adjusted_10_growth
mean = sum(test_list) / len(test_list)
res_adjusted10 = sum((i - mean) ** 2 for i in test_list) / len(test_list)

test_list = non_20_growth
mean = sum(test_list) / len(test_list)
res_non20 = sum((i - mean) ** 2 for i in test_list) / len(test_list)
test_list = measured_20_growth
mean = sum(test_list) / len(test_list)
res_measured20 = sum((i - mean) ** 2 for i in test_list) / len(test_list)
test_list = adjusted_20_growth
mean = sum(test_list) / len(test_list)
res_adjusted20 = sum((i - mean) ** 2 for i in test_list) / len(test_list)

y_error_non = [res_non10 * 20, res_non20]
y_error_measured = [res_measured10, res_measured20]
y_error_adjusted = [res_adjusted10, res_adjusted20]

# plotting graph
x = [1,2]
x_labels = ["1", "1/2"]
plt.plot(x, non_y, label = "Maximum Sensitivity Strategy", marker = 'o')
plt.plot(x, measured_y, label = "Adaptive Strategy", marker = 'o')
plt.plot(x, adjusted_y, label = "Adjusted Maximum Sensitivity Strategy", marker = 'o')
plt.xlabel("Relative Cell Stress")
plt.ylabel("Growth Rate")
plt.title("Growth")
plt.legend()
plt.xticks(x, x_labels)
plt.show()









non_10_growth = [3.865406622256547,
4.051378263338248,
4.286931270566979]
measured_10_growth = [1.8559910550674186,
1.902567017100802,
1.8306594878586457]

adjusted_10_growth = [4.195632561785881,
4.648884897475104,
4.913574807958569]

non_20_growth = [2.7555296571667593,
2.6018669766414235,
2.7188753787993427]

measured_20_growth = [2.8019614749635906,
2.745199361381614,
2.68785152908403]

adjusted_20_growth = [3.599975827476044,
3.9143978523850347,
3.6779582804114272]


non_y = []
measured_y = []
adjusted_y = []
non10sum = 0.0
measured10sum = 0.0
adjusted10sum = 0.0

non20sum = 0.0
measured20sum = 0.0
adjusted20sum = 0.0
for i in range(3):
    non10sum += non_10_growth[i]
    measured10sum += measured_10_growth[i]
    adjusted10sum += adjusted_10_growth[i]

    non20sum += non_20_growth[i]
    measured20sum += measured_20_growth[i]
    adjusted20sum += adjusted_20_growth[i]

non_y.append(non10sum/3)
non_y.append(non20sum/3)

measured_y.append(measured10sum/3)
measured_y.append(measured20sum/3)

adjusted_y.append(adjusted10sum/3)
adjusted_y.append(adjusted20sum/3)


test_list = non_10_growth
mean = sum(test_list) / len(test_list)
res_non10 = sum((i - mean) ** 2 for i in test_list) / len(test_list)
test_list = measured_10_growth
mean = sum(test_list) / len(test_list)
res_measured10 = sum((i - mean) ** 2 for i in test_list) / len(test_list)
test_list = adjusted_10_growth
mean = sum(test_list) / len(test_list)
res_adjusted10 = sum((i - mean) ** 2 for i in test_list) / len(test_list)

test_list = non_20_growth
mean = sum(test_list) / len(test_list)
res_non20 = sum((i - mean) ** 2 for i in test_list) / len(test_list)
test_list = measured_20_growth
mean = sum(test_list) / len(test_list)
res_measured20 = sum((i - mean) ** 2 for i in test_list) / len(test_list)
test_list = adjusted_20_growth
mean = sum(test_list) / len(test_list)
res_adjusted20 = sum((i - mean) ** 2 for i in test_list) / len(test_list)

y_error_non = [res_non10 * 20, res_non20]
y_error_measured = [res_measured10, res_measured20]
y_error_adjusted = [res_adjusted10, res_adjusted20]

# plotting graph
x = [1,2]
x_labels = ["1", "1/2"]
plt.plot(x, non_y, label = "Maximum Sensitivity Strategy", marker = 'o')
plt.plot(x, measured_y, label = "Adaptive Strategy", marker = 'o')
plt.plot(x, adjusted_y, label = "Adjusted Maximum Sensitivity Strategy", marker = 'o')
plt.xlabel("Relative Cell Stress")
plt.ylabel("Mutual Information [$bits$]")
plt.title("Information")
plt.legend()
plt.xticks(x, x_labels)
plt.show()

