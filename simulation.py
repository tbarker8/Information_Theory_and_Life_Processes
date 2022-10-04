import cell
import enviornment
import MICalc
import utils
import random
import math
from scipy import special
import matplotlib.pyplot as plt

class simulation(object):
    #concnetrations parameter array
    # [diffCoeff, absorb/reflect, repeatFrequency, magnitude, locationA (negative if random), locationB]
        #concentrations diffusion coefficent (0 if no diffuse)
        #conentration boundary conditions (absorb, reflect)
        #concentration degradation or depletion rate due to cells consumption

    #cells parameter array
    #cell meta stats: [ Absorbtion rate, receptorConsumptionRate, survivalCost, velocities]
    #Stats: [Arec, Brec, MaxRec, Amol, Bmol, ATP, biomass, generation, distTrav, ID]
    #cell location Stats: [array of Locations (empty if everywhere)(zero indexed), #TODO add later: cell counts at the locations]

        #cell initilization stats ( A/B rec/mol, ATP, biomass threshold)
        #cell locations and counts
        #cell MetaParameters (receptor binding rate, absorbtion rate, consumption rate, costs...)

    #Enviornemnt paramter array: [simulationLengthTime]
        #simulation length

    #simulation parameters: length, data to record (MI data: input concentrations, resulting velocities) etc..

    #...

    def VonMises(self, var, maximum, offset, length, locationStep):
        arr = []
        kappa = 1/var
        mu = (2 * math.pi / length) * offset
        bessel = (special.iv(0, kappa))
        sum = 0.0
        increment = 2 * math.pi / length
        for i in range(length):
            index = i - length/2
            x = index * increment
            numerator = math.exp(kappa * math.cos(x - mu))
            arr.append(numerator / (2 * math.pi * bessel))
            sum += numerator / (2 * math.pi * bessel)
        arr = (arr / sum)*maximum/locationStep
        return arr

    #def uniform(self, ):

    #reads in and returns an SNR value for a given set of receptors and concentrations
    def findSNRfromTable(self, receptors, concentrations):
        SNRs = []
        g = open("concentrationsSNR.txt", "r")
        allconcentrations = g.readlines()[0].split(",")
        for i in range(len(allconcentrations)):
            allconcentrations[i] = float(allconcentrations[i])
        g.close()
        f = open("allSNRtable.txt", "r")
        lines = f.readlines()
        f.close()
        print(receptors)
        #print(concentrations)
        for i in range(len(receptors)):
            if i == 1:
                print(self.SimEnviornment.cellsAlive())
                print("-------------------------------------------------------------------------------")
            if receptors[i] == 0:
                SNRs.append(0.0)
                continue
            for j in range(0, len(lines)):
                arrline = lines[j].split(",")
                if int(arrline[0]) == int(receptors[i]):
                    for k in range(len(allconcentrations)):
                        if math.floor(10000*allconcentrations[k]) == math.floor(10000*concentrations[i]):
                            SNRs.append(float(arrline[k+1]))
                            break
                    break
        print(SNRs)
        sum = 0.0
        for i in SNRs:
            sum += i
        return (sum/len(SNRs))

    #initalize the simulation
    def __init__(self, ConcParams, CellMetaStats, CellStats, CellLocationStats, EnviornmentParams, SimParams, newRand):
        length = int(SimParams[0] * (1/EnviornmentParams[4]))
        #print(length)
        #concentrations
        diffusionCoefficent = ConcParams[0]
        absorb = False

        if ConcParams[1] == "absorb":
            absorb = True
        self.concentrationFreq = ConcParams[2]
        self.concentrationMag = ConcParams[3]
        self.locationA = ConcParams[4]
        self.locationB = ConcParams[5]
        self.locationC = ConcParams[12]
        ConcentrationA = []
        ConcentrationB = []
        ConcentrationC = []
        for i in range(length):
            ConcentrationA.append(0.0)
            ConcentrationB.append(0.0)
            ConcentrationC.append(0.0)
        if self.locationA > -1:
            ConcentrationA[self.locationA] = self.concentrationMag
            ConcentrationB[self.locationB] = self.concentrationMag
            ConcentrationC[self.locationC] = self.concentrationMag
        else:
            ConcentrationA[newRand.getNewInt(0, length)] = self.concentrationMag
            ConcentrationB[newRand.getNewInt(0, length)] = self.concentrationMag
            ConcentrationC[newRand.getNewInt(0, length)] = self.concentrationMag
        if ConcParams[6] == "static":
            ConcentrationA = self.VonMises(ConcParams[7], ConcParams[8], ConcParams[9], length, EnviornmentParams[4])
            ConcentrationB = self.VonMises(ConcParams[7], ConcParams[8], ConcParams[10], length, EnviornmentParams[4])
            #ret = "concA = ["
            #ret2 = "concB = ["
            #for i in range(0, len(ConcentrationA)-1):
            #    ret = ret + str(ConcentrationA[i]) + ", "
            #    ret2 = ret2 + str(ConcentrationB[i]) + ", "
            #ret = ret + str(ConcentrationA[len(ConcentrationA)-1]) + "]"
            #ret2 = ret2 + str(ConcentrationB[len(ConcentrationB)-1]) + "]"
            #print(ret)
            #print(ret2)
            #print(ret)
            #exit(-1)
            #exit(-1)
            ConcentrationC = self.VonMises(ConcParams[7], ConcParams[8], ConcParams[11], length, EnviornmentParams[4])

        #cell locations data
        CellLocations = CellLocationStats[0]
        if len(CellLocations) == 0:
            CellLocations = range(length-1)
        largestID = CellStats[9]
        Allcells = []
        for i in range(len(CellLocations)):
            currCell = cell.cell(CellStats[0], CellStats[1], CellStats[2], CellStats[3], CellStats[4], CellStats[5], CellStats[6], CellStats[7], CellStats[8], largestID, newRand)
            currCell.AbsorbtionRate = CellMetaStats[0]
            currCell.ReceptorConsumptionRate = CellMetaStats[1]
            currCell.survivalCost = CellMetaStats[2]
            currCell.VelocityMultiplier = CellMetaStats[3]
            currCell.noise = CellMetaStats[6]
            currCell.receptor_mode = CellMetaStats[7]
            currCell.VelocityMultiplier = CellMetaStats[8]
            currCell.Combined_portion = CellMetaStats[9]
            currCell.Divide_portion = CellMetaStats[10]
            currCell.adaptive_ratio = CellMetaStats[11]
            if CellMetaStats[4] == "non":
                currCell.mutate = False
            currCell.decisiontype = CellMetaStats[5]
            Allcells.append(currCell)
            largestID += 1

        self.SimEnviornment = enviornment.enviornment(length, ConcentrationA, ConcentrationB, diffusionCoefficent, Allcells, CellLocations, largestID, newRand, ConcentrationC, EnviornmentParams[4])
        self.SimEnviornment.fullDivide = EnviornmentParams[1]
        self.SimEnviornment.fullDie = EnviornmentParams[2]
        self.SimLength = EnviornmentParams[0]
        self.SimEnviornment.timeStep = EnviornmentParams[3]
        #self.SimEnviornment.locationStep = EnviornmentParams[4]
        self.dataRecord = SimParams
        return

    #If there was diffusion
    def runBeginningConc(self, time):
        for i in range(time):
            self.SimEnviornment.runConcentrationAdjusted()

    def moveRunReceptors(self):
        for i in range(self.SimLength):
            self.SimEnviornment.runCells()
            if i % 10 == 0:
                print(i)
                utils.plotAllCells(self.SimEnviornment, i)
        return self.SimEnviornment.Areceptors, self.SimEnviornment.Breceptors, self.SimEnviornment.allLocations

    def moveRunAB(self):
        Amol = []
        Bmol = []
        Acons =  []
        Bcons = []
        for i in range(self.SimLength):
            self.SimEnviornment.runCells()
            A,B, outA, outB = self.SimEnviornment.giveCellAB()
            self.SimEnviornment.decreaseSpace()

            Amol = Amol + A
            Bmol = Bmol + B
            Acons = Acons + outA
            Bcons = Bcons + outB

            #if i % 1 == 0:
                #print(i)
                #utils.plotAllCells(self.SimEnviornment, i)
        print(self.SimEnviornment.divisions)
        return Amol, Bmol, Acons, Bcons

    #records the movement of a cell in a run
    def moveRun(self):
        locs = []
        vel = 0.0
        multiple = 1.0
        totalCells = []
        divides_per_time = []
        realCells = []
        for i in range(self.SimLength):
            self.SimEnviornment.divides_per_time = 0
            self.SimEnviornment.runCells()
            print(self.SimEnviornment.time )
            locations = self.SimEnviornment.giveALlCellLocations()
            vels = self.SimEnviornment.giveAllCellVelocities()
            cellsAlive = self.SimEnviornment.cellsAlive()
            divides_per_time.append(self.SimEnviornment.divides_per_time*multiple)
            totalCells.append(cellsAlive * multiple)
            realCells.append(cellsAlive)

            if (self.SimEnviornment.decreaseSpace()):
                multiple = multiple * (cellsAlive / self.SimEnviornment.lowSpace)
            for j in range(len(vels)):
                vel = vel + vels[j]

            for j in range(len(locations)):
                locs.append(locations[j])
        print(vel)
        return locs, self.SimEnviornment.divides, self.SimEnviornment.splitLoc, self.SimEnviornment.Amol, self.SimEnviornment.Bmol, self.SimEnviornment.Aconcs, self.SimEnviornment.Bconcs, totalCells, divides_per_time, self.SimEnviornment.mol_times, realCells, self.SimEnviornment.Aconcentrations, self.SimEnviornment.Bconcentrations

    #simulates the enviornment in static concnetration and returns the division rate and MI of the environment
    def staicConcRunTime(self, printRun):
        retDivisions = []
        MIs = []
        Hx = []
        Hxgiveny = []
        ret2 = []
        entr = []
        entr_in = []
        roll = 10
        dataX1 = []
        dataX2 = []
        dataX3 = []
        dataX4 = []
        dataX5 = []
        dataX6 = []

        dataY1 = []
        dataY2 = []
        dataY3 = []
        dataY4 = []
        dataY5 = []
        dataY6 = []


        for i in range(roll):
            dataX1.append([])
            dataX2.append([])
            dataX3.append([])
            dataX4.append([])
            dataX5.append([])
            dataX6.append([])

            dataY1.append([])
            dataY2.append([])
            dataY3.append([])
            dataY4.append([])
            dataY5.append([])
            dataY6.append([])

        totalCells = []
        multiple = 1
        for i in range(self.SimLength/self.SimEnviornment.timeStep):
            #print(i)
            #print(self.SimEnviornment.cellsAlive())
            if i == 1:
                print(self.SimEnviornment.cellsAlive())
                #print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
            rollingEntropy = 0.0
            arr = self.SimEnviornment.Bconcentrations
            self.SimEnviornment.runCells()
            #print("made it here")
            if i % 100 == 0:
                #print(str(i) + " time")
                if printRun:
                    utils.plotAllCells(self.SimEnviornment, i)
                    #utils.plotReceptors(self.SimEnviornment, i)
            cellsAlive = self.SimEnviornment.cellsAlive()
            totalCells.append(cellsAlive*multiple)

            if(self.SimEnviornment.decreaseSpace()):
                multiple = multiple *(cellsAlive/self.SimEnviornment.lowSpace)
            #if i % 5 == 0:
            #    print(str(i) + " time")
            #    print(self.SimEnviornment.cellsAlive())
            if i % (self.SimLength-1) == 0 and i != 0:
                concs = self.SimEnviornment.giveAllCurrentConc()
                sum = 0.0
                for i in range(len(totalCells)):
                    sum += totalCells[i]/(len(totalCells))

                retDivisions.append(sum)

                dataX = []
                #print(self.SimEnviornment.cellsAlive())
                dataX1[i % roll] = concs[0]
                dataX2[i % roll] = concs[1]
                dataX3[i % roll] = concs[2]
                dataX4[i % roll] = concs[3]
                dataX5[i % roll] = concs[8]
                dataX6[i % roll] = concs[9]
                datx1 = []
                datx2 = []
                datx3 = []
                datx4 = []
                datx5 = []
                datx6 = []

                for j in range(len(dataX1)):
                    for k in range(len(dataX1[j])):
                        datx1.append(dataX1[j][k])
                    for k in range(len(dataX2[j])):
                        datx2.append(dataX2[j][k])
                    for k in range(len(dataX3[j])):
                        datx3.append(dataX3[j][k])
                    for k in range(len(dataX4[j])):
                        datx4.append(dataX4[j][k])
                    for k in range(len(dataX5[j])):
                        datx5.append(dataX5[j][k])
                    for k in range(len(dataX6[j])):
                        datx6.append(dataX6[j][k])

                dataX.append(datx1)
                dataX.append(datx2)
                dataX.append(datx3)
                dataX.append(datx4)
                #dataX.append(datx5)
                #dataX.append(datx6)

                dataY = []
                dataY1[i % roll] = concs[4]
                dataY2[i % roll] = concs[5]
                dataY3[i % roll] = concs[6]
                dataY4[i % roll] = concs[7]
                dataY5[i % roll] = concs[10]
                dataY6[i % roll] = concs[11]
                daty1 = []
                daty2 = []
                daty3 = []
                daty4 = []
                daty5 = []
                daty6 = []
                for j in range(len(dataY1)):
                    for k in range(len(dataY1[j])):
                        daty1.append(dataY1[j][k])
                    for k in range(len(dataY2[j])):
                        daty2.append(dataY2[j][k])
                    for k in range(len(dataY3[j])):
                        daty3.append(dataY3[j][k])
                    for k in range(len(dataY4[j])):
                        daty4.append(dataY4[j][k])
                    for k in range(len(dataY5[j])):
                        daty5.append(dataY5[j][k])
                    for k in range(len(dataY6[j])):
                        daty6.append(dataY6[j][k])
                dataY.append(daty1)
                dataY.append(daty2)
                dataY.append(daty3)
                dataY.append(daty4)
                #dataY.append(daty5)
                #dataY.append(daty6)

                testMI = MICalc.MICalc()
                #print(dataX)
                #print(dataY)
                ret = testMI.AltMI(dataX,dataY)
                ent_ret = testMI.entropy_allX(dataY)
                entr_in_ret = testMI.entropy_allX(dataX)
                #ret = 0;
                #
                print(ret)
                ret2.append(rollingEntropy)
                MIs.append(ret[0])
                Hx.append(ret[1])
                Hxgiveny.append(ret[2])
                entr.append(ent_ret)
                entr_in.append(entr_in_ret)
            if i % 1 == 0:
                self.SimEnviornment.resetDivisions()
        #print(MIs)
        #exit(-1)
        print(totalCells)
        move = self.SimEnviornment.giveMovement()
        sum1 = 0;
        for i in range(len(move)):
            sum1 = sum1 + math.fabs(move[i])
        sum1 = sum1/len(move)
        print(str(sum1) + "------------------------------------------------asdfasdgasfdgasdfvsadvsdavsdagvsfda------------------------------------------")
        print(MIs)
        print("MIs- -----------------------------------------------------------------------------------------------------------")
        return [totalCells, MIs, entr, entr_in, Hx, Hxgiveny]

    def staicConcRunTimeVar(self, printRun):
        retDivisions = []
        MIs = []
        Hx = []
        Hxgiveny = []
        ret2 = []
        entr = []
        entr_in = []
        roll = 10
        dataX1 = []
        dataX2 = []
        dataX3 = []
        dataX4 = []
        dataX5 = []
        dataX6 = []

        dataY1 = []
        dataY2 = []
        dataY3 = []
        dataY4 = []
        dataY5 = []
        dataY6 = []

        enviornment_Stats = []


        for i in range(roll):
            dataX1.append([])
            dataX2.append([])
            dataX3.append([])
            dataX4.append([])
            dataX5.append([])
            dataX6.append([])

            dataY1.append([])
            dataY2.append([])
            dataY3.append([])
            dataY4.append([])
            dataY5.append([])
            dataY6.append([])

        totalCells = []
        multiple = 1
        MI_vars = []
        for i in range(math.floor(self.SimLength/self.SimEnviornment.timeStep)):
            #print(i)
            #print(self.SimEnviornment.cellsAlive())
            #if i == 1:
                #print(self.SimEnviornment.cellsAlive())
                #print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")
            rollingEntropy = 0.0
            arr = self.SimEnviornment.Bconcentrations
            self.SimEnviornment.runCells()
            #print("made it here")
            if i % 10 == 0:
                print(str(i) + " time")
                if printRun:
                    utils.plotAllCells(self.SimEnviornment, i)
                    utils.plotAllConcentrat(self.SimEnviornment, i)
            cellsAlive = self.SimEnviornment.cellsAlive()
            totalCells.append(cellsAlive*multiple)

            if(self.SimEnviornment.decreaseSpace()):
                multiple = multiple *(cellsAlive/self.SimEnviornment.lowSpace)
            #if i % 5 == 0:
            #    print(str(i) + " time")
            #    print(self.SimEnviornment.cellsAlive())
            if i == math.floor(self.SimLength/self.SimEnviornment.timeStep)-1:
                concs = self.SimEnviornment.giveAllCurrentConc()
                sum = 0.0
                for i in range(len(totalCells)):
                    sum += totalCells[i]/(len(totalCells))

                retDivisions.append(sum)

                dataX = []
                #print(self.SimEnviornment.cellsAlive())
                dataX1[i % roll] = concs[0]
                dataX2[i % roll] = concs[1]
                dataX3[i % roll] = concs[2]
                dataX4[i % roll] = concs[3]
                dataX5[i % roll] = concs[8]
                dataX6[i % roll] = concs[9]
                datx1 = []
                datx2 = []
                datx3 = []
                datx4 = []
                datx5 = []
                datx6 = []

                for j in range(len(dataX1)):
                    for k in range(len(dataX1[j])):
                        datx1.append(dataX1[j][k])
                    for k in range(len(dataX2[j])):
                        datx2.append(dataX2[j][k])
                    for k in range(len(dataX3[j])):
                        datx3.append(dataX3[j][k])
                    for k in range(len(dataX4[j])):
                        datx4.append(dataX4[j][k])
                    for k in range(len(dataX5[j])):
                        datx5.append(dataX5[j][k])
                    for k in range(len(dataX6[j])):
                        datx6.append(dataX6[j][k])

                dataX.append(datx1)
                dataX.append(datx2)
                dataX.append(datx3)
                dataX.append(datx4)
                #dataX.append(datx5)
                #dataX.append(datx6)

                dataY = []
                dataY1[i % roll] = concs[4]
                dataY2[i % roll] = concs[5]
                dataY3[i % roll] = concs[6]
                dataY4[i % roll] = concs[7]
                dataY5[i % roll] = concs[10]
                dataY6[i % roll] = concs[11]
                daty1 = []
                daty2 = []
                daty3 = []
                daty4 = []
                daty5 = []
                daty6 = []
                for j in range(len(dataY1)):
                    for k in range(len(dataY1[j])):
                        daty1.append(dataY1[j][k])
                    for k in range(len(dataY2[j])):
                        daty2.append(dataY2[j][k])
                    for k in range(len(dataY3[j])):
                        daty3.append(dataY3[j][k])
                    for k in range(len(dataY4[j])):
                        daty4.append(dataY4[j][k])
                    for k in range(len(dataY5[j])):
                        daty5.append(dataY5[j][k])
                    for k in range(len(dataY6[j])):
                        daty6.append(dataY6[j][k])
                dataY.append(daty1)
                dataY.append(daty2)
                dataY.append(daty3)
                dataY.append(daty4)
                #dataY.append(daty5)
                #dataY.append(daty6)

                testMI = MICalc.MICalc()
                #print(dataX)
                #print(dataY)
                print(len(self.SimEnviornment.giveAllCells()))
                ret = testMI.AltMI(dataX,dataY)
                ent_ret = testMI.entropy_allX(dataY)
                entr_in_ret = testMI.entropy_allX(dataX)
                #ret = 0;
                #
                ret2.append(rollingEntropy)
                MIs.append(ret[0])
                Hx.append(ret[1])
                Hxgiveny.append(ret[2])
                entr.append(ent_ret)
                entr_in.append(entr_in_ret)

            input_var = self.SimEnviornment.Recep_Var_Arr_input
            output_var = self.SimEnviornment.Recep_Var_Arr_Output

            if i == math.floor(self.SimLength/self.SimEnviornment.timeStep)-1:
                MIs_var_array = []
                MI_var = MICalc.MICalc()
                arr = []
                for k in range(len(input_var)):
                    if len(input_var[k][0]) == 0:
                        continue
                    input_var_MI = [input_var[k][0], input_var[k][1], input_var[k][2], input_var[k][3]]
                    output_var_MI = [output_var[k][0], output_var[k][1], output_var[k][2], output_var[k][3]]
                    var_length = len(input_var[k][0])
                    curr_MI = MI_var.AltMI(input_var_MI, output_var_MI)[0]
                    for j in range(var_length):
                        MIs_var_array.append(curr_MI)
                    arr.append(len(input_var[i][1]))
                mean = 0.0
                for k in range(len(MIs_var_array)):
                    mean = mean + MIs_var_array[k]
                mean = mean/len(MIs_var_array)
                final_var = 0.0
                for k in range(len(MIs_var_array)):
                    final_var = final_var + math.pow(math.fabs(MIs_var_array[k] - mean),2)
                final_var = math.sqrt(final_var/len(MIs_var_array))
                MI_vars.append(final_var)
                plt.grid()
                plt.hist(MIs_var_array, density = True, bins = 50)
                plt.clf()

            if i % 1 == 0:
                self.SimEnviornment.resetDivisions()

            enviornment_Stats.append(self.SimEnviornment.getAllStats())

        move = self.SimEnviornment.giveMovement()
        sum1 = 0;
        for i in range(len(move)):
            sum1 = sum1 + math.fabs(move[i])
        sum1 = sum1/len(move)
        return [totalCells, MIs, entr, entr_in, Hx, Hxgiveny, MI_vars, enviornment_Stats]



    #runs the simulation and calculates the MI, Division rate and calculated SNR
    def staicConcRunTimeSNR(self, printRun):
        retDivisions = []
        MIs = []
        SNRs = []
        for i in range(self.SimLength):
            #print(str(i) + " time")
            self.SimEnviornment.runCells()

            if i % 10 == 0:
                print(str(i) + " time")
                if printRun:
                    utils.plotAllCells(self.SimEnviornment, i)
                    utils.plotReceptors(self.SimEnviornment, i)
            if i % 1 == 0:
                concs = self.SimEnviornment.giveAllCurrentConc()
                allReceptors, concnentrations = self.SimEnviornment.giveAllReceptors()
                if self.SimEnviornment.cellsAlive() > 0:
                    retDivisions.append((self.SimEnviornment.giveDivisions()-self.SimEnviornment.giveDeaths())/self.SimEnviornment.cellsAlive())
                else:
                    retDivisions.append(0)
                dataX = []
                dataX.append(concs[0])
                dataX.append(concs[1])
                dataX.append(concs[2])
                dataX.append(concs[3])
                print(dataX)

                dataY = []
                dataY.append(concs[4])
                dataY.append(concs[5])
                dataY.append(concs[6])
                dataY.append(concs[7])
                testMI = MICalc.MICalc()
                ret = testMI.AltMI(dataX,dataY)
                MIs.append(ret)
                #SNR
                SNR = self.findSNRfromTable(allReceptors, concnentrations)
                SNRs.append(SNR)
        return [retDivisions, MIs, SNRs]


    def traditionalRun(self,time):
        self.runBeginningConc(time)
        frequency = self.concentrationFreq
        for i in range(self.SimLength):
            print("Time: " + str(i) + "-------------------------------------------------------------------------------")
            print(self.SimEnviornment.toString())
            self.SimEnviornment.runConcentrationAdjusted()
            self.SimEnviornment.runCells()
            if i % 10 == 0:
                utils.plotAllCells(self.SimEnviornment, i)
                utils.plotReceptors(self.SimEnviornment, i)
            if frequency > 0:
                if i % frequency == 0:
                    self.SimEnviornment.addAB(self.concentrationMag, self.locationA, self.locationB)

    def test(self, preset_bool, presetA, presetB, presetA_loc, presetB_loc):
        def separate_AB(alphalow, alphahigh, alphas, betas, Al, Alb, Bl, Blb, AR, ARb, BR, BRb):
            ret_al = []
            ret_alb = []
            ret_bl = []
            ret_blb = []
            ret_ar = []
            ret_arb = []
            ret_br = []
            ret_brb = []
            for i in range(len(alphas)):
                curr_alpha = 0.0
                if alphas[i] == 0 and betas[i] == 0:
                    curr_alpha = 0.0
                else:
                    curr_alpha = alphas[i] / (alphas[i] + betas[i])

                if curr_alpha < alphahigh and curr_alpha >= alphalow:
                    ret_al.append(Al[i])
                    ret_alb.append(Alb[i])
                    ret_bl.append(Bl[i])
                    ret_blb.append(Blb[i])
                    ret_ar.append(AR[i])
                    ret_arb.append(ARb[i])
                    ret_br.append(BR[i])
                    ret_brb.append(BRb[i])

            return ret_al, ret_alb, ret_bl, ret_blb, ret_ar, ret_arb, ret_br, ret_brb
        x = range(self.SimEnviornment.length)
        lam = 0.4
        lamprime = lam
        gamma = 0.2
        gammaprime = gamma
        L = math.floor(self.SimEnviornment.length)
        Aknoght = 400
        uniformConc = (lamprime * (Aknoght)) / (gammaprime * L)
        uniAconc = []
        uniBconc = []

        rat_A = []
        rat_B = []

        rat_ALin = []
        rat_ARin = []
        rat_ALout = []
        rat_ARout = []

        rat_BLin = []
        rat_BRin = []
        rat_BLout = []
        rat_BRout = []

        for i in range(len(self.SimEnviornment.Aconcentrations)):
            #uniAconc.append(uniformConc)
            #uniBconc.append(uniformConc)
            pass

        #self.SimEnviornment.Aconcentrations = uniAconc
        #self.SimEnviornment.Bconcentrations = uniBconc
        paths = []
        totalCells = []
        multiple = 1
        MI_entropy = []
        muadj = self.SimEnviornment.timeStep * lam
        #presetA = [51, 100, 240, 244, 289, 326, 374, 379]
        #presetB = [18, 40, 49, 135, 138, 163, 172, 174, 207, 233, 241, 246, 256, 346, 382]
        #presetA_loc = [187, 116, 163, 185, 171, 35, 60, 199]
        #presetB_loc = [141, 39, 14, 20, 3, 13, 157, 2, 191, 100, 55, 171, 197, 119, 155]
        dataX = []
        dataX.append([])
        dataX.append([])
        dataX.append([])
        dataX.append([])

        dataY = []
        dataY.append([])
        dataY.append([])
        dataY.append([])
        dataY.append([])

        if not preset_bool:
            presetA = []
            presetB = []
            presetA_loc = []
            presetB_loc = []

        all_Aratios = []
        for i in range(10):
            all_Aratios.append([])

        for i in range(len(self.SimEnviornment.Aconcentrations)):
            self.SimEnviornment.Aconcentrations[i] = self.SimEnviornment.Aconcentrations[i] * 20
            self.SimEnviornment.Bconcentrations[i] = self.SimEnviornment.Bconcentrations[i] * 20

        for i in range(math.ceil(self.SimLength*(1/self.SimEnviornment.timeStep))):

            print("time" + str(i))

            self.SimEnviornment.runCells()
            concs = self.SimEnviornment.giveAllCurrentConc()

            #A_ratios = self.SimEnviornment.getRatioA()

            #for i in range(len(A_ratios)):
            ##    for j in range(len(A_ratios[i])):
            #        all_Aratios[i].append(A_ratios[i][j])

            for j in range(len(concs[0])):
                dataX[0].append(concs[0][j])
                dataX[1].append(concs[1][j])
                dataX[2].append(concs[2][j])
                dataX[3].append(concs[3][j])

                dataY[0].append(concs[4][j])
                dataY[1].append(concs[5][j])
                dataY[2].append(concs[6][j])
                dataY[3].append(concs[7][j])

            #self.SimEnviornment.runConcentrationAdjusted()
            cellsAlive = self.SimEnviornment.cellsAlive()
            print(cellsAlive)

            #if (i * self.SimEnviornment.timeStep) % 0.1 == 0:
            if True:
                totalCells.append(cellsAlive * multiple)
            if (self.SimEnviornment.decreaseSpace()):
                multiple = multiple * (cellsAlive / self.SimEnviornment.lowSpace)
            rand1 = random.uniform(0, 1)
            rand2 = random.uniform(0, 1)
            #presetA = [3, 5, 83, 85, 121, 124, 125, 136, 203, 224, 246, 254, 280, 312, 324, 326, 340]
            #presetB = [4, 47, 98, 114, 125, 166, 175, 177, 191, 203, 210, 374, 380, 397]
            #preset_bool = False
            if not preset_bool:
                if rand1 < muadj:
                    presetA.append(i)
                    #print("rand1 happened")
                    locationA = random.randint(0, math.floor(self.SimEnviornment.length*(1/self.SimEnviornment.locationStep))-1)
                    presetA_loc.append(locationA)
                    #self.SimEnviornment.addA(Aknoght/self.SimEnviornment.locationStep, locationA)
                if rand2 < muadj:
                    presetB.append(i)

                    #print("rand2 happened")
                    locationB = random.randint(0, math.floor(self.SimEnviornment.length*(1/self.SimEnviornment.locationStep))-1)
                    presetB_loc.append(locationB)
                    #self.SimEnviornment.addB(Aknoght/self.SimEnviornment.locationStep, locationB)
            else:
                for j in range(len(presetA)):
                    if i == presetA[j]:
                        self.SimEnviornment.addA(Aknoght / self.SimEnviornment.locationStep, presetA_loc[j])
                for j in range(len(presetB)):
                    if i == presetB[j]:
                        self.SimEnviornment.addB(Aknoght / self.SimEnviornment.locationStep, presetB_loc[j])

            Aconc = self.SimEnviornment.Aconcentrations
            Bconc = self.SimEnviornment.Bconcentrations
            degradeCoe = gamma*self.SimEnviornment.timeStep
            #for j in range(len(Aconc)):
            #    Aconc[j] = Aconc[j]*(1-degradeCoe)
            #    Bconc[j] = Bconc[j]*(1-degradeCoe)
            #total_cells = len(self.SimEnviornment.giveAllCells()) + 0.00001
            cells = self.SimEnviornment.positionHash
            arry = []
            max = 0
            for pos in cells:
                arry.append(len(cells[pos]))
                max += len(cells[pos])
            Aconc = self.SimEnviornment.Aconcentrations
            Bconc = self.SimEnviornment.Bconcentrations
            Aconci = []
            Bconci = []

            for j in range(len(arry)):
                arry[j] = (arry[j]/max)*300

            for j in range(len(Aconc)):
                Aconci.append(Aconc[j])
                Bconci.append(Bconc[j])
            #plt.grid()
            labels = []
            for j in range(math.floor(self.SimEnviornment.length/10)):
                labels.append(j*10)
            ticks = []
            for j in range(math.floor(self.SimEnviornment.length/10)):
                ticks.append(10*j)
            print(len(Aconci))
            plt.bar(x, Aconci, width=1)
            plt.bar(x, Bconci, width=1)
            time_str = "Time: "+ str(round(i*self.SimEnviornment.timeStep, 2))
            print(len(arry))
            plt.plot(x, arry, color = 'r', label = time_str)
            plt.xticks(ticks,labels)
            plt.legend()
            #plt.ylim(0, 0.5)
            path = "sim_temp/"+ str(i) + ".png"
            plt.savefig(path)
            paths.append(path)
            plt.clf()
            self.SimEnviornment.decreaseSpace()

            #---------------------------------------------------------------------------------------
            if True:
                input_var = self.SimEnviornment.Recep_Var_Arr_input
                output_var = self.SimEnviornment.Recep_Var_Arr_Output
                #if i == math.floor(self.SimLength/self.SimEnviornment.timeStep)-1:
                MIs_var_array = []
                MI_var = MICalc.MICalc()
                arr = []
                for k in range(len(input_var)):
                    if len(input_var[k][0]) == 0:
                        continue
                    input_var_MI = [input_var[k][0], input_var[k][1], input_var[k][2], input_var[k][3]]
                    output_var_MI = [output_var[k][0], output_var[k][1], output_var[k][2], output_var[k][3]]
                    var_length = len(input_var[k][0])
                    curr_MI = MI_var.AltMI(input_var_MI, output_var_MI)[0]
                    for j in range(var_length):
                        MIs_var_array.append(curr_MI)
                    arr.append(len(input_var[i][1]))
                #mean = 0.0
                final_entropy = MI_var.returnSingleEntropy(MIs_var_array, MIs_var_array)
                #for k in range(len(MIs_var_array)):
                #    mean = mean + MIs_var_array[k]
                #mean = mean/len(MIs_var_array)
                ##final_var = 0.0
                #for k in range(len(MIs_var_array)):
                #    final_var = final_var + math.pow(math.fabs(MIs_var_array[k] - mean),2)
                #final_var = math.sqrt(final_var/len(MIs_var_array))
                MI_entropy.append(final_entropy)
                #plt.grid()
                #plt.hist(MIs_var_array, density = True, bins = 50)
                #plt.clf()
            #---------------------------------------------------------------------------------------
        testMI = MICalc.MICalc()
        alpha_beta_divide = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75,
                             0.8, 0.85, 0.9, 0.95, 1]
        MI_sums = []
        #for i in range(len(alpha_beta_divide) - 1):
        #    print(i)
        #    ave_alpha = (alpha_beta_divide[i + 1] + alpha_beta_divide[i]) / 2
        #    curr_al, curr_alb, curr_bl, curr_blb, curr_ar, curr_arb, curr_br, curr_brb = separate_AB( alpha_beta_divide[i], alpha_beta_divide[i + 1], all_Aratios[0], all_Aratios[1], all_Aratios[2], all_Aratios[6], all_Aratios[4], all_Aratios[8], all_Aratios[3], all_Aratios[7], all_Aratios[5], all_Aratios[9])

        #    weighted_MI = (1-ave_alpha) * 0.5 * testMI.AltMI([curr_al], [curr_alb])[0]
        #    weighted_MI = weighted_MI + ave_alpha * 0.5 * testMI.AltMI([curr_bl], [curr_blb])[0]
        #    weighted_MI = weighted_MI + (1-ave_alpha) * 0.5 * testMI.AltMI([curr_ar], [curr_arb])[0]
        #    weighted_MI = weighted_MI + ave_alpha * 0.5 * testMI.AltMI([curr_br], [curr_brb])[0]
        #    MI_sums.append(weighted_MI)

        print("preset A")
        print(presetA)
        print("preset B")
        print(presetB)
        print("presetA_loc")
        print(presetA_loc)
        print("presetB_loc")
        print(presetB_loc)
        import imageio
        images = []

        print("Length of datax")
        print(len(dataX))
        print("Length of datay")
        print(len(dataY))
        MI_return = (testMI.AltMI(dataX, dataY))
        for filename in paths:
            images.append(imageio.imread(filename))
        imageio.mimsave('sim_gifs/movie.gif', images)
        presets = []
        presets.append(presetA)
        presets.append(presetB)
        presets.append(presetA_loc)
        presets.append(presetB_loc)

        return totalCells, MI_entropy, MI_return, presets, MI_sums


