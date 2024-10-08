# This program should simulate and graph muons interacting with the detector

import math
import random as ran
from math import atan, degrees
import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


def gammaEvent(maxCount):
    class GrowingList(list):
        def __setitem__(self, index, value):
            if index >= len(self):
                self.extend([None] * (index + 1 - len(self)))
            list.__setitem__(self, index, value)
    count = 0
    scintillators = [0, 0, 0, 0, 0, 0, 0, 0]
    topFlag = 0
    botFlag = 0
    topEvent = 0
    botEvent = 0
    coincidence = 0
    topXpos = GrowingList()
    topYpos = GrowingList()
    botXpos = GrowingList()
    botYpos = GrowingList()
    slope = GrowingList()
    gamma = GrowingList()
    i = 0
    while i <= maxCount:
        slope[i] = 0
        gamma[i] = 0
        topXpos[i] = 0
        topYpos[i] = 0
        botXpos[i] = 0
        botYpos[i] = 0
        i += 1
    while count < maxCount:
        topgammaX = ran.randint(0, 200)
        topgammaY = ran.randint(0, 200)
        botgammaX = ran.randint(0, 200)
        botgammaY = ran.randint(0, 200)
        muDir = ran.randint(0, 100)
        topXpos[count] = topgammaX
        topYpos[count] = topgammaY
        botXpos[count] = botgammaX
        botYpos[count] = botgammaY
        if (10 <= topgammaX <= 81) & (10 <= topgammaY <= 34):
            topFlag = 1
        if topFlag == 1:
            topEvent += 1
            if 10 <= topgammaY <= 16:
                scintillators[7] += 1
            elif 16 < topgammaY <= 22:
                scintillators[6] += 1
            elif 22 < topgammaY <= 28:
                scintillators[5] += 1
            elif 28 < topgammaY <= 34:
                scintillators[4] += 1
        if (10 <= botgammaX <= 81) & (10 <= botgammaY <= 34):
            botFlag = 1
        if botFlag == 1:
            botEvent += 1
            if 10 <= botgammaY <= 16:
                scintillators[3] += 1
            elif 16 < botgammaY <= 22:
                scintillators[2] += 1
            elif 22 < botgammaY <= 28:
                scintillators[1] += 1
            elif 28 < botgammaY <= 34:
                scintillators[0] += 1
        if (botFlag & topFlag) == 1:
            coincidence += 1
        if muDir > 15:
            if (botgammaX - topgammaX) == 0:
                slope[count] = "undefined"
            elif (botgammaX - topgammaX) != 0:
                slope[count] = (botgammaY - topgammaY) / (botgammaX - topgammaX)
            gamma[count] = degrees(math.atan((botgammaY - topgammaY) / 48))
        if 0 <= muDir <= 15:
            if (topgammaX - botgammaX) == 0:
                slope[count] = "undefined"
            elif (topgammaX - botgammaX) != 0:
                slope[count] = (topgammaY - botgammaY) / (topgammaX - botgammaX)
            gamma[count] = degrees(math.atan((topgammaY - botgammaY) / 48))
        count += 1
    print("Top Events Occured:", "\n", topEvent, "\n", "Bottom Events Occured:", "\n", botEvent, "\n", "Coincidences",
          "\n", coincidence, "\n", "Hits on Each Scintillator:", "\n", scintillators, "\n", "Incident Angles:", "\n", gamma, "\n")
    
    today = datetime.datetime.now()
    txt = open("Monte Carlo Sim.txt",'w')
    txt.write("\n")
    txt.write("Current Run Time:")
    txt.write("\n")
    txt.write(str(today))
    txt.write("\n")
    txt.write("Top Events Triggered:")
    txt.write("\n")
    txt.write(str(topEvent))
    txt.write("\n")
    txt.write("Bottom Events Triggered:")
    txt.write("\n")
    txt.write(str(botEvent))
    txt.write("\n")
    txt.write("Number of Coincidences:")
    txt.write("\n")
    txt.write(str(coincidence))
    txt.write("\n")
    txt.write("Hits on Top Plane Scintillators:")
    txt.write("\n")
    i = 4
    while i <= 7:
        txt.write(str(scintillators[i]))
        txt.write("\n")
        i += 1
    txt.write("Hits on Bottom Plane Scintillators:")
    txt.write("\n")
    i = 0
    while i <= maxCount:
        txt.write(str(slope[i]))
        txt.write("\n")
        i += 1
    txt.write("Incident Angles of Cosmic Rays:")
    txt.write("\n")
    i = 0
    while i <= maxCount:
        txt.write(str(gamma[i]))
        txt.write("\n")
        i += 1
    txt.close()

    x = np.arange(8)
    plt.figure(1)
    plt.bar(x, scintillators)
    plt.xticks(x,('1', '2', '3', '4', '5', '6', '7', '8'))
    plt.xlabel('Detector')
    plt.ylabel('Events')
    plt.title('Combined')
    plt.figure(2)
    plt.hist2d(topXpos, topYpos, bins=maxCount, norm=LogNorm())
    plt.colorbar()
    plt.show()


gammaEvent(10000)

