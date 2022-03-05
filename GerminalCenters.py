from atexit import register
import deap
import numpy
import random
from deap import tools, creator, base
import bitstring
import scipy.spatial.distance
import copy

numVirus = 10
epitotesSize = 128 # the size of the epitotes and antibodies
cycletime = 8 # the time to mutate in hours
closeAntibodyPercentage = 50 #values between 0-100
diff = 25 #the difference between the antibody and the virus
def newVirusOrAntibody():
    return bitstring.BitArray(random.randint(0,2^128 - 1), length=128)

def closeAntibody(baseVirus, diff):
    anitbody = copy.deepcopy(baseVirus[random.randint(0,len(baseVirus)-1)])
    for i in range(0,len(anitbody)):
        if random.randint(100,0) <= diff:
            anitbody.invert(i)
    return anitbody

def createAntibodies(baseVirus, closePrecent, diff):
    if random.randint(100,0) <= closePrecent:
        return closeAntibody(baseVirus, diff)
    else: return newVirusOrAntibody()

#we create a base virus. for each antibody, we'll vary between new random numbers and tweaking this
baseVirus = []
for i in range(0,numVirus):
    baseVirus.append(newVirusOrAntibody())

def fitnessFunction(epitotes, antibode):
    max(map(lambda x : scipy.spatial.distance.hamming(x,antibode) * epitotesSize,epitotes))

creator.create("FitnessMax",base.Fitness,weights = (128.0,))
creator.create("BCell",bitstring.BitArray,fitness=creator.FitnessMax)

toolbox = base.Toolbox()

toolbox.register("attr_new_anti",newVirusOrAntibody,baseVirus,closeAntibodyPercentage, diff)
toolbox.register("individual", tools.initIterate, creator.BCell, toolbox.attr_new_anti)
toolbox.register("GerminalCenter", tools.initRepeat, list, toolbox.individual)

