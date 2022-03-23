import csv
import deap
import numpy
import array
import random
from deap import tools, creator, base
import bitstring
from pyrsistent import mutant
import scipy.spatial.distance
import copy

from sympy import Max

global germinalIndex
germinalIndex = 0
numVirus = 30
epitotesSize = 48  # the size of the epitotes and antibodies
cycletime = 8  # the time to mutate in hours
closeAntibodyPercentage = 50  # values between 0-100
diff = 30  # the difference between the antibody and the virus
# whether the multiple germinal centers are fighting the same virus.
sameVirus = False
numberOfGernimalCenters = 12  # The number of germinal centers
numberToTransport = 5  # the number of B Cells to share


def newVirusOrAntibody():
    base = []
    for _ in range(0, epitotesSize):
        base.append(random.randint(0, 1))
    print(base)
    return base


def closeAntibody(baseVirus, diff):
    anitbody = copy.deepcopy(baseVirus[random.randint(0, len(baseVirus)-1)])
    for i in range(0, len(anitbody)):
        if random.randint(0, 100) <= diff:
            anitbody[i] = int(not anitbody[i])
    print(anitbody)
    return anitbody


def createAntibodies(baseVirus, closePrecent, diff):
    if random.randint(0, 100) <= closePrecent:
        return closeAntibody(baseVirus[germinalIndex], diff)
    else:
        return newVirusOrAntibody()


# we create a base virus. for each antibody, we'll vary between new random numbers and tweaking this
baseViruses = []
if(sameVirus):
    baseVirus = []
    for i in range(0, numVirus):
        baseVirus.append(newVirusOrAntibody())
    for j in range(0, numberOfGernimalCenters):
        baseViruses.append(copy.deepcopy(baseVirus))
else:
    for j in range(0, numberOfGernimalCenters):
        baseVirus = []
        for i in range(0, numVirus):
            baseVirus.append(newVirusOrAntibody())
        baseViruses.append(copy.deepcopy(baseVirus))
print(baseViruses[0])

totalDataV1 = []
totalDataVAll = []
totalAveragesV1 = []
totalAveragesVAll = []
for i in range(0, numberOfGernimalCenters+1):
    totalDataV1.append([])
    totalDataVAll.append([])
    totalAveragesV1.append([])
    totalAveragesVAll.append([])
# might want to add a third parameter to account for a specific virus we are targetting.


def fitnessFunction(epitotes, antibode):
    return (singleFitness(epitotes[germinalIndex], antibode), generalFitness(epitotes[germinalIndex], antibode))
# How well against of the virus's in our study?


def singleFitness(epitotes, antibode):
    return max(map(lambda x: epitotesSize - (scipy.spatial.distance.hamming(x, antibode) * epitotesSize), epitotes))
# How well against the whole collection?


def generalFitness(epitotes, antibode):
    return sum(map(lambda x: epitotesSize - (scipy.spatial.distance.hamming(x, antibode) * epitotesSize), epitotes))/(len(epitotes))


def randomsingleMutation(antibody):
    antibody.invert(random.randint(0, len(antibody) - 1))


creator.create("FitnessMax", base.Fitness, weights=(1.0, 10.0))
creator.create("BCell", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()

toolbox.register("attr_new_anti", createAntibodies,
                 baseViruses, closeAntibodyPercentage, diff)
toolbox.register("individual", tools.initIterate,
                 creator.BCell, toolbox.attr_new_anti)
toolbox.register("GerminalCenter", tools.initRepeat, list, toolbox.individual)

toolbox.register("generalEval", lambda x: generalFitness(baseViruses, x))
toolbox.register("evaluate", lambda x: fitnessFunction(baseViruses, x))
toolbox.register("mutate", tools.mutFlipBit, indpb=0.1)
toolbox.register("select", tools.selTournament, tournsize=3)


def main():
    global germinalIndex
    pops = []
    for i in range(0, numberOfGernimalCenters):
        pops.append(toolbox.GerminalCenter(n=100))
    print("Start of evolution")
    MUTPB = 0.5
    for pop in pops:
        fitnesses = list(map(toolbox.evaluate, pop))
        for ind, fit in zip(pop, fitnesses):
            print(fit)
            ind.fitness.values = fit

        print("  Evaluated %i individuals" % len(pop))

        fits = [ind.fitness.values[0] for ind in pop]
        fits2 = [ind.fitness.values[1] for ind in pop]

    # Variable keeping track of the number of generations
    g = 0

    # Begin the evolution

    totaltime = 0
    germinalIndex = 0
    while g < 100:
        # A new generation
        g = g + 1
        totalDataV1[0].append(totaltime)
        totalDataVAll[0].append(totaltime)
        totalAveragesV1[0].append(totaltime)
        totalAveragesVAll[0].append(totaltime)
        for pop in pops:

            print("-- Generation %i -- At time %i For Center %i" %
                  (g, totaltime, germinalIndex))
            # print(pop)
            # Select the next generation individuals
            # I'll need to figure out how it choose individuals and how we can change that.
            offspring = toolbox.select(pop, len(pop))
            # Clone the selected individuals
            offspring = list(map(toolbox.clone, offspring))

            # Apply mutation
            for mutant in offspring:
                # mutate an individual with probability MUTPB
                if random.random() < MUTPB:
                    toolbox.mutate(mutant)
                    del mutant.fitness.values
            if numberOfGernimalCenters > 1:
                otherPop = pops[random.randint(0, len(pops)) % len(pops)]
                random1 = random.randint(0, len(otherPop)) % len(otherPop)
                random2 = (random1 + numberToTransport) % len(otherPop)
                print(random1)
                print(random2)
                if(random1 < random2):
                    exchangedBCells = otherPop[random1: random2]
                    otherPop[random1:random2] = pop[random1:random2]
                    offspring[random1:random2] = exchangedBCells
                    for i in range(random1, random2):
                        del otherPop[i].fitness.values
                        del offspring[i].fitness.values

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            print("  Evaluated %i individuals" % len(invalid_ind))

            # The population is entirely replaced by the offspring
            pop[:] = offspring

            # Gather all the fitnesses in one list and print the stats
            fits = [ind.fitness.values[0] for ind in pop]
            fits2 = [ind.fitness.values[1] for ind in pop]

            length = len(pop)
            mean = sum(fits) / length
            sum2 = sum(x*x for x in fits)
            std = abs(sum2 / length - mean**2)**0.5
            mean2 = sum(fits2) / length
            sum2 = sum(x*x for x in fits2)
            std2 = abs(sum2 / length - mean2**2)**0.5
            print("  Min %s %s" % (min(fits), min(fits2)))
            print("  Max %s %s" % (max(fits), max(fits2)))
            print("  Avg Vs 1 %s" % mean)
            print("  Std Vs 1 %s" % std)
            print("  Avg Vs all %s" % mean2)
            print("  Std Vs all %s" % std2)
            germinalIndex += 1
            totalDataV1[germinalIndex].append(max(fits))
            totalDataVAll[germinalIndex].append(max(fits2))
            totalAveragesV1[germinalIndex].append(mean)
            totalAveragesVAll[germinalIndex].append(mean2)
        germinalIndex = 0
        totaltime = totaltime + cycletime

    print("-- End of (successful) evolution --")
    for pop in pops:
        best_ind = tools.selBest(pop, 1)[0]
        print("Best individual is %s, %s" %
              (best_ind, best_ind.fitness.values))
    f1 = open("BestAgainstAll.csv", "w", newline="")
    writer = csv.writer(f1)
    writer.writerows(totalDataVAll)
    f1.close()
    f1 = open("AveragesAgainst1.csv", "w", newline="")
    writer = csv.writer(f1)
    writer.writerows(totalAveragesV1)
    f1.close()
    f1 = open("AverageAgainstAll.csv", "w", newline="")
    writer = csv.writer(f1)
    writer.writerows(totalAveragesVAll)
    f1.close()
    f1 = open("BestAgainst1.csv", "w", newline="")
    writer = csv.writer(f1)
    writer.writerows(totalDataV1)
    f1.close()


if __name__ == "__main__":
    main()
