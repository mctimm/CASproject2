import deap
import numpy
import array
import random
from deap import tools, creator, base
import bitstring
from pyrsistent import mutant
import scipy.spatial.distance
import copy

numVirus = 10
epitotesSize = 128 # the size of the epitotes and antibodies
cycletime = 8 # the time to mutate in hours
closeAntibodyPercentage = 50 #values between 0-100
diff = 25 #the difference between the antibody and the virus
def newVirusOrAntibody():
    base = []
    for i in range(0,128):
        base.append(random.randint(0,1))
    print(base)
    return base

def closeAntibody(baseVirus, diff):
    anitbody = copy.deepcopy(baseVirus[random.randint(0,len(baseVirus)-1)])
    for i in range(0,len(anitbody)):
        if random.randint(0,100) <= diff:
            anitbody[i] = int(not anitbody[i])
    print(anitbody)
    return anitbody

def createAntibodies(baseVirus, closePrecent, diff):
    if random.randint(0,100) <= closePrecent:
        return closeAntibody(baseVirus, diff)
    else: return newVirusOrAntibody()

#we create a base virus. for each antibody, we'll vary between new random numbers and tweaking this
baseVirus = []
for i in range(0,numVirus):
    baseVirus.append(newVirusOrAntibody())
#might want to add a third parameter to account for a specific virus we are targetting.
def fitnessFunction(epitotes, antibode):
    return (singleFitness(epitotes, antibode), generalFitness(epitotes, antibode))
#How well against of the virus's in our study?
def singleFitness(epitotes, antibode):
    return max(map(lambda x : epitotesSize - (scipy.spatial.distance.hamming(x,antibode) * epitotesSize),epitotes))
#How well against the whole collection?
def generalFitness(epitotes, antibode):
    return sum(map(lambda x : epitotesSize - (scipy.spatial.distance.hamming(x,antibode) * epitotesSize),epitotes))/(len(epitotes))

def randomsingleMutation(antibody):
    antibody.invert(random.randint(0,len(antibody) - 1))

creator.create("FitnessMax",base.Fitness,weights = (1.0, 0.1))
creator.create("BCell",list,fitness=creator.FitnessMax)

toolbox = base.Toolbox()

toolbox.register("attr_new_anti",createAntibodies,baseVirus,closeAntibodyPercentage, diff)
toolbox.register("individual", tools.initIterate, creator.BCell, toolbox.attr_new_anti)
toolbox.register("GerminalCenter", tools.initRepeat, list, toolbox.individual)

toolbox.register("generalEval",lambda x:generalFitness(baseVirus,x))
toolbox.register("evaluate",lambda x:fitnessFunction(baseVirus,x))
toolbox.register("mutate",tools.mutFlipBit,indpb=0.1)
toolbox.register("select",tools.selTournament, tournsize=3)


def main():
    pop = toolbox.GerminalCenter(n=300)
    print("Start of evolution")
    MUTPB = 0.2
    fitnesses = list(map(toolbox.evaluate, pop))
    for ind, fit in zip(pop, fitnesses):
        print(fit)
        ind.fitness.values = fit
    
    print("  Evaluated %i individuals" % len(pop))

    fits = [ind.fitness.values[0] for ind in pop]

    # Variable keeping track of the number of generations
    g = 0

    # Begin the evolution
    
    totaltime = 0
    
    while g < 1000:
        # A new generation
        g = g + 1
        totaltime = totaltime + cycletime
        print("-- Generation %i --" % g)

        # Select the next generation individuals
        #I'll need to figure out how it choose individuals and how we can change that.
        offspring = toolbox.select(pop, len(pop))
        # Clone the selected individuals
        offspring = list(map(toolbox.clone, offspring))

        # Apply mutation
        for mutant in offspring:
            # mutate an individual with probability MUTPB
            if random.random() < MUTPB:
                toolbox.mutate(mutant)
                del mutant.fitness.values

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

        length = len(pop)
        mean = sum(fits) / length
        sum2 = sum(x*x for x in fits)
        std = abs(sum2 / length - mean**2)**0.5

        print("  Min %s" % min(fits))
        print("  Max %s" % max(fits))
        print("  Avg %s" % mean)
        print("  Std %s" % std)

    print("-- End of (successful) evolution --")

    best_ind = tools.selBest(pop, 1)[0]
    print("Best individual is %s, %s" % (best_ind, best_ind.fitness.values))
if __name__ == "__main__":
    main()
