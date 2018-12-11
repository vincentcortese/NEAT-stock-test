
# coding: utf-8

# In[28]:


# import pandas as pd
# import numpy as np
import random
import neatmodel
import neatgenome
import neatgene


# In[26]:
compatibilityThreshold = 3.0
interspeciesMatingCoeff = 0.001
flipDisabledCoeff = 0.25
flipEnabledCoeff = 0.05
stagnantGenLimit = 15
preserveSpecies = True
normalizeGenome = True
disjointCoeff = 2.0
weightAvgCoeff = 2.4
survivalCoeff = 0.2
weightMutationCoeff = 0.8
totalWeightResetCoeff = 0.1
connectionMutationCoeff = 0.05
nodeMutationCoeff = 0.03
maxPerturb = 2.5
INPUT = 6
OUTPUT = 1
INIT_EQUITY = 100.0
def sortOrganisms(arr, l, r):
    # using mergeSort to sort array
    if r > l:
        m = (l + r)//2
        sortOrganisms(arr, l, m)
        sortOrganisms(arr, m+1, r)
        merge(arr, l, m, r)
def merge(arr, l, m, r):
    n1 = m-l+1
    n2 = r-m
    L = [arr[l+x] for x in range(0, n1)]
    R = [arr[m+1+x] for x in range(0, n2)]
    i = 0
    j = 0
    n = l
    while i < n1 and j < n2:
        # if L[i].equity() >= R[j].equity():
        if L[i].potential() >= R[j].potential():
            arr[n] = L[i]
            i += 1
        else:
            arr[n] = R[j]
            j += 1
        n += 1
    while i < n1:
        arr[n] = L[i]
        i += 1
        n += 1
    while j < n2:
        arr[n] = R[j]
        j += 1
        n += 1
            
class Main(object):
    def __init__(self, numIterations, numModels):
        self.__iterations = numIterations
        # unsure if this is right syntax for declaring empty arrays
        # models refers to the models themselves (the current generation of topologies)
        self.__models = []
        # history refers to every unique Gene, and each's historical marker
        self.__history = []
        # species should be a two-dimensional array, where it contains lists of species - need to recheck the paper if this 
        # approach is correct
        self.__species = []
        self.__numModels = numModels
        self.__globalInnovationNum = 0
        self.__specieshighs = []
        self.__speciesstagnantgens = []
        self.__maxOverall = 0.0
        self.__championOverall = None
        # seed the rng
        random.seed()
    def iterations(self):
        return self.__iterations
    def geneRecord(self):
        return self.__history
    def models(self):
        return self.__models
    def species(self):
        return self.__species
    def generateModels(self):
        # generate X number of uniform populations
        for a in range(0, self.__numModels):
            # write a base genome 
            newGen = neatgenome.Genome(INPUT, OUTPUT)
            # note that the number 6 is arbitrary - I have no idea how many input nodes there are 
            # because I have no idea what they should be!
            if a == 0:
                # for x in range(INPUT, INPUT + OUTPUT):
                    # self.__history.append(neatgene.Gene(x, x, 1, self.__globalInnovationNum))
                    # self.__globalInnovationNum += 1
                for x in range(0, INPUT):
                    # again 6 is the number of input nodes - is arbitrary
                    # self.__history.append(neatgene.Gene(x, x, 1, self.__globalInnovationNum))
                    # self.__globalInnovationNum += 1
                    for y in range(INPUT, INPUT + OUTPUT):
                        # I am fairly certain that we will only have two output nodes however~
                        self.__history.append(neatgene.Gene(x, y, 1, self.__globalInnovationNum))
                        self.__globalInnovationNum += 1
            self.__models.append(neatmodel.Model(newGen, INIT_EQUITY))
            # note that the number 3000 is arbitrary and will depend entirely on the amount of money we are willing to commit
    def runModels(self, data):
        # run through all models
        # change runModels(self) to runModels(self, data)
        # for DATA in DATASTREAM:
        for modelindex in range(0, len(self.__models)):
            # need to link the stock price data and pass it through the run function in proper format
            # model.run(DATA)
            # print("Running Model " + str(modelindex + 1) + " out of " + str(len(self.__models)))
            model = self.__models[modelindex]
            model.runAlt(data)
            # print(model.equity())
    
    def cleanGeneRecord(self):
        for g in self.__history:
            active = False
            for model in self.__models:
                if g in model.genome().genes():
                    active = True
            if active == False:
                self.__history.remove(g)
                
    def speciate(self):
        for s in self.__species:
            # clearing the species list except for one random genome from the previous generation
            if len(s) > 1:
                # if there are topologies in this species list
                # clear the species list except for one random genome
                x = random.randint(0, len(s) - 1)
                srep = s[x]
                newS = []
                # I'm not sure if this actually works, I think it should
                newS.append(srep)
                self.__species[self.__species.index(s)] = newS
            # print("SPECIES LENGTH(!!!): " + str(len(self.__species[0])))
        oldspeciesnum = len(self.__species)
        totalcompatDist = 0
        for m in self.__models:
            unspeciated = True
            # flag to check if m has been entered into a species or not
            for s in self.__species:
                if len(s) > 0:
                    d = 0
                    # d is the number of disjoint genes
                    w = 0
                    # w is the average weight difference between shared genes
                    i = 0
                    # i is the number of shared genes
                    # check the representative in each species list
                    for mg in m.genome().genes():
                        if mg not in s[0].genome().genes():
                            d += 1
                        else:
                            for sg in s[0].genome().genes():
                                if mg == sg:
                                    w += abs(mg.weight() - sg.weight())
                                    i += 1
                    # get the average weight difference
                    if i != 0:
                        w /= i
                    if normalizeGenome == True:
                        N = max(m.genome().glen(), s[0].genome().glen()) * 1.0
                    else:
                        N = 1.0
                    # print("d: " + str(d) + " N: " + str(N) + " w: " + str(w))
                    compatDist = disjointCoeff * d / N + weightAvgCoeff * w
                    totalcompatDist += compatDist
                    if compatDist < compatibilityThreshold:
                        # if the compatibility distance between two genomes below threshold
                        # add the model to the species
                        s.append(m)
                        unspeciated = False
                        break
            # if none of the species are compatible enough with m
            if unspeciated == True:
                newspecies = []
                newspecies.append(m)
                self.__species.append(newspecies)
                self.__specieshighs.append(0)
                self.__speciesstagnantgens.append(0)
        # print("Average Compatibility Distance: " + str(totalcompatDist / len(self.__models)))
        for sx in range(0, oldspeciesnum):
            # remove the species representatives from the last generation
            if sx >= len(self.__species):
                break
            if len(self.__species[sx]) > 1:
                self.__species[sx].pop(0)
            else:
                # if no new models fall into this species
                # remove that species
                # del self.__species[sx]
                pass
        
        # print("MODEL NUMBER IS AT: " + str(len(self.__models)))
    
    def shuffleOrganisms(self, arr):
        temp = []
        while(len(arr) > 0):
            randIndex = random.randint(0, len(arr) - 1)
            temp.append(arr.pop(randIndex))
        for x in temp:
            arr.append(x)
        temp = []
    
    
    
    #
    #def getParentAlt(self, arr):
     #   fitness = 0
      #  curr = 0
        #minimum = 0
       # for m in arr:
         #   if m.equity() < minimum:
          #      minimum = m.equity()
        #for m in arr:
         #   if minimum < 0: 
          #      fitness += m.equity() - minimum + 1
           # else:
            #    fitness += m.equity()
        #chance = random.random() * fitness
        #for m in arr:
        #    if minimum < 0:
        #        curr += m.equity() - minimum + 1
        #    else:
        #        curr += m.equity()
        #    if chance < curr:
        #        return m
    #
    
    def getParentAlt(self, arr):
        fitness = 0
        curr = 0
        minimum = 0
        for m in arr:
            if m.potential() < minimum:
                minimum = m.potential()
        for m in arr:
            if minimum < 0: 
                fitness += m.potential() - minimum + 1
            else:
                fitness += m.potential()
        chance = random.random() * fitness
        for m in arr:
            if minimum < 0:
                curr += m.potential() - minimum + 1
            else:
                curr += m.potential()
            if chance < curr:
                return m
    def getParent(self, arr):
        # print("LENGTH OF ARR: " + str(len(arr)))
        totalfitness = 0
        minimum = 0
        for m in arr:
            if m.potential() < minimum:
                minimum = m.potential()
        for m in arr:
            if minimum <= 0:
                totalfitness += m.potential() - minimum + 1
            else:
                totalfitness += m.potential()
        while True:
            for m in arr:
                chance = 0
                if minimum <= 0:
                    chance = (m.potential() - minimum + 1) * 1.0 / totalfitness
                else:
                    chance = (m.potential() * 1.0) / totalfitness
                # print("CHANCE: " + str(chance))
                if random.random() < chance:
                    return m
    def stagnantGens(self):
        return self.__speciesstagnantgens
    def evolve(self):
        # evolve the next generation: after calling runModels(), go through each model and call equity() - this is fitness function
        # when adding connections a few rules apply:
        # connections cannot cause cycles (check ultimate root node)
        # connections cannot end with input nodes
        # connections cannot start with output nodes
        # connections cannot start and end with the same node (this would make them Node Genes, not Connection Genes)
        # all mutations (genes) are stored globally (at this level)
        # all genes have a historical marker
        # when adding a new gene, check to see if gene is already present in self.__history
        # note that mutation can mean disabling a gene or reenabling a disabled gene
        
        
        
        
        self.speciate()
        # make sure that the current population of models is speciated
        
        # for g in self.__history:
            # print(str(g.inputNum()) + " -" + str(g.weight()) +"-> " + str(g.outputNum()) + " enabled: " + str(g.isEnabled()))
        # print("MODEL NUMBER IS AT: " + str(len(self.__models)))
        mutationCoeff = 1.0
        totalfitness = 0
        minfitness = 1
        speciesfitness = [0 for _ in range(0, len(self.__species))]
        activespecies = 0
        for spe in self.__species:
            if len(spe) > 0:
                activespecies += 1
        addInChampion = [False for _ in range(0, len(self.__species))]
        maxOrganism = [None for _ in range(0, len(self.__species))]
        maxFitness = [0.0 for _ in range(0, len(self.__species))]
        
        
        for sx in range(0, len(self.__species)):
            for ssx in range(0, len(self.__species[sx])):
                speciesfitness[sx] += self.__species[sx][ssx].potential()
                if maxFitness[sx] < self.__species[sx][ssx].potential():
                    maxFitness[sx] = self.__species[sx][ssx].potential()
                    maxOrganism[sx] = self.__species[sx][ssx]
            
            # calculate the fitness of a species by taking the average fitness
            if len(self.__species[sx]) != 0:
                speciesfitness[sx] /= len(self.__species[sx])
            
            if len(self.__species[sx]) > 0:
                if maxFitness[sx] > self.__specieshighs[sx]:
                    self.__specieshighs[sx] = maxFitness[sx]
                    self.__speciesstagnantgens[sx] = 0
                else:
                    self.__speciesstagnantgens[sx] += 1
                    if activespecies == 1 and self.__speciesstagnantgens[sx] >= stagnantGenLimit * 2:
                        mutationCoeff = 10.0
                    else:
                        mutationCoeff = 1.0
            activespecies = 0
            for spe in self.__species:
                if len(spe) > 0:
                    activespecies += 1
            if self.__speciesstagnantgens[sx] >= stagnantGenLimit and activespecies > 1:
                self.__species[sx] = []
                self.__speciesstagnantgens[sx] = 0
            # find the minimum fitness value
            if speciesfitness[sx] < minfitness:
                minfitness = speciesfitness[sx]
            
        activespecies = 0
        for spe in self.__species:
            if len(spe) > 0:
                activespecies += 1
        # print("THERE ARE " + str(activespecies) + " SPECIES")
        for sx in range(0, len(self.__species)):
            # if the minimum fitness value is negative
            # then shift all species fitnesses so minimum fitness value is at one
            if minfitness <= 0:
                speciesfitness[sx] -= minfitness
            # remember to calculate total fitness
            if self.__speciesstagnantgens[sx] < stagnantGenLimit or activespecies == 1:
                totalfitness += speciesfitness[sx]
            
        # calculate the new populations of each species
        speciesnewpops = [0 for _ in range(0, len(self.__species))]
        # print("TOTAL FITNESS: " + str(totalfitness))
        for x in range(0, len(self.__species)):
            # each species is assigned a number of offspring proportionate to the species fitness as a whole
            if totalfitness == 0:
                totalfitness = 1
            # print("SPECIES FITNESS: " + str(speciesfitness[x]))
            # print("SPECIES FITNESS OUT OF TOTAL FITNESS: " + str(speciesfitness[x] * 1.0 / totalfitness))
            if self.__championOverall != None:    
                speciesnewpops[x] = int(round(speciesfitness[x] * 1.0 / totalfitness * (self.__numModels - 1)))
            else:
                speciesnewpops[x] = int(round(speciesfitness[x] * 1.0 / totalfitness * self.__numModels))
            activespecies = 0
            
            # if self.__speciesstagnantgens[x] >= stagnantGenLimit and activespecies > 1:
                # speciesnewpops[x] = 0
            if len(self.__species) == 1:
                speciesnewpops[x] = max(self.__numModels, speciesnewpops[x])
            if len(self.__species[x]) >= 5:
                speciesnewpops[x] -= 1
                addInChampion[x] = True
        if self.__championOverall != None:
            # print("Overall Champion: (" + str(self.__maxOverall) + ")")
            # for g in self.__championOverall.genome().genes():
                # print(str(g.inputNum()) + " -[" + str(g.weight()) +"]-> " + str(g.outputNum()) + " enabled: " + str(g.isEnabled()))
            if self.__championOverall not in self.__models:
                self.__championOverall.resetEquity()
                self.__models.append(self.__championOverall)
        for s in self.__species:
            specieshighscore = 0
            for m in s:
                if specieshighscore < m.potential():
                    specieshighscore = m.potential()
            # if len(s) > 0:
                # print("Species " + str(self.__species.index(s)) + " Population: " + str(len(s)) + " with Highest Fitness: " + str(specieshighscore) + " and Average Fitness: " + str(speciesfitness[self.__species.index(s)]) + " (anticipating "+str(speciesnewpops[self.__species.index(s)] + (addInChampion[self.__species.index(s)])) +" members in next generation) at " + str(self.__speciesstagnantgens[self.__species.index(s)]) + " Stagnant Generations")
                # if addInChampion[self.__species.index(s)] == True:
                    # if maxOrganism[self.__species.index(s)] != None:
                        # print("Species " + str(self.__species.index(s)) + " Champion: (" + str(maxOrganism[self.__species.index(s)].potential()) + ")")
                        # for g in maxOrganism[self.__species.index(s)].genome().genes():
                            # print(str(g.inputNum()) + " -[" + str(g.weight()) +"]-> " + str(g.outputNum()) + " enabled: " + str(g.isEnabled()))
                # randomSpecimen = random.randint(0, len(s) - 1)
                # print("Example of Species " + str(self.__species.index(s)) + " (" + str(s[randomSpecimen].equity()) + ")")
                # for g in s[randomSpecimen].genome().genes():
                    # print(str(g.inputNum()) + " -[" + str(g.weight()) +"]-> " + str(g.outputNum()) + " enabled: " + str(g.isEnabled()))
        # sorting species
        # for s in self.__species:
        #    self.shuffleOrganisms(s)
        for s in self.__species:
            # only keep the fittest organisms
            if len(s) > 0:
                sortOrganisms(s, 0, len(s)-1)
                # output = ""
                # for m in s:
                    # output += str(m.equity()) + " "
                # print(output)
                minimumSurvivalNumber = 1
                if preserveSpecies == False:
                    minimumSurvivalNumber = 0
                keepTill = int(max(minimumSurvivalNumber, survivalCoeff * len(s)))
                s = s[0:keepTill]
                if keepTill == 0:
                    s = []
            # print("MODEL NUMBER IS AT: " + str(len(self.__models)))
        for m in self.__models:
            if m.potential() > self.__maxOverall:
                self.__maxOverall = m.potential()
                self.__championOverall = m
        self.__models = []
        # print("MODEL NUMBER IS AT: " + str(len(self.__models)))
        # print("START OF REPOPULATION")
        for i in range(0, len(self.__species)):
            # print("SPECIES NEW POP: " + str(speciesnewpops[i]))
            if len(self.__species[i]) > 0:
                for j in range(0, speciesnewpops[i]):
                    # print("Generating Child " + str(j + 1) + " out of " + str(speciesnewpops[i]))
                    # breed from surviving organisms
                    newGenome = neatgenome.Genome(INPUT, OUTPUT)
                    # note that parents are selected from species lists
                    # print("GETTING PARENTS")
                    father = self.getParent(self.__species[i])
                    mother = self.getParent(self.__species[i])
                    # print("PARENTS GOTTEN")
                    # there is a small chance that the other parent will be of a different species
                    activespecies = 0
                    for spe in self.__species:
                        if len(spe) > 0:
                            activespecies += 1
                    if activespecies > 1 and random.random() < interspeciesMatingCoeff:
                        motherSpeciesIndex = random.randint(0, len(self.__species) - 1)
                        # print("INTERSPECIES BREEDING")
                        # print("i: " + str(i) + " motherSpeciesIndex: " + str(motherSpeciesIndex) + " len(self.__species[motherSpeciesIndex]): " + str(len(self.__species[motherSpeciesIndex])))
                        while motherSpeciesIndex == i or len(self.__species[motherSpeciesIndex]) == 0:
                            # print("motherSpeciesIndex: " + str(motherSpeciesIndex) + " len(self.__species[motherSpeciesIndex]): " + str(len(self.__species[motherSpeciesIndex])))
                            motherSpeciesIndex = random.randint(0, len(self.__species) - 1)
                        mother = self.getParent(self.__species[motherSpeciesIndex])
                        # print("INTERSPECIES PARENT GOTTEN")
                    # note that the parents can be the same organism
                    # print("BREEDING PARENTS")
                    if father.potential() == mother.potential():
                        # if parent fitness values are the same
                        # the loops below are very inefficient
                        for fg in father.genome().genes():
                            # checking each gene in father's genome
                            if fg not in mother.genome().genes():
                                # disjoint genes are inherited randomly in equal fitness parents
                                if random.random() < 0.5:
                                    # 50-50 whether or not it is inherited
                                    newGenome.addGene(fg.inputNum(), fg.outputNum(), fg.weight(), fg.historical())    
                            for mg in mother.genome().genes():
                                # checking each gene in mother's genome
                                if mg not in father.genome().genes():
                                    if random.random() < 0.5:
                                        newGenome.addGene(mg.inputNum(), mg.outputNum(), mg.weight(), mg.historical())
                                if fg == mg:
                                    if random.random() < 0.5:
                                        newGenome.addGene(fg.inputNum(), fg.outputNum(), fg.weight(), fg.historical())
                                    else:
                                        newGenome.addGene(mg.inputNum(), mg.outputNum(), mg.weight(), mg.historical())
                    else:
                        # for parents with unequal fitness values
                        if father.potential() > mother.potential():
                            # if father is more fit
                            for fg in father.genome().genes():
                                newGenome.addGene(fg.inputNum(), fg.outputNum(), fg.weight(), fg.historical())
                        else:
                            # if mother is more fit
                            for mg in mother.genome().genes():
                                newGenome.addGene(mg.inputNum(), mg.outputNum(), mg.weight(), mg.historical())
                    for ng in newGenome.genes():
                        if ng.isEnabled() == False:
                            if random.random() < flipDisabledCoeff:
                                ng.setEnabled(True)
                        else:
                            if random.random() < flipEnabledCoeff * mutationCoeff:
                                ng.setEnabled(False)
                    # connection mutation
                    inputLayerSize = len(newGenome.inputSet()) 
                    outputLayerSize = len(newGenome.outputSet())
                    hiddenLayerSize = len(newGenome.nodeSet()) - len(newGenome.inputSet()) - len(newGenome.outputSet())
                    startingConnections = inputLayerSize * outputLayerSize
                    hiddenConnections = hiddenLayerSize * (hiddenLayerSize - 1)
                    connectionsFromInputs = inputLayerSize * hiddenLayerSize
                    connectionsToOutputs = outputLayerSize * hiddenLayerSize
                    possibleConnections = startingConnections + hiddenConnections + connectionsFromInputs + connectionsToOutputs
                    remainingConnections = possibleConnections - len(newGenome.genes())
                    # print("CONNECTION MUTATION CODE")
                    if remainingConnections > 0 and random.random() < connectionMutationCoeff:
                        # print("CONNECTION MUTATION TRIGGERED")
                        newGene = None
                        inGenome = True
                        # flag to determine whether or not new gene mutated is really a new gene
                        n = 0
                        # absolute limit to retries
                        # print("INSIDE CONNECTION CODE")
                        while inGenome == True and n < remainingConnections:
                            inGenome = False
                            # print(n)
                            n += 1
                            newInput = newGenome.nodeSet()[random.randint(0, len(newGenome.nodeSet()) - 1)]
                            # while newInput in newGenome.outputSet():
                                # make sure that the newInput is not one of the outputs
                                # newInput = newGenome.nodeSet()[random.randint(0, len(newGenome.nodeSet()) - 1)]
                            newOutput = newGenome.nodeSet()[random.randint(newGenome.inputSet()[-1], len(newGenome.nodeSet()) - 1)]
                            while newOutput in newGenome.inputSet() or (newInput in newGenome.outputSet() and newOutput in newGenome.outputSet()):
                                newInput = newGenome.nodeSet()[random.randint(0, len(newGenome.nodeSet()) - 1)]
                                newOutput = newGenome.nodeSet()[random.randint(newGenome.inputSet()[-1], len(newGenome.nodeSet()) - 1)]
                            
                            newGene = neatgene.Gene(newInput, newOutput, round(random.random() * 2.0 - 1.0, 3), self.__globalInnovationNum)
                        
                            if newGene in newGenome.genes():
                                inGenome = True
                            else:
                                if newGene not in self.__history:
                                    self.__history.append(newGene)
                                    self.__globalInnovationNum += 1
                                newGenome.addGene(newGene.inputNum(), newGene.outputNum(), newGene.weight(), newGene.historical()) 
                                break
                        # print("FINISH CONNECTION MUTATION")
                    else:
                        # node mutation
                        # print("NODE MUTATION CODE")
                        if random.random() < nodeMutationCoeff * mutationCoeff:
                            # print("NODE MUTATION TRIGGERED")
                            inGenome = True
                            # flag to determine whether or not new gene mutatedis really a new gene
                            n = 0
                            # absolute limit to retries
                            while inGenome == True:
                                inGenome = False
                                # print(n)
                                n += 1
                                splicedGene = newGenome.genes()[random.randint(0, len(newGenome.genes()) - 1)]
                                # while splicedGene.inputNum() == splicedGene.outputNum():
                                #    splicedGene = newGenome.genes()[random.randint(0, len(newGenome.genes()) - 1)]
                                # here we may face some issues down the road depending on the node number we choose for the new node
                                # choosing new node number based on min node num not already in or surpassed by the node numbers in
                                # the topology
                                newNodeNum = -1
                                for g in newGenome.genes():
                                    # find the max number in newGenome.nodeSet()
                                    if newNodeNum < g.inputNum():
                                        newNodeNum = g.inputNum()
                                    if newNodeNum < g.outputNum():
                                        newNodeNum = g.outputNum()
                                # increment the new node number to get a new number
                                newNodeNum += 1
                                # make new genes by splitting an existing connection into two
                                newGeneA = neatgene.Gene(splicedGene.inputNum(), newNodeNum, 1.0, self.__globalInnovationNum)
                                newGeneB = neatgene.Gene(newNodeNum, splicedGene.outputNum(), splicedGene.weight(), self.__globalInnovationNum + 1)
                                if newGeneA in newGenome.genes() or newGeneB in newGenome.genes():
                                    inGenome = True
                                else:
                                    # if the genes are unique then add them in
                                    splicedGene.disable()
                                    newGenome.addGene(newGeneA.inputNum(), newGeneA.outputNum(), newGeneA.weight(), newGeneA.historical())
                                    newGenome.addGene(newGeneB.inputNum(), newGeneB.outputNum(), newGeneB.weight(), newGeneB.historical())
                                    if newGeneA not in self.__history and newGeneB not in self.__history:
                                        self.__history.append(newGeneA)
                                        self.__globalInnovationNum += 1
                                        self.__history.append(newGeneB)
                                        self.__globalInnovationNum += 1
                                        break
                            # print("NODE MUTATION DONE")
                        else:
                            # weight coefficient mutation - note that this only changes one gene's weighting!
                            # print("WEIGHT MUTATION CODE")
                            for g in newGenome.genes():
                                if random.random() < weightMutationCoeff: 
                                    # print("WEIGHT MUTATION TRIGGERED")
                                    oldg = g
                                    newGenome.genes().remove(g)
                                    if random.random() < totalWeightResetCoeff:
                                        newGenome.addGene(oldg.inputNum(), oldg.outputNum(), round(random.random() * 2 - 1.0, 3), oldg.historical())
                                    else:
                                        perturbance = round(random.random() * maxPerturb * 2.0 - maxPerturb, 3) + oldg.weight()
                                        if perturbance > 1.0:
                                            perturbance = 1.0
                                        if perturbance < -1.0:
                                            perturbance = -1.0
                                        newGenome.addGene(oldg.inputNum(), oldg.outputNum(), perturbance, oldg.historical())
                            # print("WEIGHT MUTATION DONE")
                    # add the model to the population
                    # print("ADDING IN MODEL TO POPULATION")
                    # for g in newGenome.genes():
                        # print(str(g.inputNum()) + " -[" + str(g.weight()) +"]-> " + str(g.outputNum()) + " enabled: " + str(g.isEnabled()))
                    newModel = neatmodel.Model(newGenome, INIT_EQUITY)
                    self.__models.append(newModel)
                    # print("MODEL ADDED IN")
            if addInChampion[i] == True:
                # print("ADDING IN CHAMPION")
                # for species with 5 or more members, add in the most successful organism unchanged
                
                # print("Species " + str(i) + " Champion (" + str(maxOrganism[i].equity()) + ")")
                if maxOrganism[i] != None:
                    maxOrganism[i].resetEquity()
                # for g in maxOrganism[i].genome().genes():
                    # print(str(g.inputNum()) + " -[" + str(g.weight()) +"]-> " + str(g.outputNum()) + " enabled: " + str(g.isEnabled()))
                    self.__models.append(maxOrganism[i])
                # print("CHAMPION ADDED IN")
        if self.__championOverall != None:
            #print("Overall Champion (" + str(self.__maxOverall) + ")")
            self.__championOverall.resetEquity()
            
            
            self.__models.append(self.__championOverall)   
        
        # print("TOTAL POPULATION: " + str(len(self.__models)))

