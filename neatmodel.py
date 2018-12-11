
# coding: utf-8

# In[3]:


import neatgenome
import neatneuron

# In[2]:

import math

maxOutput = 1.0
minOutput = -1.0
slippageCoeff = 0.05
def neat_paper_sigmoid(x):
    return 1.0 / (1.0 + math.pow(math.e, -4.9 * x))
def fast_sigmoid(x):
    output = x * 1.0 / (1.0 + abs(x))
    if output > maxOutput:
        return maxOutput
    if output < minOutput:
        return minOutput
    return round(output, 3)
def activation_function(x):
    # this is the sigmoid function used in the neat paper
    # return neat_paper_sigmoid(x)
    # this is fast sigmoid function
    return fast_sigmoid(x)
    # this is simple clamp, do not use
    if x > maxOutput:
        return maxOutput
    if x < minOutput:
        return minOutput
    return x

class Model(object):
    # this whole section of code needs to be thoroughly checked for errors
    # imo this is where the most fatal of them will be
    # where genome is a genomelen x genomelen size array
    # where genomelen is the length of the genome
    def __init__(self, genome, initValue):
        # print("NEW MODEL")
        self.__genome = genome
        self.__brain = []
        self.__layerMap = []
        # initValue refers to the amount of money we are willing to commit to a single trading strategy
        self.__initValue = initValue
        self.__equity = initValue
        self.__stocks = 0
        self.__lastPrice = 0.0
        # note that model networks now have their neurons divided into layers and that genes are translated into neurons
        self.__inputLayer = []
        self.__outputLayer = []
        for g in self.__genome.genes():
            if g.isEnabled() == True:
                newNeuronA = neatneuron.Neuron(g.inputNum())
                newNeuronB = neatneuron.Neuron(g.outputNum())
                if newNeuronA not in self.__brain:
                    newNeuronA.addNout(neatneuron.Connection(g.outputNum(), g.weight()))
                    self.__brain.append(newNeuronA)
                    if newNeuronA.node() in self.__genome.inputSet():
                        self.__inputLayer.append(newNeuronA)
                    if newNeuronA.node() in self.__genome.outputSet():
                        self.__outputLayer.append(newNeuronA)
                else:
                    self.__brain[self.__brain.index(newNeuronA)].addNout(neatneuron.Connection(g.outputNum(), g.weight()))
                if newNeuronB not in self.__brain:
                    self.__brain.append(newNeuronB)
                    if newNeuronB.node() in self.__genome.outputSet():
                        self.__outputLayer.append(newNeuronB)
                    if newNeuronB.node() in self.__genome.inputSet():
                        self.__inputLayer.append(newNeuronB)
        self.__layerMap.append(self.__inputLayer)
        self.separateIntoLayers()
        self.__layerMap.append(self.__outputLayer)
    def separateIntoLayers(self):
        done = False
        while done == False:
            tentativeLayer = []
            for no in self.__layerMap[-1]:
                for out in no.nouts():
                    outn = self.__brain[self.__brain.index(neatneuron.Neuron(out.output))]
                    if outn not in self.__outputLayer and outn not in tentativeLayer:
                        if outn.layer() == -1 or outn.layer() > no.layer():
                            tentativeLayer.append(outn)
                            outn.setLayer(len(self.__layerMap))
            if len(tentativeLayer) > 0:
                self.__layerMap.append(tentativeLayer)
            else:
                done = True
    def stocks(self):
        return self.__stocks
    def potential(self):
        return self.__equity
        # return self.__equity + self.__lastPrice * self.__stocks
    def resetEquity(self):
        self.__equity = self.__initValue
        self.__stocks = 0
    def genomeNodeSet(self):
        # note that genome's nodeSet is used to 
        return self.__genome.nodeSet()
    def genomeOutputSet(self):
        return self.__genome.outputSet()
    def genomeInputSet(self):
        return self.__genome.inputSet()
    def layerMap(self):
        return self.__layerMap
    
    def buy(self, value):
        self.__stocks += math.floor(self.__equity * 1.0 / (value * (1.00 + slippageCoeff)))
        self.__equity -= (value * (1.00 + slippageCoeff)) * math.floor(self.__equity * 1.0 / (value * (1.00 + slippageCoeff)))
        pass
    def sell(self, value):
        self.__equity += (value * (1.00 - slippageCoeff)) * self.__stocks
        self.__stocks = 0
        pass
       
    
    def equity(self):
        return self.__equity
    def genome(self):
        return self.__genome
    def runAlt(self, data):
        for d in data:
            # Data Point d - set of attributes that are deemed necessary for trading strategies
            for neuron in self.__brain:
                neuron.resetVal()
            self.__lastPrice = d[0]
            for i in range(0, len(d)):
                if i < len(self.__inputLayer):
                    self.__inputLayer[i].resetVal()
                    self.__inputLayer[i].passVal(d[i])
            for layer in self.__layerMap:
                for neuron in layer:
                    for output in neuron.nouts():
                        nextNeuron = self.__brain[self.__brain.index(neatneuron.Neuron(output.output))]
                        nextNeuron.passVal(activation_function(neuron.val() * output.weight))
            for neuronIndex in range(len(self.__outputLayer)):
                if self.__outputLayer[neuronIndex].val() > 1.0:
                    
                    self.buy(d[0])
                if self.__outputLayer[neuronIndex].val() < -1.0:
                    
                    self.sell(d[0])
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            for g in self.genome().genes():
                if g not in other.genome().genes():
                    return False
            for g in other.genome().genes():
                if g not in self.genome().genes():
                    return False
            return True
        return False
    def __ne__(self, other):
        if isinstance(other, self.__class__):
            for g in self.genome().genes():
                if g not in other.genome().genes():
                    return True
            for g in other.genome().genes():
                if g not in self.genome().genes():
                    return True
            return False
        return True
    def forceEquity(self, value):
        self.__equity = value


