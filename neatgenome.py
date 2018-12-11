
# coding: utf-8

# In[ ]:


import neatgene


# In[1]:


class Genome(object):
    def __init__(self, inputNum, outputNum):
        # nodes is the array containing all the Genes
        self.__nodes = []
        self.__inputNums = inputNum
        self.__outputNums = outputNum
        self.__inputSet = []
        self.__outputSet = []
        self.__nodeSet = []
        for i in range(0, inputNum + outputNum):
            # need a link from each node to itself to ensure that everything gets mapped properly
            # note that these links are functionally useless and should be treated as such
            # weight has to be one
            if i < inputNum:
                self.__inputSet.append(i)
            else:
                self.__outputSet.append(i)
            self.__nodeSet.append(i)
            # self.__nodes[i] = neatgene.Gene(i, i, 1, -1)
        # attach links between every input node and every output node
        for x in range(0, inputNum):
            for y in range(inputNum, inputNum + outputNum):
                self.__nodes.append(neatgene.Gene(x, y, 1.0, (y - inputNum) + x * 2))
    def inputSet(self):
        return self.__inputSet
    def outputSet(self):
        return self.__outputSet
    def nodeSet(self):
        return self.__nodeSet
    def inputs(self):
        return self.__inputNums
    def outputs(self):
        return self.__outputNums
    def glen(self):
        # returns the number of nodes
        return len(self.__nodeSet)
    def genes(self):
        # return the array of Genes
        return self.__nodes
    def addGene(self, iNum, oNum, weight, historical):
        # appends another Gene to the array
        tempg = neatgene.Gene(iNum, oNum, weight, historical)
        if tempg not in self.__nodes:
            if iNum not in self.__nodeSet:
                self.__nodeSet.append(iNum)
            if oNum not in self.__nodeSet:
                self.__nodeSet.append(oNum)
            self.__nodes.append(neatgene.Gene(iNum, oNum, weight, historical))
