from collections import namedtuple
Connection = namedtuple("Connection", "output weight")
class Neuron(object):
    def __init__(self, node):
        self.__node = node
        self.__nouts = []
        self.__layer = -1
        self.__val = 0.0
    def addNout(self, value):
        isin = False
        for nout in self.__nouts:
            if value.output == nout.output:
                isin = True
        if isin == False:
            self.__nouts.append(value)
    def setLayer(self, value):
        self.__layer = value
    def passVal(self, value):
        self.__val += value
    def resetVal(self):
        self.__val = 0.0
    def val(self):
        return self.__val
    def layer(self):
        return self.__layer
    def node(self):
        return self.__node
    def nouts(self):
        return self.__nouts
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if self.__node == other.node():
                return True
        return False
    def __ne__(self, other):
        if isinstance(other, self.__class__):
            return False
        if other.node() == self.__node:
            return False
        return True
