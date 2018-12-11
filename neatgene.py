
# coding: utf-8

# In[2]:


class Gene(object):
    def __init__(self, inputNum, outputNum, weight, historical):
        # inputNum refers to the input node number, outputNum refers to the output node number, weight
        # refers to the weight of the connection (if it is a connection), islink refers to whether or not
        # it is a connection, and historical is the historical marker
        self.__historical = historical
        self.__nodeNum = inputNum
        self.__inputNum = inputNum
        self.__outputNum = outputNum
        self.__weight = weight
        self.__enabled = True
    def __eq__(self, other):
        # overrides '==' to make equality/comparison an easier deal when deciding to assign a new historical marker and
        # whether to add to history or not
        if isinstance(other, self.__class__):
            return self.__inputNum == other.__inputNum and self.__outputNum == other.__outputNum
        return False
    def __ne__(self, other):
        # see the comment on __eq__
        return self.__inputNum != other.__inputNum or self.__outputNum != other.__outputNum
    def inputNum(self):
        # returns the input node number
        return self.__inputNum
    def outputNum(self):
        # returns the output node number
        return self.__outputNum
    def historical(self):
        # returns the historical marker
        return self.__historical
    def setEnabled(self, value):
        self.__enabled = value
    def weight(self):
        # returns the weight of the connection, note that this is always 1 if the Gene is a node not a link
        return self.__weight
    def isEnabled(self):
        # returns whether or not this Gene is enabled
        return self.__enabled
    def disable(self):
        # this will disable the Gene
        self.__enabled = False
    def reenable(self):
        # this will enable the Gene
        self.__enabled = True
