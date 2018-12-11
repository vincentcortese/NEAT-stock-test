import math
import random

import alfa_vantage
import neatmain


data = alfa_vantage.data
DATA_LENGTH = len(data)
IN_SAMPLE = 90
OUT_SAMPLE = 45
OPTIMIZE_GENS = 5
MODELNUM = 1000

print("Test Started (Daily Stock Data)")
neat = neatmain.Main(int(math.floor((DATA_LENGTH - IN_SAMPLE) * 1.0/(OUT_SAMPLE))), MODELNUM)
neat.generateModels()
print("Generated Models Successfully")
for i in range(neat.iterations()):
    # print("ITERATION " + str(i + 1) + " OUT OF " + str(neat.iterations()))
    for j in range(OPTIMIZE_GENS):
        print("GENERATION " + str(j + 1) + " OUT OF " + str(OPTIMIZE_GENS) + " ON ITERATION " + str(i + 1) + " OUT OF " + str(neat.iterations()))
        for m in neat.models(): 
            m.resetEquity()
        neat.runModels(data[i * OUT_SAMPLE:i * OUT_SAMPLE + IN_SAMPLE])
        neat.cleanGeneRecord()
        neat.evolve()
    for m in neat.models():
        m.resetEquity()
    # print("OUT OF SAMPLE RUN")
    neat.runModels(data[((i + 1) * IN_SAMPLE):((i + 1) * IN_SAMPLE + OUT_SAMPLE)])
    maxFitness = 0
    averageFitness = 0
    maxEquity = 0
    averageEquity = 0
    maxStocks = 0
    averageStocks = 0
    realMaxEquity = 0
    for m in neat.models():
        averageFitness += m.potential()
        averageEquity += m.equity()
        averageStocks += m.stocks()
        if maxFitness < m.potential():
            maxFitness = m.potential()
            maxEquity = m.equity()
            maxStocks = m.stocks()
        if realMaxEquity < m.equity():
            realMaxEquity = m.equity()
    averageFitness = averageFitness * 1.0 / len(neat.models())
    averageEquity = averageEquity * 1.0 / len(neat.models())
    averageStocks = averageStocks * 1.0 / len(neat.models())
    print("Max Fitness for Iteration " + str(i + 1) + ": " + str(maxFitness) + " with Equity: " + str(maxEquity) + " and Stocks: " + str(maxStocks))
    print("Average Fitness for Iteration " + str(i + 1) + ": " + str(averageFitness) + " and Average Equity: " + str(averageEquity) + " and Average Stocks: " + str(averageStocks))
    print("Max Equity for Iteration " + str(i + 1) + ": " + str(realMaxEquity))
    randm = random.randint(0, len(neat.models()) - 1)
    # print("Genome for Model " + str(randm) + ": ")
    # for g in neat.models()[randm].genome().genes():
        # print(str(g.inputNum()) + " -" + str(g.weight()) +"-> " + str(g.outputNum()) + " enabled: " + str(g.isEnabled()))

for m in neat.models():
    m.resetEquity()
neat.runModels(data[0:DATA_LENGTH])
maxFitness = 0
averageFitness = 0
maxEquity = 0
averageEquity = 0
maxStocks = 0
averageStocks = 0
realMaxEquity = 0
for m in neat.models():
    averageFitness += m.potential()
    averageEquity += m.equity()
    averageStocks += m.stocks()
    if maxFitness < m.potential():
        maxFitness = m.potential()
        maxEquity = m.equity()
        maxStocks = m.stocks()
    if realMaxEquity < m.equity():
        realMaxEquity = m.equity()
averageFitness = averageFitness * 1.0 / len(neat.models())
averageEquity = averageEquity * 1.0 / len(neat.models())
averageStocks = averageStocks * 1.0 / len(neat.models())
print("Max Fitness: " + str(maxFitness) + " with Equity: " + str(maxEquity) + " and Stocks: " + str(maxStocks))
print("Average Fitness: " + str(averageFitness) + " and Average Equity: " + str(averageEquity) + " and Average Stocks: " + str(averageStocks))
print("Max Equity: " + str(realMaxEquity))
print("TEST COMPLETED SUCCESSFULLY") 
    