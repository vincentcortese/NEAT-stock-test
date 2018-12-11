import csv

'''Data dates from 9:30 am on June, 8, 2018 with every minute up until
June, 21 till 2:39 pm. Before market close

This will give a lot of specific data for micro trends, I am considering
moving averages and comparing daily data during this range. 
'''
data = []
with open('GE_daily.csv', 'r') as GE_15:
    GE = csv.reader(GE_15, quoting=csv.QUOTE_NONNUMERIC)
    for line in GE:
        data.append(line)
# print('[Open, High, Low, Close, Volume]')
# print(data[0:5])
