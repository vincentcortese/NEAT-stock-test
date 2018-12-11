'''
Created on Jul 27, 2018

@author: vmaaaaan
'''

"""data format: open, high, low, close, volume
"""
import alpha_vantage.timeseries as avts


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
ts = avts.TimeSeries(key = 'UQV2D45VFXJGX00j')
print("Fetching Data...")

raw_data, meta_data = ts.get_daily('ge', outputsize="full")
#print(raw_data)
data_p1 = str(raw_data).split("': '")
data_p2 = []
for d in data_p1:
    for ds in d.split("', '"):
        data_p2.append(ds)
data_p2b = []
for d in data_p2:
    for ds in d.split("'},"):
        data_p2b.append(ds)
data_rev = []
new_seg = []
for d in data_p2b:
    if is_number(d.strip()):
        new_seg.append(float(d.strip()))
    if len(new_seg) == 5:
        data_rev.append(new_seg)
        new_seg = []
                
data = []
for dx in range(len(data_rev)-1, -1, -1):
    data.append(data_rev[dx])
shortened = data[-1260:]

no_vol = [list(row) for row in shortened]
for x in no_vol:
    del x[4]
print(len(shortened))
print("Data Fetched")
#print(raw_data)
#print(data_p1)
#print(data_p2)
print(shortened[0:2])
print(no_vol[0:2])



