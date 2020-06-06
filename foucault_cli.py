from foucault_test import FoucaultTest
import json
import sys
import matplotlib.pyplot as plt

spec = {}
with open(sys.argv[1]) as f:
    spec = json.load(f)

to_meters = lambda x: x * spec['meter_conversion_factor']
from_meters = lambda x: x / spec['meter_conversion_factor']

model = FoucaultTest(to_meters(spec['apperture']), spec['f_ratio'],
                               [to_meters(i) for i in spec['outer_zone_radii']])
data = [0] * len(spec['test_data'][0])
for i in spec['test_data']:
    for j in range(len(i)):
        data[j] += i[j]
data = [i / len(spec['test_data']) for i in data]
res = model.surface_error([to_meters(i) for i in data])
plt.plot([0] + spec['outer_zone_radii'], [0] + [from_meters(i) for i in res])
plt.show()
