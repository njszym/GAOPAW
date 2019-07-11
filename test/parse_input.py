import json
from types import SimpleNamespace


with open('input.json') as input:
    input_settings = json.load(input,object_hook=lambda d: SimpleNamespace(**d))

cmpd_list = input_settings.compounds

for cmpd in cmpd_list:
    cmpd = cmpd.__dict__
    formula = cmpd['formula']
    property_list = [property for property in cmpd.keys() if property != 'formula']
    for property in property_list:
