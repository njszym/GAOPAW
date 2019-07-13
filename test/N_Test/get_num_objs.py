import json
from types import SimpleNamespace


def main():
    """
    Parse input.json and return total number of objective functions
    """
    with open('input.json') as input:
        input_settings = json.load(input,object_hook=lambda d: SimpleNamespace(**d))
    cmpd_list = input_settings.compounds
    element_list = []
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        element_list.extend(parse_elems(formula))
    element_list = unique(element_list)
    num_elems = len(element_list)
    if num_elems == 1 and len(cmpd.keys()) == 1:
        if elem in ['N','P']:
            print('\n'+'Number of objective functions: 2\n')
        else:
            print('\n'+'Number of objective functions: 3\n')
        return
    cmpd_diff_dict = form_cmpd_dict(cmpd_list)
    num_properties = 0
    for cmpd in cmpd_diff_dict.keys():
        for lat_type in cmpd_diff_dict[cmpd].keys():
            for property in cmpd_diff_dict[cmpd][lat_type].keys():
                num_properties += 1
    num_obj_fns = 0
    for elem in element_list:
        if elem in ['N','P']:
            num_obj_fns += 2
        else:
            num_obj_fns += 3
    num_obj_fns += num_properties
    print('\n'+'Number of objective functions: '+str(num_obj_fns)+'\n')


def form_cmpd_dict(cmpd_list):
    """
    Constructs empty dictionary of correct length for testing
    """
    cmpd_diff_dict = {}
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        lat_type = cmpd['lattice_type']
        property_list = [property for property in cmpd if property not in ['formula','lattice_type']]
        cmpd_diff_dict[formula] = {}
        cmpd_diff_dict[formula][lat_type] = {}
        for property in property_list:
            cmpd_diff_dict[formula][lat_type][property] = {}
    return cmpd_diff_dict

def test_cmpd_list(cmpd_list,cmpd_diff_dict,cmpd_template_dir):
    """
    Perform property tests for all compounds
    """
    for cmpd in cmpd_list:
        cmpd = cmpd.__dict__
        formula = cmpd['formula']
        lat_type = cmpd['lattice_type']
        property_list = [property for property in cmpd if property not in ['formula','lattice_type']]
        for property in property_list:
            ae_value = cmpd[property]
            cmpd_diff_dict[formula][lat_type][property], error_check = test_property(formula,lat_type,property,ae_value,cmpd_template_dir)
            if error_check:
                bad_run(cmpd_diff_dict)
                return cmpd_diff_dict, True
    return cmpd_diff_dict, False

def unique(value_list): 
    """
    Get list of unique elements to be tested
    """
    try: ## if list of numbers
        value_list = [round(float(value),3) for value in value_list]
    except ValueError: ## if list of strings
        list = [str(value) for value in value_list]
    unique_list = []
    for value in value_list:
        if value not in unique_list:
            unique_list.append(value)
    return unique_list

def parse_elems(formula):
    """
    Parse compound formula to obtain constituent elements
    """
    letters_only = ''.join([letter for letter in formula if not letter.isdigit()])
    index = -1
    elems = []
    for letter in letters_only:
        if letter.isupper():
            elems.append(letter)
            index += 1
        else:
            elems[index] += letter

    return elems

if __name__=='__main__':
    main()

