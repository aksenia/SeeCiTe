#!/usr/bin/env python3
'''
General utility functions concerning dictionaries, lists etc.
'''

def get_range(d, begin, end):
    '''
    Subset a dictionary, returning the elements in the given index range(treat dict as array)
    Args:
        d: python dictionary
        begin: first key element
        end: end key element

    Returns:
        subdictionary with keys 'ranging' from 'begin' to 'end'

    '''
    return dict(e for i, e in enumerate(d.items()) if begin <= i <= end)


def flatten_list_of_lists(nested_list):
    '''
    Utility function. Flatten list of lists
    :param nested_list: list with 1-level deep nested structure
    :return: flat list
    '''
    result =[]
    for sublist in nested_list:
        for item in sublist:
            result.append(item)
    return(result)


def write_list(filename, listvar):
    '''
    Write list to a file, one word per row
    Args:
        filename: name of the file out
        listvar: list of strings

    Returns: NA, writes out a file

    '''
    out = open(filename, "w")
    for s in listvar:
        out.write("%s\n" % s)
    out.close()
