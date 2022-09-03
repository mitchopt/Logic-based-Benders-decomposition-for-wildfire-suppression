# -*- coding: utf-8 -*-
"""
Instance parameters for the experiments

"""

def ParametersSmall():
    Parameters = []
    for n1 in range(8):
        for n2 in range(10):
            Parameters.append(("10", n1, n2))
    for n1 in range(8, 16):
        for n2 in range(10):
            Parameters.append(("20", n1, n2))
    for n1 in range(16, 24):
        for n2 in range(10):
            Parameters.append(("30", n1, n2))
    return Parameters

def ParametersLarge():
    Parameters = []
    for Size in ["a", "b"]:
        for n in range(8):
            Parameters.append((Size, n))
    return Parameters
            
