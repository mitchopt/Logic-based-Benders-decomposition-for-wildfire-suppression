# -*- coding: utf-8 -*-
"""
Solve instance(s) using a MIP

"""

# Packages
from parameters import ParametersSmall, ParametersLarge
from model_MIP import FireMIP
from ast import literal_eval
from pathlib import Path
import json
import sys
import csv


def Load(Instance):
    with open(Instance, 'r') as file:
        Data = json.load(file)
        Delay = Data["Delay"]
        ArrivalTimeTarget = Data["ArrivalTimeTarget"]
        ResAtTime = {literal_eval(t): Data["ResAtTime"][t] for t in Data["ResAtTime"]}
        Ignitions = [tuple(n) for n in Data["Ignitions"]]
        N = [tuple(n) for n in Data["Nodes"]]
        A = {literal_eval(a): Data["Arcs"][a] for a in Data["Arcs"]}
    
    return N, A, Ignitions, Delay, ArrivalTimeTarget, ResAtTime


# Parameters
Folder = "small"
l = "S"
Parameters = ParametersSmall()
(Size, n1, n2) = Parameters[int(sys.argv[1])]
Inst = f"instances/{Folder}/{Size}/{l}{n1}_{n2}.json"
    
# Instance location
Instance = Path(__file__).parent / Inst

# Solution location
SolutionFile = Path(__file__).parent / f"solutions/{Folder}/{Size}/Sol_{l}{n1}.csv"


# Add header
if not Path.exists(SolutionFile):
    with open(SolutionFile, 'a', newline='') as file:
        writer = csv.writer(file)
        Row = ["Folder", "Size", "ID", "Rep", "Nodes", "Arcs", "Method", "Obj val", 
                "Obj Bound", "Time (s)", "Opt Cuts", "Feas Cuts"]
        writer.writerow(Row)


# Retrieve instance
N, A, Ignitions, Delay, ArrivalTimeTarget, ResAtTime = Load(Instance)

# Gurobi parameters
TimeLimit = 7200

# Seed
GurobiSeed = 0

# Solve with greedy LBBD
ModelMIP = FireMIP(N, A, ResAtTime, Ignitions, Delay, ArrivalTimeTarget, TimeLimit, GurobiSeed)


# Store greedy solution info
Row = [Folder, Size, n1, n2, len(N), len(A), "MIP"]
Row.append(round(ModelMIP.objVal))
Row.append(round(ModelMIP.ObjBound, 2))
Row.append(round(ModelMIP.RunTime, 2))
Row.append("null")
Row.append("null")


with open(SolutionFile, 'a', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(Row)

