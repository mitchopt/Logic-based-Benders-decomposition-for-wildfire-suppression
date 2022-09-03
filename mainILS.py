# -*- coding: utf-8 -*-
"""
Solve instance(s) using the ILS

"""

# Packages
from parameters import ParametersSmall, ParametersLarge
from IteratedLocalSearch import FireILS
from ast import literal_eval
from pathlib import Path
import json
import time
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
(Size, n1, n2) = Parameters[int(sys.argv[0])]
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

# ILS Parameters
MultiStarts = 50
MaxCandidates = 5
p1 = 0.075
p2 = 0.025
MaxNeighbours = 20
MaxModifications = 5
MaxFailures = 100
MaxNoImprovements = 50


# Solve
StartTime = time.time()
ZSol, ObjVal = FireILS(
    N, A, ResAtTime, Ignitions, Delay, ArrivalTimeTarget,
    MultiStarts, p1, p2, MaxNeighbours, MaxModifications, 
    MaxFailures, MaxNoImprovements, MaxCandidates)
Runtime = time.time() - StartTime


# Store ILS solution info
Row = [Folder, Size, n1, n2, len(N), len(A), "ILS"]
Row.append(round(ObjVal))
Row.append("null")
Row.append(round(Runtime, 2))
Row.append("null")
Row.append("null")


with open(SolutionFile, 'a', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(Row)




