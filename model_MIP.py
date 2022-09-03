# -*- coding: utf-8 -*-
"""
Main MIP model

"""

# Packages
from shortest_paths import ShortestPaths
import gurobipy as gp


# ----------------------- #
# --- MIP Formulation --- #
# ----------------------- #
def FireMIP(N, A, ResAtTime, Ignitions, Delay, ArrivalTimeTarget, TimeLimit, GurobiSeed):

    
    # Epsilon
    EPS = 0.0001
    
    # Generate in- and outarcs
    InArcs = {n: [] for n in N}
    OutArcs = {n: [] for n in N}
    for a in A:
        OutArcs[a[0]].append((a[0], a[1]))
        InArcs[a[1]].append((a[0], a[1]))
        
    # Time periods
    T = list(ResAtTime.keys())
    

    # Time periods + target
    T1 = T + [ArrivalTimeTarget]


    # Retrieve the root node
    assert len(Ignitions) == 1
    RootNode = Ignitions[0]
    
    
    # Model
    Model = gp.Model()
    Model.setParam("Threads", 1)
    Model.setParam("TimeLimit", TimeLimit)
    
    
    # Seed if given
    if type(GurobiSeed) is int:
        Model.setParam("Seed", GurobiSeed)
            
    
    # Flow variables
    # X[a] the flow on arc a
    X = {a: Model.addVar() for a in A}

    # Dual variables
    # Lambda[n] is the arrival time of fire at node n
    Lambda = {n: Model.addVar(lb=-gp.GRB.INFINITY) for n in N if n != RootNode}
    Lambda[RootNode] = Model.addVar(lb=0, ub=0)
    
    # Slack variables 
    S = {a: Model.addVar() for a in A}

    # Resource allocation
    Z = {(n, t): Model.addVar(
        vtype=gp.GRB.BINARY) for n in N for t in T}
    
    # Binary variables
    # Q[a] = 1 if arc a appears in the tree
    Q = {a: Model.addVar(vtype=gp.GRB.BINARY) for a in A}
    for a in A:
        Model.addConstr(X[a] <= (len(N) - 1)*Q[a])
        
    # Indicator variables
    # Y[n, t] = 1 if node n has burned by time t
    Y = {(n, t): Model.addVar(vtype=gp.GRB.BINARY) for n in N for t in T1}
    
    
    # Ensure len(Nodes) - 1 paths leave the root node
    RootCon = Model.addConstr(gp.quicksum(
        X[a] for a in OutArcs[RootNode]) == len(N) - 1)
    
    # Exactly one path ends at every other node
    FlowConservation = {n: Model.addConstr(
        gp.quicksum(X[a] for a in InArcs[n])
        - gp.quicksum(X[a] for a in OutArcs[n]) == 1)
        for n in N if n != RootNode}
    
    # Dual constraints
    DualConstraint = {a: Model.addConstr(
        Lambda[a[1]] - Lambda[a[0]] + S[a] == A[a] + Delay * gp.quicksum(
            Z[a[0], t] for t in T)) for a in A}

    # Calculate Big M
    BigM = (len(N) - 1)*max(A[a] for a in A) + \
        (sum(ResAtTime[t] for t in ResAtTime) - 1)*Delay + EPS

    # Bound the slack on
    for a in A:  # arcs that don't belong to the tree
        Model.addConstr(S[a] <= BigM*(1 - Q[a]))

    # At most one resource per node
    AtMostOneResPernode = {n: Model.addConstr(
        gp.quicksum(Z[n, t] for t in T) <= 1) for n in N}
    
    # Up to ResAtTime[t] resources at time t
    MaxRes = {t: Model.addConstr(gp.quicksum(
        Z[n, t] for n in N) <= ResAtTime[t]) for t in T}

    # Only put resources on unburned nodes
    ResOnlyIfNotBurned = {(n, t): Model.addConstr(
        Z[n, t] <= 1 + (Lambda[n] - t)/t ) for n in N for t in T}

    # Lambdas force the Y variables
    DoesBurn = {(n, t): Model.addConstr(
        Y[n, t] >= (t - Lambda[n])/t) for n in N for t in T1}
            
    # Objective is the number of nodes burned by the target time
    Model.setObjective(gp.quicksum(
        Y[n, ArrivalTimeTarget] for n in N), gp.GRB.MINIMIZE)

    # Solve problem
    Model.optimize()
    
    
    # Save results
    Model._Optimal = set(
        n for n in N for t in T if Z[n, t].x > .1)
    Model._ArrivalTime, Model._FirePath, Model._Pred = ShortestPaths(
        N, A, InArcs, OutArcs, Model._Optimal, Ignitions, Delay)
    Model._Burned = {n: Model._ArrivalTime[n] < ArrivalTimeTarget for n in N}
    
    
    # Return model
    return Model
