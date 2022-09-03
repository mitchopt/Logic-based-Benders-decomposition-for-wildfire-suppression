# -*- coding: utf-8 -*-
"""
Main LBBD model

"""

# Packages
from shortest_paths import ShortestPaths
import gurobipy as gp
import math


# ------------------------ #
# --- LBBD Formulation --- #
# ------------------------ #
def FireLBBD(Nodes, Arcs, ResAtTime, Ignitions, Delay, ArrivalTimeTarget, 
             TimeLimit, GurobiSeed, ZStart, Heuristic):

    # Epsilon
    EPS = 0.0001
    
    # Gurobi model
    Model = gp.Model()
    Model.setParam("Threads", 1)
    # Model.setParam("OutputFlag", 0)
    Model.setParam("LazyConstraints", 1)
    
    # Statistics
    Model._ShortestPathProblemsSolved = 0
    Model._OptimalityCuts = 0
    Model._FeasibilityCuts = 0
    
    
    # Seed if one is given
    if GurobiSeed is not None:
        Model.setParam("Seed", GurobiSeed)
            
        
    # Set time limit if given
    if TimeLimit is not None:
        Model.setParam("TimeLimit", TimeLimit)
    
    
    # Generate in- and outarcs
    InArcs = {n: [] for n in Nodes}
    OutArcs = {n: [] for n in Nodes}
    for a in Arcs:
        OutArcs[a[0]].append((a[0], a[1]))
        InArcs[a[1]].append((a[0], a[1]))
    

    # Decision variables
    # IsResource[n, t] = 1 if we put 
    # a resource on node n at time t
    IsResource = {(n, t): Model.addVar(
        vtype=gp.GRB.BINARY) for n in Nodes for t in ResAtTime}


    # Benders variables
    # DoesBurn[n] = 1 if node n burns before the target time
    DoesBurn = {n: Model.addVar(vtype=gp.GRB.BINARY) for n in Nodes}


    # Constraints
    # Up to Res[t] resources at time t
    ResPerTime = {t: Model.addConstr(gp.quicksum(
        IsResource[n, t] for n in Nodes) <= ResAtTime[t]) for t in ResAtTime}

    # Up to one resource per node
    ResPerNode = {n: Model.addConstr(gp.quicksum(
        IsResource[n, t] for t in ResAtTime) <= 1) for n in Nodes}

    # No resources at ignition node
    for (n, t) in IsResource:
        if n in Ignitions:
            IsResource[n, t].ub = 0

    #  The objective is the number of burned nodes
    Model.setObjective(gp.quicksum(DoesBurn[n] for n in DoesBurn),
                       gp.GRB.MINIMIZE)
    

    # Ass initial cuts
    Model._ShortestPathProblemsSolved += 1
    ArrNone, PathNone, Pred = ShortestPaths(
        Nodes, Arcs, InArcs, OutArcs, set(), Ignitions, Delay)

    # Cut on each node
    for n in DoesBurn:

        # Cut on burned nodes
        if ArrNone[n] < ArrivalTimeTarget:
    
            # Minimum interdictions needed
            Gap = ArrivalTimeTarget - ArrNone[n]
            Required = math.ceil(Gap / Delay)
                
            # Fire path expression
            FirePathExpr = gp.quicksum(
                IsResource[nn, t] / Required for nn in PathNone[n][:-1]
                for t in ResAtTime if t <= ArrNone[nn] + (Required - 1) * Delay)
    
            # Add initial cut
            Model.addConstr(DoesBurn[n] >= 1 - FirePathExpr)  
    

    # Callback
    def Callback(model, where):
        if where == gp.GRB.Callback.MIPSOL:
        
            # Retrieve incumbent solution
            IsResourceGet = model.cbGetSolution(IsResource)
            DoesBurnGet = model.cbGetSolution(DoesBurn)
            IsResourceV = {n: round(IsResourceGet[n]) for n in IsResourceGet}
            DoesBurnV = {n: round(DoesBurnGet[n]) for n in DoesBurnGet}

            
            Incumbent = set(  # Set of nodes with a resource
                n for n in Nodes for t in ResAtTime if IsResourceV[n, t] > .5)


            # Solve shortest paths problem
            Model._ShortestPathProblemsSolved += 1
            ArrivalTime, FirePath, Pred = ShortestPaths(
                Nodes, Arcs, InArcs, OutArcs, Incumbent, Ignitions, Delay)


            # Cut on nodes
            for n in DoesBurn:

                # FEASIBILITY
                if any(IsResourceV[n, t] > .5 for t in ResAtTime):
                    tt = min(t for t in ResAtTime if IsResourceV[n, t] > .5)

                    if ArrivalTime[n] < tt:
                        
                        # Minimum interdictions needed
                        Gap = tt - ArrNone[n]
                        Required = math.ceil(Gap / Delay)
    
                        # Fire path expression
                        FirePathExpr = gp.quicksum(
                            IsResource[nn, t] for nn in FirePath[n][:-1] 
                            for t in ResAtTime if t <= ArrivalTime[nn] +
                            (Required - 1) * Delay if IsResourceV[nn, t] < .5)
    
                        # Feasibility cut
                        Model.cbLazy(1 - IsResource[n, tt] + FirePathExpr / Required  >= 1)
                        Model._FeasibilityCuts += 1
    
    
                # OPTIMALITY
                if DoesBurnV[n] < 1 - EPS:
                    if ArrivalTime[n] < ArrivalTimeTarget - EPS:
                    
                        # Find the minimum interdictions required
                        Gap = ArrivalTimeTarget - ArrivalTime[n]
                        Required = Gap // Delay + 1*int(Gap % Delay != 0)
                        
                        # Fire path expression
                        FirePathExpr = gp.quicksum(
                            IsResource[nn, t] for nn in FirePath[n][:-1]
                            for t in ResAtTime if t <= ArrivalTime[nn] +
                            (Required - 1) * Delay if IsResourceV[nn, t] < .5)
    
                        # Optimality cut
                        Model.cbLazy(DoesBurn[n] >= 1 - FirePathExpr / Required)
                        Model._OptimalityCuts += 1
    

    # Starting solution
    if ZStart is not None:
        for _ in IsResource:
            IsResource[_].Start = 0
            if _ in ZStart:
                IsResource[_].Start = 1
        
        # Evaluate starting solution
        StartRes = set(n for (n, t) in ZStart)
        ArrStart, _, __ = ShortestPaths(
            Nodes, Arcs, InArcs, OutArcs, StartRes, Ignitions, Delay)
        for n in DoesBurn:
            if ArrStart[n] < ArrivalTimeTarget:
                DoesBurn[n].Start = 1
            else:
                DoesBurn[n].Start = 0
    
    
    

    # Solve heuristically
    if Heuristic:
        Model.setParam("OutputFlag", 0)

        i = 0
        T = list(ResAtTime.keys())
        while i < len(ResAtTime):
            t = T[i]
            
            # Fix right hand sides
            for n in Nodes:
                for tt in ResAtTime:
                    if tt > t:
                        ResPerTime[tt].RHS = 0
            
            # Solve iteration t
            Model.optimize(Callback)
            
            # Fix solution for t
            for n in Nodes:
                if IsResource[n, t].x > 0.5:
                    Model.addConstr(IsResource[n, t] == 1)
            
            # Reset right hand sides
            for n in Nodes:
                for tt in ResAtTime:
                    if tt > t:
                        ResPerTime[tt].RHS = ResAtTime[tt]
            
            # Iterate
            i += 1
            
    else:

        # Otherwise solve exactly
        Model.optimize(Callback)


    # Save results
    Model._Optimal = set(
        n for n in Nodes for t in ResAtTime if IsResource[n, t].x > .1)
    Model._ArrivalTime, Model._FirePath, Model._Pred = ShortestPaths(
        Nodes, Arcs, InArcs, OutArcs, Model._Optimal, Ignitions, Delay)
    Model._Burned = {n: Model._ArrivalTime[n] < ArrivalTimeTarget for n in Nodes}
    Model._Theta = {n: DoesBurn[n].x for n in DoesBurn}


    # Return model
    return Model, [(n, t) for (n, t) in IsResource if IsResource[n, t].x > 0.5]
