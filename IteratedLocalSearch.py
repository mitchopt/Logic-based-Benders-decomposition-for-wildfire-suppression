# -*- coding: utf-8 -*-
"""
Our implementation of the Iterated Local Search (ILS) metaheuristic from
the related paper,

    Iterated Local Search for the Placement of Wildfire Suppression Resources
    Mendes and Alvelos, 2022, European Journal of Operations Research


"""

# Packages
from shortest_paths import ShortestPaths
import numpy as np


# ------------------------------------------- #
# --- Iterated Local Search Metaheuristic --- #
# ------------------------------------------- #
def FireILS(Nodes, Arcs, ResAtTime, Ignitions, Delay, ArrivalTimeTarget,
            MultiStarts, p1, p2, MaxNeighbours, MaxModifications, 
            MaxFailures, MaxNoImprovements, MaxCandidates):
    EPS = 0.0001

    # Generate in- and outarcs
    InArcs = {n: [] for n in Nodes}
    OutArcs = {n: [] for n in Nodes}
    for a in Arcs:
        OutArcs[a[0]].append((a[0], a[1]))
        InArcs[a[1]].append((a[0], a[1]))
    
    
    #
    # ----------------------------------- #
    # ----- Construct random solution --- #
    # ----------------------------------- #
    def ConstructRandomSolution(MaxCandidates):
        
        # Initial zero solution
        ZSol = {(n, t): 0 for n in Nodes for t in ResAtTime}
        
        # Add resources until all added
        while sum(ZSol.values()) < sum(ResAtTime.values()):
            
            # Get the earliest release time  of an available resource
            tt = min(t for (n, t) in ZSol if sum(
                ZSol[n, t] for n in Nodes) < ResAtTime[t])
            
            # Get nodes with res
            Incumbent = set(n for n in Nodes for t in ResAtTime if ZSol[n, t] > .5)
            
            # Run Dijkstra to get arrival times
            ArrivalTime, _, __ = ShortestPaths(
                Nodes, Arcs, InArcs, OutArcs, Incumbent, Ignitions, Delay)
            
            # Identify unburned nodes at time tt in order of arrival time
            Unburned = sorted([(ArrivalTime[n], n) for n in Nodes
                               if ArrivalTime[n] >= tt and not
                               any(ZSol[n, t] > 0.5 for t in ResAtTime)])
            
            
            # Get candidate nodes
            Candidates = [n for (a, n) in Unburned][:MaxCandidates]
            # Get the nodes that will burn first if there is no change
            # minArrival = min(_[0] for _ in Unburned)
            # Candidates = [n for (a, n) in Unburned if a < minArrival + EPS]
            
            # Chose a random candidate node
            Choice = Candidates[np.random.randint(len(Candidates))]
            
            # Add a resource to the chosen node at the chosen
            ZSol[Choice, tt] = 1
        
        return ZSol
    
    
    #
    # ---------------------------------------------------------- #
    # ----- Generate many random solutions and keep the best --- #
    # ---------------------------------------------------------- #
    def MultiStartConstructiveHeuristic(Iterations, MaxCandidates):
        
        # Starting best initial solution
        ZSolBest = ConstructRandomSolution(MaxCandidates)
        
        # Evaluate initial solution
        ArrivalTime, _, __ = ShortestPaths(
            Nodes, Arcs, InArcs, OutArcs, set(
                n for n in Nodes for t in ResAtTime 
                if ZSolBest[n, t] > .5), Ignitions, Delay)
        
        # Starting best objective value
        objBest = len([n for n in Nodes if ArrivalTime[n] < ArrivalTimeTarget])

        for _ in range(Iterations - 1):
            
            # Construct a random solution
            ZSol = ConstructRandomSolution(MaxCandidates)
            
            # Evaluate random solution
            ArrivalTime, _, __ = ShortestPaths(
                Nodes, Arcs, InArcs, OutArcs, set(
                    n for n in Nodes for t in ResAtTime 
                    if ZSol[n, t] > .5), Ignitions, Delay)
        
            # Calculate objective value
            objNew = len([n for n in Nodes if ArrivalTime[n] < ArrivalTimeTarget])
            
            # Check if better
            if objNew < objBest:
                ZSolBest = ZSol
                objBest = objNew
        
        # Return best
        return ZSolBest
    


    # --- Neighbourhood of a single node --- #
    def NeighOneNode(n):
        (i, j) = n
        return {(i - 1, j), (i + 1, j), (i, j - 1), 
                (i, j + 1), (i - 1, j - 1), (i + 1, j - 1), 
                (i - 1, j + 1), (i + 1, j + 1)}.intersection(set(Nodes))
    
    

    # -------------------- #
    # --- Local Search --- #
    # -------------------- #
    def LocalSearch(ZSol):
        Improvement = 1
        
        
        # Temporary copy
        TempZSol = ZSol.copy()
            
        # Go until no improvement
        while Improvement > 0.5:
            Improvement = 0
        
            # The current set of nodes with a resource
            HasRes = set(n for n in Nodes for t in ResAtTime if ZSol[n, t] > .5)
        
            # Current arrival times and objective values
            ArrivalTime, _, __ = ShortestPaths(Nodes, Arcs, InArcs, OutArcs, HasRes, Ignitions, Delay)
            
            # Current best objective value
            BestObj = len([n for n in Nodes if ArrivalTime[n] < ArrivalTimeTarget])
            
            # For each node with a resource we try moving that resource
            # to another node. Here n is the node we removed a resource from
            for n in list(HasRes):
                
                # Get the time period of the resource we are removing
                tt = min(t for t in ResAtTime if ZSol[n, t] > 0.5)
                
                # Set of nodes with resources after we remove n
                RemovedHasRes = HasRes - {n}
                
                # Get arrival times after removal of n
                RemovedArrivalTime, _, __ = ShortestPaths(
                    Nodes, Arcs, InArcs, OutArcs, RemovedHasRes, Ignitions, Delay)
                
                # Create extended neighbourhood of 
                # the nodes that still have resources
                Neighbourhood = set()
                for i in list(RemovedHasRes):
                    Neighbourhood |= NeighOneNode(i)
                Neighbourhood -= RemovedHasRes
                
                # Get the sorted neighbours that are not burned yet at time tt
                SortedUnburned = sorted([(RemovedArrivalTime[i], i) for i in Neighbourhood
                                   if RemovedArrivalTime[i] >= tt])
                
                # Keep only the next nodes to burn in the neihbourhood
                Neighbourhood = [i for (_, i) in SortedUnburned]
                Neighbourhood = Neighbourhood[:MaxNeighbours]
                
                # Iterate the extended neighbours
                # nn is the node we try adding to
                for nn in Neighbourhood:
                    
                    # Nodes with resource after adding one to nn
                    AddHasRes = RemovedHasRes | {nn}
                    
                    # Get arrival times after adding the resource
                    AddArrivalTime, _, __ = ShortestPaths(
                        Nodes, Arcs, InArcs, OutArcs, AddHasRes, Ignitions, Delay)
                    
                    # Evaluate feasibility of the move
                    INFEASIBLE = False
                    for (nnn, t) in ZSol:
                        if nnn != nn:  # nn is fine by construction
                            if ZSol[nnn, t] > 0.5:
                                if AddArrivalTime[nnn] < t:
                                    INFEASIBLE = True
                                    break
                    
                    # Skip nn if the move is infeasible
                    if INFEASIBLE:
                        continue
                    
                    # Otherwise check the objective value
                    ObjAfterMove = len([i for i in Nodes if AddArrivalTime[i] < ArrivalTimeTarget])
                    
                    # If the objective is better then update the solution
                    NextImprovement = BestObj - ObjAfterMove
                    if NextImprovement > Improvement:
                        BestObj = ObjAfterMove
                        Improvement = NextImprovement
                        TempZSol = ZSol.copy()
                        TempZSol[n, tt] = 0
                        TempZSol[nn, tt] = 1
                
            # If we improved then update ZSol for the next iteration
            if Improvement > 0:
                ZSol = TempZSol
                
        # Return the output
        return ZSol

    
    # --------------------
    # --- Pertubations ---
    # --------------------
    #
    # 1. Remove a resource from the node 
    # with the largest deployment time
    def Pertubation1(ZSol):
        ZSolNew = ZSol.copy()
        
        tMax = max(_[1] for _ in ZSolNew if ZSolNew[_] > 0.5)
        tNodes = [n for (n, t) in ZSolNew if ZSolNew[n, t] > 0.5 and t == tMax]
        nChoice = tNodes[np.random.randint(len(tNodes))]
        ZSolNew[nChoice, tMax] = 0
        return ZSolNew
    
    
    # 2. Add a resource to a node
    #  Basically identical to ConstructRandomSolution
    #  except that we only add one resource
    def Pertubation2(ZSol):
        ZSolNew = ZSol.copy()
        
        # If there are no resources available, do nothing
        if sum(ZSol.values()) >= sum(ResAtTime.values()) - EPS:
            return ZSol
        
        # Get the earliest release time  of an available resource
        tt = min(t for (n, t) in ZSolNew if sum(
            ZSolNew[n, t] for n in Nodes) < ResAtTime[t])
        
        # Get nodes with resources
        Incumbent = set(n for n in Nodes for t in ResAtTime if ZSolNew[n, t] > .5)
        
        # Run Dijkstra to get arrival times
        ArrivalTime, _, __ = ShortestPaths(
            Nodes, Arcs, InArcs, OutArcs, Incumbent, Ignitions, Delay)
        
        # Identify candidate nodes for a new resource
        Unburned = sorted([(ArrivalTime[n], n) for n in Nodes
                           if ArrivalTime[n] >= tt and not
                           any(ZSolNew[n, t] > 0.5 for t in ResAtTime)])

        # Get the nodes that will burn first if there is no change
        minArrival = min(_[0] for _ in Unburned)
        Candidates = [n for (a, n) in Unburned if a < minArrival + EPS]
        
        # Chose a random candidate node
        Choice = Candidates[np.random.randint(len(Candidates))]
        
        # Add a resource to the chosen node at the chosen
        ZSolNew[Choice, tt] = 1
        return ZSolNew

    # Make random modifications until reaching either the maximum number
    # of modifications or the maximum number of failures
    def Pertubation3(ZSol, MaxModifications, MaxFailures):
        Mod = 0
        Fail = 0
        ZSolNew = ZSol.copy()
        while Mod < MaxModifications and Fail < MaxFailures:
            
            
            
            # Get nodes with resources
            HasResource = set(n for n in Nodes for t in ResAtTime if ZSolNew[n, t] > .5)
            
            # print(HasResource)
            n = list(HasResource)[np.random.randint(len(HasResource))]
            tt = min(t for (nn, t) in ZSolNew if nn == n and ZSolNew[n, t] > 0.5)

            # Set of nodes with resources after we remove n
            RemovedHasRes = HasResource - {n}
            
            # Get arrival times after removal of n
            RemovedArrivalTime, _, __ = ShortestPaths(
                Nodes, Arcs, InArcs, OutArcs, RemovedHasRes, Ignitions, Delay)
            
            # Create a broader neighbourhood of candidate nodes
            Neighbourhood = [nn for nn in Nodes if nn != n and ZSolNew[nn, tt] < 0.5]
            
            # Get the sorted neighbours that are not burned yet at time tt
            SortedUnburned = sorted([(RemovedArrivalTime[i], i) for i in Neighbourhood
                               if RemovedArrivalTime[i] >= tt])
            
            # Get the "MaxNeighbours" next nodes that burn if no change
            Candidates = [n for (a, n) in SortedUnburned]
            Candidates = Candidates[:MaxNeighbours]
            
            # Chose a random candidate node
            Choice = Candidates[np.random.randint(len(Candidates))]
            
            # Add the resource to HasRes
            AddHasRes = RemovedHasRes | {Choice}
            
            # Get arrival times after adding the resource
            AddArrivalTime, _, __ = ShortestPaths(
                Nodes, Arcs, InArcs, OutArcs, AddHasRes, Ignitions, Delay)
            
            # Evaluate feasibility of the move
            INFEASIBLE = False
            for (nn, t) in ZSolNew:
                if nn != Choice:  # Choice is fine by construction
                    if ZSolNew[nn, t] > 0.5:
                        if AddArrivalTime[nn] < t:
                            INFEASIBLE = True
                            break
            
            # Step
            if INFEASIBLE:
                Fail += 1
            else:
                Fail = 0
                Mod += 1
                ZSolNew[n, tt] = 0
                ZSolNew[Choice, tt] = 1
                
        
        return ZSolNew
    
    # ------------------------------------------------------- #
    # --- Choose a pertubation according to probabilities --- #
    # ------------------------------------------------------- #
    def Perturbate(ZSol, p1, p2, MaxModifications, MaxFailures):
        
        # p1 is 0 if there are no resources left
        if sum(ZSol.values()) < 0.5:
            p1 = 0
        
        # p2 is 0 if there are no available resources
        if sum(ZSol.values()) < sum(ResAtTime.values()):
            p2 = 0
            
        # Chose a pertubation
        Rand = np.random.random()
        if Rand < p1:
            return Pertubation1(ZSol)
        else:
            if Rand < p1 + p2:
                return Pertubation2(ZSol)
            else:
                return Pertubation3(ZSol, MaxModifications, MaxFailures)
            
            
            
    # ------------------------------------------------
    # --- Complete Iterated Local Search Heuristic ---
    # ------------------------------------------------
    def IteratedLocalSearch(MultiStarts, p1, p2, MaxModifications, 
                            MaxFailures, MaxNoImprovements, MaxCandidates):
        NoImprovements = 0
        
        # Get starting solution
        print("Generating multistart solution")
        ZSolBest = MultiStartConstructiveHeuristic(MultiStarts, MaxCandidates)
        
        # Evaluate starting solution
        HasRes = set(n for (n, t) in ZSolBest if ZSolBest[n, t] > 0.5)
        ArrivalTime, _, __ = ShortestPaths(
            Nodes, Arcs, InArcs, OutArcs, HasRes, Ignitions, Delay)
        ObjVal = len([n for n in Nodes if ArrivalTime[n] < ArrivalTimeTarget])
        
        
        
        # Initial local search
        print("Initial local search")
        ZSolBest = LocalSearch(ZSolBest)

        # Start ILS
        Stop = False
        while not Stop:
            print("Best:", ObjVal)
            
            # Perturb solution
            ZSolTemp = Perturbate(ZSolBest, p1, p2, MaxModifications, MaxFailures)
            
            # Local search from perturbed solution
            ZSolTemp = LocalSearch(ZSolTemp)        
            
            # Evaluate new solution
            HasRes = set(n for (n, t) in ZSolTemp if ZSolTemp[n, t] > 0.5)
            ArrivalTime, _, __ = ShortestPaths(
                Nodes, Arcs, InArcs, OutArcs, HasRes, Ignitions, Delay)
            NewObjVal = len([n for n in Nodes if ArrivalTime[n] < ArrivalTimeTarget])
            
            
            # Check if better
            if NewObjVal < ObjVal:
                
                # Update best
                NoImprovements = 0
                ZSolBest = {(n, t): ZSolTemp[n, t] for (n, t) in ZSolTemp}
                ObjVal = NewObjVal
                
            # Otherwise iterate the
            else:  # stopping criteria
                NoImprovements += 1
                
            # Check the stopping criterion
            if NoImprovements >= MaxNoImprovements:
                Stop = True

        return ZSolBest, ObjVal
            
            
    return IteratedLocalSearch(
        MultiStarts, p1, p2, MaxModifications, MaxFailures,
        MaxNoImprovements, MaxCandidates)


