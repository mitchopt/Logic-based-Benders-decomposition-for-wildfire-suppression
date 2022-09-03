# -*- coding: utf-8 -*-
"""
Implementation of Dijkstra's algorithm using heapq

"""

# Packages
import heapq

# Dijkstra's Algorithm
def ShortestPaths(Nodes, Arcs, InArcs, OutArcs, PlacedRes, Ignitions, Delay):
    
    # Weighted graph encoding
    graph = {u: {v: Arcs[u, v] + (Delay * int(u in PlacedRes))
        for (u, v) in OutArcs[u]} for u in Nodes}
    
    # Distances
    ArrivalTime = {n: float("inf") for n in Nodes}

    # FirePaths
    FirePath = {n: None for n in Nodes}

    # Ignitions
    for n in Ignitions:
        ArrivalTime[n] = 0
        FirePath[n] = []

    # Predecessors
    Pred = {n: None for n in Nodes}

    # Priority queue
    Queue = [(0, n) for n in Ignitions]

    # Dijkstra
    while len(Queue) > 0:
        CurrentDist, CurrentNode = heapq.heappop(Queue)
        if CurrentDist > ArrivalTime[CurrentNode]:
            continue

        for Neigh, Len in graph[CurrentNode].items():
            Dist = CurrentDist + Len
            if Dist < ArrivalTime[Neigh]:
                ArrivalTime[Neigh] = Dist
                FirePath[Neigh] = FirePath[CurrentNode] + [Neigh]
                Pred[Neigh] = CurrentNode
                heapq.heappush(Queue, (Dist, Neigh))

    # Return distances and predecessors
    return ArrivalTime, FirePath, Pred