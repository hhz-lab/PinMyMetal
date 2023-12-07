from math import *
import numpy as np
import psycopg2 as pg
import itertools
from itertools import combinations
import os, sys
import getopt


data_dir = "/zinc_prediction/script/"

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid
DBname="zinc"+str(pdbid)

conn = pg.connect("dbname='"+DBname+"' user='sg_display' password='' port='5432'")
cur = conn.cursor()


sql="select pdbid,conc_comma(id_a||'_'||id_b) as id_ab from pre_pre_dist where dist <=2.5 group by pdbid"

cur.execute(sql)
data = cur.fetchall()

print_log = open(data_dir + '/' +pdbid+ '_nb_result' + '/' + 'contact_ab.csv','w')

def bfs(graph, start):
    visited, queue = set(), [start]
    while queue:
        vertex = queue.pop(0)
        if vertex not in visited:
            visited.add(vertex)
            queue.extend(graph[vertex] - visited)
    return visited

def connected_components(G):
    seen = set()
    for v in G:
        if v not in seen:
            c = set(bfs(G, v))
            yield c
            seen.update(c)

def graph(edge_list):
    result = {}
    for source, target in edge_list:
        result.setdefault(source, set()).add(target)
        result.setdefault(target, set()).add(source)
    return result

def concat(l):
    edges = []
    s = list(map(set, l))
    for i, j in combinations(range(len(s)), r=2):
        if s[i].intersection(s[j]):
            edges.append((i, j))
    G = graph(edges)
    result = []
    unassigned = list(range(len(s)))
    for component in connected_components(G):
        union = set().union(*(s[i] for i in component))
        result.append(sorted(union))
        unassigned = [i for i in unassigned if i not in component]
    result.extend(map(sorted, (s[i] for i in unassigned)))
    return result

cid = 0
for i in data:
    pdbid=i[0]
    id_ab=i[1].split(',')
   
    id_list=[]
    for r in id_ab:
        ab_list=[int(r.split('_')[0]),int(r.split('_')[1])]
        id_list.append(ab_list)
    id_list2 = concat(id_list)
    for l in id_list2:
        cid += 1
        for x in l:
            print("%s,%s,%s"%(pdbid,cid,x),file=print_log)
    
print_log.close()

