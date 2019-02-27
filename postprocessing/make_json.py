''' MAKE A JSON FILE SUITABLE FOR READING BY D3 NETWORK BROWSER'''

#========================
#=======  JSON  =========
#========================

import os
import json
import networkx as nx
from networkx.readwrite import json_graph
#-----------------------------
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-wd', '--workingdirectory',    help='full coexpression results directory')
args = parser.parse_args()
#-----------------------------
def write_json(edgelist, outputfile, nodenames={}, theseclasses={}):
	G=nx.Graph()
	for i, edge in enumerate(edgelist):
		G.add_edge(edge[0],edge[1],weight=float(edge[2])) # specify edge data
	#put in a blank edge-free node for each one that is missing
	allnodes = list(G.nodes(data=False))
	if len(allnodes)>0:
		for missing in set(range(max(allnodes)+1))-set(allnodes):
			G.add_node(missing)
		for n in G:
			G.node[n]['nodesize'] = 5000
			try:
				G.node[n]['name'] = nodenames[n]
			except:
				G.node[n]['name'] = 'none'
			try:
				G.node[n]['group'] = theseclasses[nodenames[n]]
			except:
				G.node[n]['group'] = 'none'
	d = json_graph.node_link_data(G)
	jo = open(outputfile,'w')
	jo.write('<? header("Access-Control-Allow-Origin: https://baillielab.net") ?>\n')
	json.dump(d, jo)
	jo.close()
#-----------------------------
indscoresfile = os.path.join(args.workingdirectory, "coex.individual_scores")
f=open(indscoresfile)
lines = [x.strip().split('\t') for x in f.readlines()]
f.close()
sig = {x[0].split(" ")[0]:(float(x[-1])<=0.05)*3 for x in lines[1:]}
indices = {x[0].split(" ")[0]:i for i,x in enumerate(lines[1:])}
namesofnodes = {i:x[0].split(" ")[0] for i,x in enumerate(lines[1:])}
promids = {x[2]:x[0].split(" ")[0] for x in lines[1:]}

networkfile = os.path.join(args.workingdirectory, "coex.network.php")
f=open(networkfile)
lines = [x.strip().split('\t') for x in f.readlines()]
f.close()
edges=[]
for line in lines:
	try:
		promids[line[0]]
	except:
		continue
	try:
		promids[line[1]]
	except:
		continue
	if float(line[2]) > 0:
		edges.append([indices[promids[line[0]]],indices[promids[line[1]]],float(line[2])])
jsonout = os.path.join(args.workingdirectory, "hitsnetwork.php")
write_json(edges, jsonout, nodenames=namesofnodes, theseclasses=sig)






