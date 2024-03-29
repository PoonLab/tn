import math

# settings for Seattle Example
cutoff = 0.0160
newyear = 2012
code = 'sttn93'
nodecolor = 'dodgerblue'
boldcolor = '#0051a5'
origin = 1999
scale_factor = 20.
padding = .9

handle = open('{}.txt'.format(code))
header = next(handle)

# pass through once to calculate degree size
dotfile = open('{}.dot'.format(code), 'w')
dotfile.write("""
graph seattle
{
    node [label="" shape="circle" style="filled" color="white"];
    edge [color="grey" penwidth=5];
    graph [outputorder="edgesfirst"];
""")

nodes = {}
for line in handle:
    id1, id2, dist = line.strip().split(',')
    dist = float(dist)
    
    for node in [id1, id2]:
        if node not in nodes:
            year = int(node.split('_')[-1])
            nodes.update({node: {'year': year, 'degree': 0}})
        
    if dist < cutoff:
        # draw the edge
        nodes[id1]['degree'] += 1
        nodes[id2]['degree'] += 1
        dotfile.write("\t{}--{} [len={}];\n".format(
            id1, id2, 100*dist + padding
        ))

# write out the nodes
count = 0
newcount = 0

for node in nodes:
    if nodes[node]['degree'] == 0:
        # do not display singletons
        count += 1
        if nodes[node]['year'] == newyear:
            newcount += 1
        continue
    dotfile.write("\t{} [fillcolor=\"{}\" width={}];\n".format(
        node, 
        boldcolor if nodes[node]['year']==newyear else nodecolor,
        (nodes[node]['year']-origin) / scale_factor
    ))    

dotfile.write("}\n")
dotfile.close()

quit()

print (count)
print (newcount)

