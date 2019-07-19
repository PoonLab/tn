import math

if False:
    # settings for Seattle
    cutoff = 0.0144
    newyear = 2012
    code = 'stD'
    nodecolor = 'dodgerblue'
    boldcolor = '#0051a5'
    origin = 1999
    scale_factor = 20.
    padding = 0.2
else:
    # settings for Northern Alberta
    cutoff = 0.0112
    newyear = 2013
    code = 'naD'
    nodecolor = 'orange2'
    boldcolor = '#b47300'
    origin = 2006
    scale_factor = 10.
    padding = 0.5


handle = open('/home/art/work/maup/data/pub1/{}.txt'.format(code))
header = next(handle)

# pass through once to calculate degree size
dotfile = open('/home/art/papers/maup/images/{}.dot'.format(code), 'w')
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

print (count)
print (newcount)

