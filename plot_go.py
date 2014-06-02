import numpy as np, sys, re, os.path
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.pyplot import *
filename = sys.argv[1]
rcParams.update({'figure.autolayout': True})
xlabels = []
ylabels = []
types= []
fich = open(filename, "r")
for line in fich:
    if line == "": continue
    desc, val, type = line.strip().split("\t")
    if len(desc) > 40:
        desc = desc[0:39] + "..."        
    xlabels.append(desc)
    ylabels.append(int(val))
    types.append(type)


fig = plt.figure()
width = 0
ind = np.arange(len(ylabels))
barlist=plt.bar(ind, ylabels)
plt.xticks(ind+0.5, xlabels, ha='right',rotation=70, fontsize=14)
for i, type in enumerate(types):
    if type == "sgl":
        barlist[i].set_color('orange')
    elif type == "apa":
        barlist[i].set_color('violet')    
    elif type == "MF":
        barlist[i].set_color('b')
        
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5,10.5)
fd = {'fontsize':24}
title = re.sub("_go_plot_", " ", os.path.basename((filename)))
title = re.sub("\.txt", "", title)
#title = "Plant-only expressed genes GO ontologies"
plt.title(title, y=1.08,fontdict=fd)
l1 = Line2D([], [], linewidth=3, color="orange") 
l2 = Line2D([], [], linewidth=3, color="violet") 
l3 = Line2D([], [], linewidth=3, color="b") 
legend([l1, l2], ["Single cut", "APA"], fontsize=20) 
plt.savefig(filename.replace("txt", "pdf"))
