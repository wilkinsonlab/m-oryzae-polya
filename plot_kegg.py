# coding: utf-8

import numpy as np, sys, re, os.path, math
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
    term, val, pvalue, desc, type = line.strip().split("\t")
    if len(desc) > 60:
        desc = desc[0:60] + "..."        
    xlabels.append(desc)
    ylabels.append(1-math.log((float(pvalue))))
    types.append(type)


fig = plt.figure()
width = 0
ind = np.arange(len(ylabels))
barlist=plt.barh(ind, ylabels)
plt.yticks(ind+0.5, xlabels, ha='right',rotation=0, fontsize=14)
for i, type in enumerate(types):
    if type == "UP":
        barlist[i].set_color('r')
    elif type == "DOWN":
        barlist[i].set_color('b')    

        
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(18.5,10.5)
fd = {'fontsize':24}
title = re.sub("__", " ", os.path.basename((filename)))
title = re.sub("_plot", "", title)
title = re.sub("_vs_", u" -> ", title)
title = re.sub("2D4", u"Δrbp35", title)
title = re.sub("--N", u"-MM-N", title)
title = re.sub("--C", u"-MM-C", title)
plt.title(title, y=1.08,fontdict=fd)

l1 = Line2D([], [], linewidth=3, color="r") 
l2 = Line2D([], [], linewidth=3, color="b") 
l3 = Line2D([], [], linewidth=3, color="b") 
#legend([l1, l2], ["Up-regulated", "Down-regulated"], fontsize=20,loc='upper center',bbox_to_anchor=(0.5, 1.5)) 
plt.savefig(filename.replace("txt", "pdf"))
