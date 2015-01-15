from pysvg.filter import *
from pysvg.gradient import *
from pysvg.linking import *
from pysvg.script import *
from pysvg.shape import *
from pysvg.structure import *
from pysvg.style import *
from pysvg.text import *
from pysvg.builders import *
import pysvg
import sys
import re




infile = open(sys.argv[1], 'r')

colors = ['red', 'yellow', 'blue', 'cyan', 'magenta', 'pink', 'navy', 'olivedrab', 'gray', 'orange', 'black', 'green']
entry = {}
labels = {}
complex = ""
for line in infile:
    items = line.strip().split('\t')
    if items[0][0] == ">":
         complex = items[0].replace(">", "")   
         entry[complex] = []  
         labels[complex] = [] 
    else:
         entry[complex].append(items[1:])
         labels[complex].append(items[0])

x = 10
y = 10
for complex, lst in entry.items():
    s = svg()
    oh = ShapeBuilder()
    style=StyleBuilder()
    style.setFontSize('1em') #no need for the keywords all the time
    t1=text(complex + " " ,x-10,y)
    t1.set_style(style.getStyle())
    s.addElement(t1)
    y += 10  
    x = 20
    for species in zip(*lst):
        k = 0
        for protein in species:
            if protein == "1":
                 s.addElement(oh.createCircle(x-10, y, 2,  fill=colors[k]))
            k += 1     
            x += 10
        x = 20          
        y += 10   
    for label in labels[complex]:
        style=StyleBuilder()
        style.setFontSize('0.5em') #no need for the keywords all the time
        t1=text(label,x-10,y)
        t1.set_style(style.getStyle())
        t1.set_transform("rotate(45 " + str(x-10) + " " + str(y) + ")")
        s.addElement(t1)    
        x += 10
    t1=text(" ", 10 ,y+30)    
    s.addElement(t1)    
    s.save('proteins_' + complex.replace("/", "") + '.svg')
    x = 10    
    y = 10     
