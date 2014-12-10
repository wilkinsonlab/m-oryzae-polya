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


z = 1.0
s = svg()
oh = ShapeBuilder()

infile = open(sys.argv[1], 'r')
orderfile = open(sys.argv[2], 'r')
homofile = open(sys.argv[3], 'r')
name = sys.argv[4]
column = int(sys.argv[5])

domain_color = {}
colors = ['red', 'yellow', 'blue', 'green', 'cyan', 'magenta', 'pink', 'navy', 'fucsia', 'gray', 'orange']
colors_k = 0
entries = {}
lengths = {}

for line in infile:
    items = line.strip().split('\t')
    species = items[0]#[0:10]
    start = int(items[6])
    end = int(items[7])
    if len(items) == 11:
        domain = (items[4])
    else:    
        domain = (items[column])
    length = int(items[2])
    if length != 0: 
      if not domain_color.has_key(domain):
         domain_color[domain] = colors[colors_k]
         colors_k += 1
    if not entries.has_key(species):
        entries[species] = []
    entries[species].append((domain, start, end))
    if not lengths.has_key(species):
        lengths[species] = length
        
homologs = {}
for line in homofile:
    sp, ln = line.strip().split("\t")
    homologs[sp] = int(ln)

order = []
for line in orderfile:
    order.append(line.strip())

style=StyleBuilder()
style.setFontSize('3em') #no need for the keywords all the time
t1=text(name,0,30)
t1.set_style(style.getStyle())
s.addElement(t1)
 
x = 170
y = 70
for species in order:
    s.addElement(text(species, 0, y))  
    if species in homologs:
        s.addElement(oh.createLine(x, y , x + homologs[species], y, strokewidth=1, stroke="black"))
    if lengths.has_key(species) and lengths[species] != 0:
        for (domain, start, end) in sorted(entries[species], key=lambda x: x[1]):
            s.addElement(oh.createRect(x + start, y-3, end - start, 6, 1, 1,
                         strokewidth=1, stroke='black', fill=domain_color[domain]))
    y += 20
    
for domain, color in domain_color.items():
    s.addElement(oh.createRect(10, y, 100, 20, 5, 5,
                 strokewidth=1, stroke='black', fill=domain_color[domain]))
    t = text(domain, 120, y + 12)
    s.addElement(t)
    y += 40
    
s.save(sys.argv[1] + '.domains.svg')