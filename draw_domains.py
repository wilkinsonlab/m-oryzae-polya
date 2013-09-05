from pysvg.filter import *
from pysvg.gradient import *
from pysvg.linking import *
from pysvg.script import *
from pysvg.shape import *
from pysvg.structure import *
from pysvg.style import *
from pysvg.text import *
from pysvg.builders import *

import sys
import re


z = 1.0
s = svg()
oh = ShapeBuilder()

infile = open(sys.argv[1], 'r')
orderfile = open(sys.argv[2], 'r')

domain_color = {}
colors = ['red', 'yellow', 'blue', 'green',
          'cyan', 'magenta', 'pink', 'navy', 'fucsia']
colors_k = 0
entries = {}
lengths = {}

for line in infile:
    items = line.split('\t')
    species = items[0]
    start = int(items[6])
    end = int(items[7])
    e_value = items[8]
    domain = items[12]
    length = int(items[2])
    if not domain_color.has_key(domain):
        domain_color[domain] = colors[colors_k]
        colors_k += 1
    if not entries.has_key(species):
        entries[species] = []
    entries[species].append((domain, start, end, e_value))
    if not lengths.has_key(species):
        lengths[species] = length

y = 10
for line in orderfile:
    pass
for species in lengths.keys():
    # species = line.strip()
    s.addElement(oh.createLine(0, y + 3, lengths[
                 species], y + 3, strokewidth=1, stroke="black"))
    for (domain, start, end, e_value) in entries[species]:
        s.addElement(oh.createRect(start, y, end - start, 6, 1, 1,
                     strokewidth=1, stroke='black', fill=domain_color[domain]))
    y += 10
s.save('HRP1_domains.svg')

s = svg()
oh = ShapeBuilder()
y = 10
for domain, color in domain_color.items():
    s.addElement(oh.createRect(10, y, 100, 20, 5, 5,
                 strokewidth=1, stroke='black', fill=domain_color[domain]))
    t = text(domain, 120, y + 15)
    s.addElement(t)
    y += 40
s.save('legend.svg')
