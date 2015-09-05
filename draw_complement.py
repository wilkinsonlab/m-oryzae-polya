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

colors = ['red', 'yellow', 'blue', 'cyan', 'magenta', 'pink', 'navy', 'olive', 'gray', 'orange', 'black', 'green', 'aliceblue', 'burlywood', 'darkgray', 'darkslategray', 'darkorange', 'darkkhaki', 'darkslateblue', 'deeppink', 'dimgray', 'gainsboro', 'dodgerblue', 'darkolivegreen', 'darkblue', 'yellow', 'blue', 'cyan', 'magenta', 'pink', 'navy', 'olive', 'gray', 'orange', 'black', 'green', 'aliceblue', 'burlywood', 'darkgray', 'darkslategray', 'yellow', 'blue', 'cyan', 'magenta', 'pink', 'navy', 'olive', 'gray', 'orange', 'black', 'green', 'aliceblue', 'burlywood', 'darkgray', 'darkslategray', 'yellow', 'blue', 'cyan', 'magenta', 'pink', 'navy', 'olive', 'gray', 'orange', 'black', 'green', 'aliceblue', 'burlywood', 'darkgray', 'darkslategray', 'yellow', 'blue', 'cyan', 'magenta', 'pink', 'navy', 'olive', 'gray', 'orange', 'black', 'green', 'aliceblue', 'burlywood', 'darkgray', 'darkslategray']
colors = {
0: 'black',
1: 'blue',
2: 'green',
3: 'yellow',
4: 'orange',
5: 'magenta',
6: 'pink',
7: 'red',
8: 'red'

}
entries = []
proteins = []

proteins = infile.readline().strip().split('\t')
for line in infile.readlines():
	items = line.strip().split('\t')
	entries.append(items[1:])

x = 10
y = 10
s = svg()
oh = ShapeBuilder()
style=StyleBuilder()
style.setFontSize('1em') #no need for the keywords all the time
t1=text(sys.argv[1] + " " ,x-10,y)
t1.set_style(style.getStyle())
s.addElement(t1)
y += 10  
x = 20
for entry in entries:
	k = 0
	for e in entry:
		#if e == "1":
			 #s.addElement(oh.createCircle(x-10, y, 2,  fill=colors[k]))
		#	 s.addElement(oh.createCircle(x-10, y, 2,  fill=('grey', 'green')[k%2]))
		t1=text(e,x-10,y)
		style.setFilling(colors[int(e)])
		t1.set_style(style.getStyle())				
		s.addElement(t1)    
		k += 1     
		x += 10
	x = 20          
	y += 10   
for protein in proteins:
	style=StyleBuilder()
	style.setFontSize('0.5em') #no need for the keywords all the time
	t1=text(protein,x-10,y)
	t1.set_style(style.getStyle())
	t1.set_transform("rotate(45 " + str(x-10) + " " + str(y) + ")")
	s.addElement(t1)    
	x += 10
t1=text(" ", 10 ,y+30)    
s.addElement(t1)    
s.save('proteins_' + sys.argv[1].replace("/", "") + '.svg')
x = 10    
y = 10     
