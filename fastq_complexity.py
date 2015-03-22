#!/usr/bin/python

import sys

readDic= {}

Nbre_reads = 0

Nbre_lines = 0

F = open(sys.argv[1])

for line in F:

 Nbre_lines += 1

 if Nbre_lines % 4 == 2:

  Nbre_reads += 1

  readDic[line] = readDic.get(line, 0) + 1

F.close()

#print "%s reads" % Nbre_reads

#print "%s distinct sequences" % (len(readDic))

print sys.argv[1], "%f complexity" % (len(readDic)/float(Nbre_reads))
