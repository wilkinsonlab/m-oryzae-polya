
import sys
import re

fp = sys.stdin
seq = sys.argv[1]
score = ""
s = []
flag = False
for line in fp:
    if line[0] == "\n":
        if not flag:
            score = ""
            s = []
            continue
        print score
        for x in s:
            print x
        print
        score = ""
        s = []
        flag = False
        continue
    line = line.strip()

    if line[0] == "#":
        continue

    if line[0] == "a":
        score = line
    if line[0] == "s":
        if re.match("s " + seq, line):
            s.insert(0, line)
            flag = True
        else:
            s.append(line)
