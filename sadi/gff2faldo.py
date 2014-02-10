


print """

# baseURI: http://example.org/insdc
# imports: file:/Users/bolleman/git/INSDC/insdc.ttl
# imports: http://biohackathon.org/resource/faldo
# imports: http://spinrdf.org/sp
# imports: http://spinrdf.org/spin

@prefix :        <http://example.org/insdc#> .
@prefix faldo:   <http://biohackathon.org/resource/faldo#> .
@prefix insdc:   <http://insdc.org/owl/> .
@prefix owl:     <http://www.w3.org/2002/07/owl#> .
@prefix rdf:     <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#> .
@prefix sp:      <http://spinrdf.org/sp#> .
@prefix spin:    <http://spinrdf.org/spin#> .
@prefix xsd:     <http://www.w3.org/2001/XMLSchema#> .

"""

s = """
<_:%d> a faldo:Region ;
           faldo:begin <_:%db> ;
           faldo:end <_:%de> ;
	   faldo:reference ddbj:%s .

<_:%db> a faldo:Position ; 
           a faldo:ExactPosition ;
           a faldo:%s ;
            faldo:position "%d"^^xsd:integer .
            

<_:%de> a faldo:Position ; 
           a faldo:ExactPosition ;
           a faldo:%s ;
            faldo:position "%d"^^xsd:integer 


"""
import sys
i = 0
for l in open(sys.argv[1], "r"):
	(chrx, x, x, start, end, x, strand, x, info) = l.strip().split("\t")
	gene = info.replace("transcript=", "").replace(";value=0","")
	if strand == "+":
		strand = "ForwardStrandPosition"
	else:
		strand = "ReverseStrandPosition"
	print s % (i,i,i,gene,i,strand,int(start),i,strand,int(end))
	i += 1 
