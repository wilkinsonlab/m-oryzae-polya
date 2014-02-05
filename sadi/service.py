import sadi,os
from rdflib import *

fastq_mcf=Namespace("file:///home/marco/Downloads/_sadi/fastq_mcf.owl#")
foaf=Namespace("http://xmlns.com/foaf/0.1/")

class ExampleService(sadi.Service):

    def getInputClass(self):
	print "input"
        return fastq_mcf.FASTQ_raw

    def getOutputClass(self):
	print "output"
        return fastq_mcf.FASTQ_processed

    def process(self, input, output):
	print "process"
	output.set(fastq_mcf.fastq_mcf, Literal("running fastq-mcf on " + input.value(foaf.name).value))



resource = ExampleService()

if __name__ == "__main__":
    sadi.serve(resource, port=9090)
