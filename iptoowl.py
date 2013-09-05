import sys
import os
import xml.etree.ElementTree as ET

tag = '{http://www.ebi.ac.uk/Tools/common/schema}'

inroot = ET.parse(sys.argv[1]).getroot()
outfile = sys.argv[2]
tempfile = sys.argv[2] + ".temp"
outroot = ET.Element('this')
tree = ET.ElementTree(outroot)

for ip in inroot.iter(tag + 'protein'):
    protein = ip.attrib["id"]
    a = ET.SubElement(outroot, 'Declaration')
    b = ET.SubElement(a, "NamedIndividual", {"IRI": "#" + protein})
    a = ET.SubElement(outroot, 'ClassAssertion')
    b = ET.SubElement(a, "Class", {"IRI": "#Protein"})
    c = ET.SubElement(a, "NamedIndividual", {"IRI": "#" + protein})
    for ip in inroot.iter(tag + 'interpro'):
        if ip.attrib["type"] == "Domain" or ip.attrib["type"] == "Repeat":
            a = ET.SubElement(outroot, 'ClassAssertion')
            b = ET.SubElement(a, 'ObjectSomeValuesFrom')
            c = ET.SubElement(b, 'ObjectProperty', {"IRI": '#hasDomain'})
            d = ET.SubElement(b, 'Class', {"IRI": '#' + ip.attrib['id']})
            e = ET.SubElement(a, 'NamedIndividual', {"IRI": '#' + protein})

tree.write(tempfile)
os.system("sed 's/<this>//g' " + tempfile +
          "| sed 's/<\/this>//g' >> " + outfile)
os.system("rm " + tempfile)


"""
     <Declaration>
        <NamedIndividual IRI="#prova"/>
    </Declaration>


<ClassAssertion>
        <Class IRI="#Protein"/>
        <NamedIndividual IRI="#prova"/>
    </ClassAssertion>


    <ClassAssertion>
        <ObjectSomeValuesFrom>
            <ObjectProperty IRI="#hasDomain"/>
            <Class IRI="#Beta-Casp"/>
        </ObjectSomeValuesFrom>
        <NamedIndividual IRI="#prova"/>
    </ClassAssertion>

"""
