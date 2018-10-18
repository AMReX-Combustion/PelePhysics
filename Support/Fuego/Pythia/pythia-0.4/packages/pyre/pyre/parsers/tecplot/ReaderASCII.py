#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        (C) 1998-2001 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 

class ReaderASCII(object):


    def read(self, input):
        self._info.log("reading from '%s' -- assuming ASCII Tecplot format" % input.name)

        while 1:
            line = input.readline()
            if not line:
                break
            
            tokens = line.split()

            name = tokens[0].lower()
            if name == "title":
                value = line.split(" = ")[1]
                self.title = value.strip()
                self._info.log("title = {%s}" % self.title)
            elif name == "variables":
                value = line.split(" = ")[1]
                self.variables = value.split()
                self._info.log("variables = %r" % self.variables)
            elif name == "zone":
                self._zone(line, input)
            else:
                import journal
                firewall = journal.firewall("tecplot")
                firewall.log("unsupported tecplot tag '%s'" % name)
                return
                
        return
                

    def __init__(self):

        self.title = ""
        self.variables = ""
        self.nodes = 0
        self.simplices = 0
        
        self.mesh = None
        self.fields = None

        self._typesET = ["TRIANGLE", "TETRAHEDRON"]

        return


    def _zone(self, line, file):

        tokens = line[4:].strip().split(",")

        nodes = 0
        simplices = 0

        for token in tokens:
            if not token:
                continue
            
            variable, value = token.split("=", 1)
            variable = variable.strip().upper()
            if variable == "T":
                title = value.strip()
            elif variable == "N":
                nodes = int(value)
            elif variable == "E":
                simplices = int(value)
            elif variable == "F":
                if value.upper() != "FEPOINT":
                    import journal
                    error = journal.error("tecplot")
                    error.log("F type '%s' not supported" % value)
            elif variable == "ET":
                value = value.upper()
                if value not in self._typesET:
                    import journal
                    error = journal.error("tecplot")
                    error.log("ET type '%s' not supported" % value)

        self._info.log("found %d nodes and %d triangles" % (nodes, simplices))

        self.nodes = nodes
        self.simplices = simplices

        self.fields = self._readFields(file)
        self.mesh = self._readMesh(file)

        return


    def _readFields(self, file):
        from Nodal import Nodal

        nodal = Nodal()
        nodal.read(self.nodes, self.variables, file)
        return nodal.fields


    def _readMesh(self, file):
        from Mesh import Mesh

        mesh = Mesh()
        mesh.read(self.simplices, file)
        return mesh.simplices


    import journal
    _info = journal.debug("tecplot")


# version
__id__ = "$Id$"

#  End of file 
