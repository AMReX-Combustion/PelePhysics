#!/usr/bin/env python
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003 All Rights Reserved
# 
#  <LicenseText>
# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 


from LineMill import LineMill

class PythonMill(LineMill):


    names = ["python"]


    def __init__(self):
        LineMill.__init__(self, "#", "#!/usr/bin/env python")
        return


    def _versionId(self):
        format = self.inventory.versionId

        if format:
            self._rep += [
                '', self.line(" version"), '__id__ = "' + format.strip() + '"'
                ]

        return


# version
__id__ = "$Id$"

#  End of file 
