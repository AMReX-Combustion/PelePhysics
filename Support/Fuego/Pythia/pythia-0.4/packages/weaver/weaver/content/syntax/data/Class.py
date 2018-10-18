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


from Type import Type


class Class(Type):


    def identify(self, auth):
        return auth.onClass(self)


    # accessors to the various components
    def extends(self):
        return tuple(self._extends)
    

    def implements(self):
        return tuple(self._implements)


    def classAttributes(self):
        return tuple(self._classAttributes)
    

    def instanceAttributes(self):
        return tuple(self._instanceAttributes)
    

    def extend(self, ancestor):
        self._extends.append(ancestor)
        return
    

    def implement(self, ancestor):
        self._implements.append(ancestor)
        return


    def classAttribute(self, variable):
        return self._classAttributes.append(variable)
    

    def instanceAttribute(self, variable):
        return self._instanceAttributes.append(variable)
    

    def __init__(self, name):
        Type.__init__(self, name)

        self._extends = []
        self._implements = []

        self._metaMethods = []
        self._classMethods = []
        self._instanceMethods = []

        self._classAttributes = []
        self._instanceAttributes = []

        return


    


# version
__id__ = "$Id$"

#  End of file 
