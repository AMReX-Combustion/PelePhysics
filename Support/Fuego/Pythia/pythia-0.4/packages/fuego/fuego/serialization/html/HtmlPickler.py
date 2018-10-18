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

from weaver.mills.HTMLMill import HTMLMill


class HtmlPickler(HTMLMill):


    def __init__(self):
        HTMLMill.__init__(self)
        self.__html = None;
        self.__body = None;
        self.__head = None;
        return

    def _html(self):
        if not self.__html:
            self._write('<HTML>')
            self.__html = 1;
        else:
            self._write('</HTML>')
            self.__html = None;


    def _head(self):
        if not self.__head:
            self._write('<HEAD>')
            self.__head = 1;
        else:
            self._write('</HEAD>')
            self.__head = None;

    def _body(self):
        if not self.__body:
            self._write('<BODY aLink=#9933ff bgColor=#000000 link=#ffffff text=#ffffff vLink=#33eeff>')
            self.__body = 1;
        else:
            self._write('</BODY>')
            self.__body = None;

    def _style(self):
        import style
        self._write( style.string() )

    def _speciesScript(self):
        import speciesScript
        self._write( speciesScript.string() )

    def _renderDocument(self, mechanism, options=None):

        self._rep += ['']

        self._html()
        self._head()
        self._style()
        self._speciesScript()        
        self._head()
        self._body()
        self._indent()

        self._renderElements(mechanism)
        self._renderSpecies(mechanism)
        self._renderReactions(mechanism)
        self._outdent()
   
        self._body()
        self._html()
        self._rep += ['']

        return


    def _renderElements(self, mechanism):
        self._write(self.line('elements'))

        self._write( '<H2>Element Section</H2>' )

        self._write('<table width="400" border = "1">')
        for element in mechanism.element():
            symbol = element.symbol
            id = element.id+1
            self._indent()
            self._write('<tr>')           
            self._indent()
            self._write('<td width="20"> %d </td>' % id) 
            self._write('<td> %s </td>' % symbol) 
            self._outdent()
            self._write('</tr>')           
            self._outdent()

        self._write('</table>')
        self._rep.append('')
            
        return


    def _speciesEvent(self, species):
        """ Only for NASA polynomials with tlow, tmid, thigh """

        str = ''
        str += 'onClick=speciesWin("%s",' % species.symbol
        for element, coefficient in species.composition:
           str += '"%s",%d,' % (element, coefficient)
        
        thermo = species.thermo
        if thermo:
            for model in thermo:
                for p in model.parameters:
                    str += '%g,' % p
                str += '%g,%g,' % (model.lowT, model.highT)

        
        str = str[:-1] + ")"
        return str

    def _renderSpecies(self, mechanism):
        nColumns = 6;
        self._write(self.line('species'))

        self._write( '<H2>Species Section</H2>' )
        
        self._write('<table align=center border=0 cellSpacing=1')
        self._write('<TR>')
        for dum in range(0,nColumns):
            self._write('<TD width=80></TD>')
        self._write('<TR>')

        col = 1
        for species in mechanism.species():
            symbol = species.symbol
            phase = species.phase
            id = species.id+1

            if col == 1:
                self._write('<tr valign=top>')

            self._indent()
            event = self._speciesEvent(species)
            self._write('<td width=80 bgcolor=#ffffff align=left height=30 vAlign=center><small>')
            self._write('<font color="#000000">%d</font color></small>' % id)
            self._write('<center> <a %s > <H2><font color="#0000FF"> %s</font color></H2> </a></center></td>' % (event, symbol) )
            self._outdent()
            
            if col == 6:
                self._write('</tr>')
                col = 1
            else:
                col = col+1

        self._write('</table>')


    def _renderReactions(self, mechanism):
        self._write(self.line('reactions'))

        self._write( '<H2>Reactions Section</H2>' )
        
        self._write('<table width="600" border = "1">')
        
        i = 0
        for reaction in mechanism.reaction():

            i += 1
            self._indent()
            self._write('<tr>')
            self._indent()
            self._write('<td width="20">%d</td>' % i)
            self._write('<td>%s</td>' % reaction.equation())
            self._outdent()
            self._write('</tr>')
            self._outdent()
            
        self._write('</table>')    

        return

    def _br(self):
        self_write('<br>')
    
    
#  End of file 
