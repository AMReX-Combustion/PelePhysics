#!/usr/bin/env python

__doc__ = """ JavaScript for the species section 
for use by HtmlPickler """

def string():
    return scriptString

scriptString = """
<!-- from speciesScript.py -->
<script language="JavaScript">
<!-- hide

function speciesObject(arg) {
 

   // constructor should be called by
   // eg: "H2O", "H", 2, "O", 1, 7 nums, lowT, highT, 7 nums, lowT, highT

   //var arg = speciesObject.arguments;
   
   this._rep = ""
   this.speciesName = arg[0];
   nel = (arg.length-19)/2;
   this.AtomN = new Array(nel);
   this.AtomQ = new Array(nel);
   ind = 1;
   for (i=0; i<nel; i++) {
      this.AtomN[i] = arg[ind++];
      this.AtomQ[i] = arg[ind++];
   }
   this.nel = nel;
   this.a1 = arg[ind++];
   this.a2 = arg[ind++];
   this.a3 = arg[ind++];
   this.a4 = arg[ind++];
   this.a5 = arg[ind++];
   this.a6 = arg[ind++];
   this.a7 = arg[ind++];
   this.aL = arg[ind++];
   this.aH = arg[ind++];
   this.b1 = arg[ind++];
   this.b2 = arg[ind++];
   this.b3 = arg[ind++];
   this.b4 = arg[ind++];
   this.b5 = arg[ind++];
   this.b6 = arg[ind++];
   this.b7 = arg[ind++];
   this.bL = arg[ind++];
   this.bH = arg[ind++];
    

   function make() {
      this._rep = "<H2> Species Name: " + this.speciesName + "</H2>";
      this._rep += "Composition";
      this._rep += "<ul>";
      for (i=0; i<this.nel; i++) {
         this._rep += "<li>"+this.AtomN[i] + ": " + this.AtomQ[i] ;
      }
      this._rep += "</ul>";

      this._rep += "<H3>Coefficients for NASA Polynomial</H3>";
      this._rep += "T = [" + this.aL + "," + this.aH + "]"
      this._rep += "<ul>"
      this._rep += "<li> a1 : " + this.a1
      this._rep += "<li> a2 : " + this.a2
      this._rep += "<li> a3 : " + this.a3
      this._rep += "<li> a4 : " + this.a4
      this._rep += "<li> a5 : " + this.a5
      this._rep += "<li> a6 : " + this.a6
      this._rep += "<li> a7 : " + this.a7
      this._rep += "</ul>"
      this._rep += "T = [" + this.bL + "," + this.bH + "]"
      this._rep += "<ul>"
      this._rep += "<li> a1 : " + this.b1
      this._rep += "<li> a2 : " + this.b2
      this._rep += "<li> a3 : " + this.b3
      this._rep += "<li> a4 : " + this.b4
      this._rep += "<li> a5 : " + this.b5
      this._rep += "<li> a6 : " + this.b6
      this._rep += "<li> a7 : " + this.b7
      this._rep += "</ul>"
   }


   function write() {
        this._make(); 
        return this._rep;
   }

   this._make = make;
   this.write = write;
}


function speciesWin() {

// forward arguments to species constructor
var species = new speciesObject(speciesWin.arguments);

myWin = open("", "displayWindow", "width=500");
myWin.document.open();
myWin.document.write(species.write());
myWin.document.close();
}

// -->
</script>
"""

# End of file
