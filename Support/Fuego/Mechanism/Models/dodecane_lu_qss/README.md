# Analytically Reduced Mechanism
This mechanism is the QSS version.




## Create mechanism

`bash make_mechanism.sh`


## Integration method

Verified with GMRES, denseDirect. Does not work yet with analytical Jacobian

## Remove quadratic coupling

3 methods are available for treatment of the quadratic coupling between QSS species
Result of 0D validation are displayed over the range of applicability of the QSS assumption (![equation](http://www.sciweavers.org/tex2img.php?eq=%5Cphi%20%5Cin%20%5B0.5%2C%202%5D%2C%20T_0%20%5Cin%20%5B800K%2C%201600K%5D%2C%20P%20%5Cin%20%5B1atm%2C%2050atm%5D%20%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0))
Left plot shows the correlation with ignition time obtained from 343 calculations spanning the QSS assumption applicability interval.
Right plot shows the ignition time at conditions similar to the KPP problem.

### 1) Species method

The smallest possible set of species that needs to be removed to eliminate quadratic coupling is identified and removed from the QSS species list. To activate it, uncomment the lines under Method 1 in `make_mechanism.sh`

<p float="left">
  <img src="Figures/correlationMethod1.png" width="350"/> 
  <img src="Figures/KPPMethod1.png" width="350"/> 
</p>


### 2) Reaction method

All reactions that create a quadratic coupling are removed by eliminating both the forward and backward reaction, even if only the forward or the backward reaction create a quadratic coupling. To activate it, uncomment the lines under Method 2 in `make_mechanism.sh`


<p float="left">
  <img src="Figures/correlationMethod2.png" width="350"/> 
  <img src="Figures/KPPMethod2.png" width="350"/> 
</p>

### 3) Some reaction method

All reactions that create a quadratic coupling are removed by eliminating the forward and/or the backward reaction that generate quadratic coupling. To activate it, uncomment the lines under Method 3 in `make_mechanism.sh`


<p float="left">
  <img src="Figures/correlationMethod3.png" width="350"/> 
  <img src="Figures/KPPMethod3.png" width="350"/> 
</p>


## Make the dependency graph

`bash make-viztool.sh`

This creates a dependency graph of the QSS species (left) as well as a graph of the reactions that induce quadratic coupling (right). The graphs are outputed in  `qssaTools/output/`


<p float="left">
  <img src="Figures/directedGraphQSS.png" width="350"/>
  <img src="Figures/quadGraphQSS.png" width="350"/>
</p>


