# Analytically Reduced NC12H26 Mechanism with QSS

## Create mechanism

4 files are needed

- Skeletal mechanism (here `skeletal.inp`)
- List of non QSS species (here `non_qssa_list.txt`)
- Thermodynamic file (here `therm.dat`)
- Transport file (here `tran.dat`)

Link the file in `make_mechanism.sh`

Execute `bash make_mechanism.sh`

## Integration method

Verified with GMRES, denseDirect. Does not work yet with analytical Jacobian

## Remove quadratic coupling

3 methods are available for treatment of the quadratic coupling between QSS species

### 1) Species method

The smallest possible set of species that needs to be removed to eliminate quadratic coupling is identified and removed from the QSS species list. To activate it, uncomment the lines under Method 1 in `make_mechanism.sh`

### 2) Reaction method

All reactions that create a quadratic coupling are removed by eliminating both the forward and backward reaction, even if only the forward or the backward reaction create a quadratic coupling. To activate it, uncomment the lines under Method 2 in `make_mechanism.sh`

### 3) Some reaction method

All reactions that create a quadratic coupling are removed by eliminating the forward and/or the backward reaction that generate quadratic coupling. To activate it, uncomment the lines under Method 3 in `make_mechanism.sh`

## Make the dependency graph

`bash make-viztool.sh`

This creates a dependency graph of the QSS species as well as a graph of the reactions that induce quadratic coupling. The graphs are outputed in  `qssaTools/output/`


