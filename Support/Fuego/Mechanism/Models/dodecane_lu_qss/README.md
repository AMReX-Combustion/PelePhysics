# Analytically Reduced NC12H26 Mechanism with QSS

## Limitations

Does not support yet analytical Jacobian for the integration method and will abort if analytical Jacobian is used.

## Create mechanism

1. Add the 4 necessary files
  - Skeletal mechanism (here `skeletal.inp`)
  - List of non QSS species (here `non_qssa_list.txt`)
  - Thermodynamic file (here `therm.dat`)
  - Transport file (here `tran.dat`)

2. Link the 4 files in `make-mechanism.sh`

3. Choose a linearization technique (See [QSS Documentation](https://pelephysics.readthedocs.io/en/latest/QSS.html)) by commenting and uncommenting the code under `#~~~~ Method 1`, `#~~~~ Method 2` or `# ~~~~ Method 3` in `make-mechanism.sh`

4. Execute `bash make-mechanism.sh`


## Make the dependency graph

This creates a dependency graph of the QSS species as well as a graph of the reactions that induce quadratic coupling. The graphs are outputed in  `qssaTools/output/`

1. Add the 4 necessary files
  - Skeletal mechanism (here `skeletal.inp`)
  - List of non QSS species (here `non_qssa_list.txt`)
  - Thermodynamic file (here `therm.dat`)
  - Transport file (here `tran.dat`)

2. Link the 4 files in `make-viztool.sh`

3. Execute `bash make-viztool.sh`



