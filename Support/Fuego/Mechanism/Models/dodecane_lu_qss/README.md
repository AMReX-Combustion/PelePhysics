# Analytically Reduced Mechanisms
!!This mechanism is the QSS version. !!


## Remove quadratic coupling

3 methods are available for treatment of the quadratic coupling between QSS species

### 1) Species method

The smallest possible set of species that needs to be removed to eliminate quadratic coupling is identified and removed from the QSS species list. To activate it, uncomment the lines under Method 1 in `make_mechanism.sh`

### 2) Reaction method

All reactions that create a quadratic coupling are removed by eliminating both the forward and backward reaction, even if only the forward or the backward reaction create a quadratic coupling. To activate it, uncomment the lines under Method 2 in `make_mechanism.sh`

### 3) Some reaction method

All reactions that create a quadratic coupling are removed by eliminating the forward and/or the backward reaction that generate quadratic coupling. To activate it, uncomment the lines under Method 3 in `make_mechanism.sh`

