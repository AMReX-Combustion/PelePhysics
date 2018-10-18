# Analytically Reduced Mechanisms
This mechanism is in a format that requires the use of the ckwc.f routines to evaluate the species production rates. These should NOT be evaluated via the reactions provided in the chemkin reaction data file.
Note that if a .c file is generated at this point, it will implement nothing because no reactions are provided in the .inp file.

NOTE: Old files have been kept, but should be used with caution since the transport part of the file is fishy. For example, fr the 17th species, WT & NLIN indicates it is N but no N is present in neither the .inp or ckwc.f files. Also it seems to require 19 specs in the ckwc.f routine but only 17 are considered in the transport routines... So: the potential user should figure out how to use those files properly.
