import ceptr.writer as cw


def sparsity(fstream, species_info):
    nSpecies = species_info.nSpecies

    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("compute the sparsity pattern of the chemistry Jacobian"),
    )
    cw.writer(
        fstream,
        "void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::GpuArray<amrex::Real,%d> Jac = {0.0};" % (nSpecies + 1) ** 2,
    )
    cw.writer(
        fstream, "amrex::GpuArray<amrex::Real,%d> conc = {0.0};" % (nSpecies)
    )
    cw.writer(fstream, "for (int n=0; n<%d; n++) {" % (nSpecies))
    cw.writer(fstream, "    conc[n] = 1.0/ %f ;" % (nSpecies))
    cw.writer(fstream, "}")
    cw.writer(fstream, "aJacobian(&Jac[0], &conc[0], 1500.0, *consP);")
    cw.writer(fstream)

    cw.writer(fstream, "int nJdata_tmp = 0;")
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "if(Jac[ %d * k + l] != 0.0){" % (nSpecies + 1))

    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "*nJdata = NCELLS * nJdata_tmp;")

    cw.writer(fstream, "}")
    cw.writer(fstream)
    cw.writer(fstream)

    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment("compute the sparsity pattern of the system Jacobian"),
    )
    cw.writer(
        fstream,
        "void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::GpuArray<amrex::Real,%d> Jac = {0.0};" % (nSpecies + 1) ** 2,
    )
    cw.writer(
        fstream, "amrex::GpuArray<amrex::Real,%d> conc = {0.0};" % (nSpecies)
    )
    cw.writer(fstream, "for (int n=0; n<%d; n++) {" % (nSpecies))
    cw.writer(fstream, "    conc[n] = 1.0/ %f ;" % (nSpecies))
    cw.writer(fstream, "}")
    cw.writer(fstream, "aJacobian(&Jac[0], &conc[0], 1500.0, *consP);")
    cw.writer(fstream)

    cw.writer(fstream, "int nJdata_tmp = 0;")
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "if(k == l){")

    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1;")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "if(Jac[ %d * k + l] != 0.0){" % (nSpecies + 1))

    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "*nJdata = NCELLS * nJdata_tmp;")

    cw.writer(fstream, "}")
    cw.writer(fstream)
    cw.writer(fstream)

    cw.writer(fstream)
    cw.writer(
        fstream,
        cw.comment(
            "compute the sparsity pattern of the simplified (for preconditioning) system Jacobian"
        ),
    )
    cw.writer(
        fstream,
        "void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::GpuArray<amrex::Real,%d> Jac = {0.0};" % (nSpecies + 1) ** 2,
    )
    cw.writer(
        fstream, "amrex::GpuArray<amrex::Real,%d> conc = {0.0};" % (nSpecies)
    )
    cw.writer(fstream, "for (int n=0; n<%d; n++) {" % (nSpecies))
    cw.writer(fstream, "    conc[n] = 1.0/ %f ;" % (nSpecies))
    cw.writer(fstream, "}")
    cw.writer(fstream, "aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);")
    cw.writer(fstream)

    cw.writer(fstream, "int nJdata_tmp = 0;")
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "if(k == l){")

    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1;")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "if(Jac[ %d * k + l] != 0.0){" % (nSpecies + 1))

    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream)
    cw.writer(fstream, "nJdata[0] = nJdata_tmp;")

    cw.writer(fstream, "}")
    cw.writer(fstream)
    cw.writer(fstream)

    cw.writer(
        fstream,
        cw.comment(
            "compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0"
        ),
    )
    cw.writer(
        fstream,
        "void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::GpuArray<amrex::Real,%d> Jac = {0.0};" % (nSpecies + 1) ** 2,
    )
    cw.writer(
        fstream, "amrex::GpuArray<amrex::Real,%d> conc = {0.0};" % (nSpecies)
    )
    cw.writer(fstream, "for (int n=0; n<%d; n++) {" % (nSpecies))
    cw.writer(fstream, "    conc[n] = 1.0/ %f ;" % (nSpecies))
    cw.writer(fstream, "}")
    cw.writer(fstream, "aJacobian(&Jac[0], &conc[0], 1500.0, *consP);")
    cw.writer(fstream)

    cw.writer(fstream, "colPtrs[0] = 0;")
    cw.writer(fstream, "int nJdata_tmp = 0;")
    cw.writer(fstream, "for (int nc=0; nc<NCELLS; nc++) {")

    cw.writer(fstream, "int offset_row = nc * %d;" % (nSpecies + 1))
    cw.writer(fstream, "int offset_col = nc * %d;" % (nSpecies + 1))
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "if(Jac[%d*k + l] != 0.0) {" % (nSpecies + 1))

    cw.writer(fstream, "rowVals[nJdata_tmp] = l + offset_row; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream, "colPtrs[offset_col + (k + 1)] = nJdata_tmp;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream)

    cw.writer(
        fstream,
        cw.comment(
            "compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0"
        ),
    )
    cw.writer(
        fstream,
        "void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::GpuArray<amrex::Real,%d> Jac = {0.0};" % (nSpecies + 1) ** 2,
    )
    cw.writer(
        fstream, "amrex::GpuArray<amrex::Real,%d> conc = {0.0};" % (nSpecies)
    )
    cw.writer(fstream, "for (int n=0; n<%d; n++) {" % (nSpecies))
    cw.writer(fstream, "    conc[n] = 1.0/ %f ;" % (nSpecies))
    cw.writer(fstream, "}")
    cw.writer(fstream, "aJacobian(&Jac[0], &conc[0], 1500.0, *consP);")
    cw.writer(fstream)

    cw.writer(fstream, "if (base == 1) {")

    cw.writer(fstream, "rowPtrs[0] = 1;")
    cw.writer(fstream, "int nJdata_tmp = 1;")
    cw.writer(fstream, "for (int nc=0; nc<NCELLS; nc++) {")

    cw.writer(fstream, "int offset = nc * %d;" % (nSpecies + 1))
    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "if(Jac[%d*k + l] != 0.0) {" % (nSpecies + 1))

    cw.writer(fstream, "colVals[nJdata_tmp-1] = k+1 + offset; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream, "rowPtrs[offset + (l + 1)] = nJdata_tmp;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "rowPtrs[0] = 0;")
    cw.writer(fstream, "int nJdata_tmp = 0;")
    cw.writer(fstream, "for (int nc=0; nc<NCELLS; nc++) {")

    cw.writer(fstream, "int offset = nc * %d;" % (nSpecies + 1))
    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "if(Jac[%d*k + l] != 0.0) {" % (nSpecies + 1))

    cw.writer(fstream, "colVals[nJdata_tmp] = k + offset; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream, "rowPtrs[offset + (l + 1)] = nJdata_tmp;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream)

    cw.writer(
        fstream,
        cw.comment("compute the sparsity pattern of the system Jacobian"),
    )
    cw.writer(fstream, cw.comment("CSR format BASE is user choice"))
    cw.writer(
        fstream,
        "void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::GpuArray<amrex::Real,%d> Jac = {0.0};" % (nSpecies + 1) ** 2,
    )
    cw.writer(
        fstream, "amrex::GpuArray<amrex::Real,%d> conc = {0.0};" % (nSpecies)
    )
    cw.writer(fstream, "for (int n=0; n<%d; n++) {" % (nSpecies))
    cw.writer(fstream, "    conc[n] = 1.0/ %f ;" % (nSpecies))
    cw.writer(fstream, "}")
    cw.writer(fstream, "aJacobian(&Jac[0], &conc[0], 1500.0, *consP);")
    cw.writer(fstream)

    cw.writer(fstream, "if (base == 1) {")

    cw.writer(fstream, "rowPtr[0] = 1;")
    cw.writer(fstream, "int nJdata_tmp = 1;")
    cw.writer(fstream, "for (int nc=0; nc<NCELLS; nc++) {")

    cw.writer(fstream, "int offset = nc * %d;" % (nSpecies + 1))
    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "if (k == l) {")

    cw.writer(fstream, "colVals[nJdata_tmp-1] = l+1 + offset; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "if(Jac[%d*k + l] != 0.0) {" % (nSpecies + 1))

    cw.writer(fstream, "colVals[nJdata_tmp-1] = k+1 + offset; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream, "rowPtr[offset + (l + 1)] = nJdata_tmp;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "rowPtr[0] = 0;")
    cw.writer(fstream, "int nJdata_tmp = 0;")
    cw.writer(fstream, "for (int nc=0; nc<NCELLS; nc++) {")

    cw.writer(fstream, "int offset = nc * %d;" % (nSpecies + 1))
    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "if (k == l) {")

    cw.writer(fstream, "colVals[nJdata_tmp] = l + offset; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "if(Jac[%d*k + l] != 0.0) {" % (nSpecies + 1))

    cw.writer(fstream, "colVals[nJdata_tmp] = k + offset; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream, "rowPtr[offset + (l + 1)] = nJdata_tmp;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream)

    cw.writer(
        fstream,
        cw.comment(
            "compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU"
        ),
    )
    cw.writer(fstream, cw.comment("BASE 0"))
    cw.writer(
        fstream,
        "void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::GpuArray<amrex::Real,%d> Jac = {0.0};" % (nSpecies + 1) ** 2,
    )
    cw.writer(
        fstream, "amrex::GpuArray<amrex::Real,%d> conc = {0.0};" % (nSpecies)
    )
    cw.writer(fstream, "for (int n=0; n<%d; n++) {" % (nSpecies))
    cw.writer(fstream, "    conc[n] = 1.0/ %f ;" % (nSpecies))
    cw.writer(fstream, "}")
    cw.writer(fstream, "aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);")
    cw.writer(fstream)

    cw.writer(fstream, "colPtrs[0] = 0;")
    cw.writer(fstream, "int nJdata_tmp = 0;")
    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "if (k == l) {")

    cw.writer(fstream, "rowVals[nJdata_tmp] = l; ")
    cw.writer(fstream, "indx[nJdata_tmp] = %d*k + l;" % (nSpecies + 1))
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "if(Jac[%d*k + l] != 0.0) {" % (nSpecies + 1))

    cw.writer(fstream, "rowVals[nJdata_tmp] = l; ")
    cw.writer(fstream, "indx[nJdata_tmp] = %d*k + l;" % (nSpecies + 1))
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream, "colPtrs[k+1] = nJdata_tmp;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream)

    cw.writer(
        fstream,
        cw.comment(
            "compute the sparsity pattern of the simplified (for precond) system Jacobian"
        ),
    )
    cw.writer(fstream, cw.comment("CSR format BASE is under choice"))
    cw.writer(
        fstream,
        "void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)",
    )
    cw.writer(fstream, "{")

    cw.writer(
        fstream,
        "amrex::GpuArray<amrex::Real,%d> Jac = {0.0};" % (nSpecies + 1) ** 2,
    )
    cw.writer(
        fstream, "amrex::GpuArray<amrex::Real,%d> conc = {0.0};" % (nSpecies)
    )
    cw.writer(fstream, "for (int n=0; n<%d; n++) {" % (nSpecies))
    cw.writer(fstream, "    conc[n] = 1.0/ %f ;" % (nSpecies))
    cw.writer(fstream, "}")
    cw.writer(fstream, "aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);")
    cw.writer(fstream)

    cw.writer(fstream, "if (base == 1) {")

    cw.writer(fstream, "rowPtr[0] = 1;")
    cw.writer(fstream, "int nJdata_tmp = 1;")
    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "if (k == l) {")

    cw.writer(fstream, "colVals[nJdata_tmp-1] = l+1; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "if(Jac[%d*k + l] != 0.0) {" % (nSpecies + 1))

    cw.writer(fstream, "colVals[nJdata_tmp-1] = k+1; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream, "rowPtr[l+1] = nJdata_tmp;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "rowPtr[0] = 0;")
    cw.writer(fstream, "int nJdata_tmp = 0;")
    cw.writer(fstream, "for (int l=0; l<%d; l++) {" % (nSpecies + 1))

    cw.writer(fstream, "for (int k=0; k<%d; k++) {" % (nSpecies + 1))

    cw.writer(fstream, "if (k == l) {")

    cw.writer(fstream, "colVals[nJdata_tmp] = l; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "} else {")

    cw.writer(fstream, "if(Jac[%d*k + l] != 0.0) {" % (nSpecies + 1))

    cw.writer(fstream, "colVals[nJdata_tmp] = k; ")
    cw.writer(fstream, "nJdata_tmp = nJdata_tmp + 1; ")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
    cw.writer(fstream, "rowPtr[l+1] = nJdata_tmp;")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")

    cw.writer(fstream, "}")
