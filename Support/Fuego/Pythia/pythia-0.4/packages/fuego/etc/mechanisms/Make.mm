# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

PROJECT = fuego
PACKAGE = mechanisms

#--------------------------------------------------------------------------
#

all: export


#--------------------------------------------------------------------------
# export

EXPORT_ETC = \
    therm.dat \
    GRIMech-1.2.ck2 \
    GRIMech-2.11.ck2 \
    GRIMech-3.0.ck2 \
    HMX-CIT-20010831.ck2 \
    HMX_20010423b.ck2 \
    HydrogenOxygen.ck2 \
    Konnov-0.4.ck2 \
    Leeds-CH4-1.2.ck2 \
    Leeds-CH4-1.3.ck2 \
    Leeds-S.ck2 \
    Maas-Pope.ck2 \
    Maas-Warnatz.ck2 \
    Maas-Warnatz_noN.ck2 \
    OxygenDissociation.ck2 \
    RDX.ck2 \
    allen.ck2 \
    arho.ck2 \
    bad-0001.ck2 \
    baulch.ck2 \
    benzene.ck2 \
    curran.ck2 \
    fujii.ck2 \
    glassman.ck2 \
    lindst.ck2 \
    lindstedt.ck2 \
    lutz88.ck2 \
    mill-bow89A.ck2 \
    mill-bow89AB.ck2 \
    mill-bow89B.ck2 \
    mill_smook.ck2 \
    oran.ck2 \
    reactions-0001.ck2 \
    reactions-0002.ck2 \
    reactions-0003.ck2 \
    reactions-0004.ck2 \
    reactions-0005.ck2 \
    reactions-0006.ck2 \
    reactions-0007.ck2 \
    reactions-0008.ck2 \
    reactions-0009.ck2 \
    reactions-0010.ck2 \
    reactions-0011.ck2 \
    reactions-0012.ck2 \
    reactions-0013.ck2 \
    simple.ck2 \
    test-mech.ck2 \
    wang_99.ck2 \
    ybco.ck2 \
    yetter.ck2 \


export:: export-etc

# version
# $Id$

#
# End of file
