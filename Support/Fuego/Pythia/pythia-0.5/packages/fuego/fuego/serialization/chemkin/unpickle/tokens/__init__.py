#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#


from __future__ import absolute_import


def sectionsTokenClasses():
    from .Comments import Comments
    from .ElementSection import ElementSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ReactionSection import ReactionSection
    from .SpeciesSection import SpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .Whitespace import Whitespace

    tokenClasses = [
        Whitespace,
        Comments,
        ElementSection,
        SpeciesSection,
        QssSpeciesSection,
        ThermoSection,
        TransSection,
        ReactionSection,
    ]

    return tokenClasses


def elementTokenClasses():
    from .Comments import Comments
    from .ElementName import ElementName
    from .ElementSection import ElementSection
    from .EndSection import EndSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ReactionSection import ReactionSection
    from .SpeciesSection import SpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .Whitespace import Whitespace

    tokenClasses = [
        Whitespace,
        Comments,
        EndSection,
        ElementSection,
        SpeciesSection,
        QssSpeciesSection,
        ThermoSection,
        TransSection,
        ReactionSection,
        ElementName,
    ]

    return tokenClasses


def speciesTokenClasses():

    from .Comments import Comments
    from .ElementSection import ElementSection
    from .EndSection import EndSection
    from .ReactionSection import ReactionSection
    from .SpeciesName import SpeciesName
    from .SpeciesSection import SpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .Whitespace import Whitespace

    tokenClasses = [
        Whitespace,
        Comments,
        EndSection,
        ElementSection,
        SpeciesSection,
        ThermoSection,
        TransSection,
        ReactionSection,
        SpeciesName,
    ]

    return tokenClasses


def qss_speciesTokenClasses():

    from .Comments import Comments
    from .ElementSection import ElementSection
    from .EndSection import EndSection

    # from SpeciesName import SpeciesName
    from .QssSpeciesName import QssSpeciesName

    # from SpeciesSection import SpeciesSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ReactionSection import ReactionSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .Whitespace import Whitespace

    tokenClasses = [
        Whitespace,
        Comments,
        EndSection,
        ElementSection,
        QssSpeciesSection,
        ThermoSection,
        TransSection,
        ReactionSection,
        QssSpeciesName,
    ]

    return tokenClasses


def thermoTokenClasses():

    from .Comments import Comments
    from .ElementSection import ElementSection
    from .EndSection import EndSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ReactionSection import ReactionSection
    from .SpeciesSection import SpeciesSection
    from .TemperatureRange import TemperatureRange
    from .ThermoLine import ThermoLine
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .Whitespace import Whitespace

    tokenClasses = [
        EndSection,
        ElementSection,
        SpeciesSection,
        QssSpeciesSection,
        ThermoSection,
        TransSection,
        ReactionSection,
        TemperatureRange,
        ThermoLine,
        Whitespace,
        Comments,
    ]

    return tokenClasses


def transTokenClasses():

    from .Comments import Comments
    from .ElementSection import ElementSection
    from .EndSection import EndSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ReactionSection import ReactionSection
    from .SpeciesSection import SpeciesSection
    from .ThermoSection import ThermoSection
    from .TransLine import TransLine
    from .TransSection import TransSection
    from .Whitespace import Whitespace

    tokenClasses = [
        EndSection,
        ElementSection,
        SpeciesSection,
        QssSpeciesSection,
        ThermoSection,
        TransSection,
        ReactionSection,
        TransLine,
        Whitespace,
        Comments,
    ]

    return tokenClasses


def reactionTokenClasses():

    from .Comments import Comments
    from .ElementSection import ElementSection
    from .EndSection import EndSection
    from .QssSpeciesSection import QssSpeciesSection
    from .Reaction import Reaction
    from .ReactionDuplicate import ReactionDuplicate
    from .ReactionEfficiencies import ReactionEfficiencies
    from .ReactionFORD import ReactionFORD
    from .ReactionHV import ReactionHV
    from .ReactionLOW import ReactionLOW
    from .ReactionLT import ReactionLT
    from .ReactionReverse import ReactionReverse
    from .ReactionRLT import ReactionRLT
    from .ReactionSection import ReactionSection
    from .ReactionSRI import ReactionSRI
    from .ReactionTROE import ReactionTROE
    from .ReactionPLOG import ReactionPLOG
    from .ReactionUnitsA import ReactionUnitsA
    from .ReactionUnitsE import ReactionUnitsE
    from .SpeciesSection import SpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .Whitespace import Whitespace

    tokenClasses = [
        Comments,
        Whitespace,
        EndSection,
        ElementSection,
        SpeciesSection,
        QssSpeciesSection,
        ThermoSection,
        TransSection,
        ReactionSection,
        ReactionDuplicate,
        ReactionLOW,
        ReactionSRI,
        ReactionTROE,
        ReactionPLOG,
        ReactionUnitsA,
        ReactionUnitsE,
        ReactionHV,
        ReactionLT,
        ReactionRLT,
        ReactionReverse,
        ReactionFORD,
        ReactionEfficiencies,
        Reaction,
    ]

    return tokenClasses


def parameterTokenClasses():

    from .Comments import Comments
    from .ParameterList import ParameterList
    from .Whitespace import Whitespace

    tokenClasses = [Comments, Whitespace, ParameterList]

    return tokenClasses


# version
__id__ = "$Id$"

#
# End of file
