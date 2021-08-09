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
    from .Whitespace import Whitespace
    from .Comments import Comments

    from .ElementSection import ElementSection
    from .SpeciesSection import SpeciesSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .ReactionSection import ReactionSection

    tokenClasses = [
        Whitespace, Comments,
        ElementSection, SpeciesSection, QssSpeciesSection, ThermoSection, TransSection, ReactionSection
        ]

    return tokenClasses


def elementTokenClasses():
    from .Whitespace import Whitespace
    from .Comments import Comments
    from .EndSection import EndSection

    from .ElementSection import ElementSection
    from .SpeciesSection import SpeciesSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .ReactionSection import ReactionSection

    from .ElementName import ElementName

    tokenClasses = [
        Whitespace, Comments,
        EndSection, ElementSection, SpeciesSection, QssSpeciesSection, ThermoSection, TransSection, ReactionSection,
        ElementName
        ]

    return tokenClasses


def speciesTokenClasses():

    from .Whitespace import Whitespace
    from .Comments import Comments
    from .EndSection import EndSection

    from .ElementSection import ElementSection
    from .SpeciesSection import SpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .ReactionSection import ReactionSection

    from .SpeciesName import SpeciesName

    tokenClasses = [
        Whitespace, Comments,
        EndSection, ElementSection, SpeciesSection, ThermoSection, TransSection, ReactionSection, SpeciesName
        ]

    return tokenClasses

def qss_speciesTokenClasses():

    from .Whitespace import Whitespace
    from .Comments import Comments
    from .EndSection import EndSection

    from .ElementSection import ElementSection
    #from SpeciesSection import SpeciesSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .ReactionSection import ReactionSection

    #from SpeciesName import SpeciesName
    from .QssSpeciesName import QssSpeciesName

    tokenClasses = [
        Whitespace, Comments,
        EndSection, ElementSection, QssSpeciesSection, ThermoSection, TransSection, ReactionSection, QssSpeciesName
        ]

    return tokenClasses


def thermoTokenClasses():

    from .Whitespace import Whitespace
    from .Comments import Comments
    from .EndSection import EndSection

    from .ElementSection import ElementSection
    from .SpeciesSection import SpeciesSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .ReactionSection import ReactionSection

    from .TemperatureRange import TemperatureRange
    from .ThermoLine import ThermoLine

    tokenClasses = [
        EndSection, ElementSection, SpeciesSection, QssSpeciesSection, ThermoSection, TransSection, ReactionSection,
        TemperatureRange, ThermoLine,
        Whitespace, Comments,
        ]

    return tokenClasses


def transTokenClasses():

    from .Whitespace import Whitespace
    from .Comments import Comments
    from .EndSection import EndSection

    from .ElementSection import ElementSection
    from .SpeciesSection import SpeciesSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .ReactionSection import ReactionSection

    from .TransLine import TransLine

    tokenClasses = [
        EndSection, ElementSection, SpeciesSection, QssSpeciesSection, ThermoSection, TransSection, 
        ReactionSection,
        TransLine,
        Whitespace, Comments,
        ]

    return tokenClasses


def reactionTokenClasses():

    from .Whitespace import Whitespace
    from .Comments import Comments
    from .EndSection import EndSection

    from .ElementSection import ElementSection
    from .SpeciesSection import SpeciesSection
    from .QssSpeciesSection import QssSpeciesSection
    from .ThermoSection import ThermoSection
    from .TransSection import TransSection
    from .ReactionSection import ReactionSection

    from .Reaction import Reaction
    from .ReactionUnitsA import ReactionUnitsA
    from .ReactionUnitsE import ReactionUnitsE
    from .ReactionDuplicate import ReactionDuplicate
    from .ReactionEfficiencies import ReactionEfficiencies
    from .ReactionHV import ReactionHV
    from .ReactionLOW import ReactionLOW
    from .ReactionFORD import ReactionFORD
    from .ReactionLT import ReactionLT
    from .ReactionRLT import ReactionRLT
    from .ReactionReverse import ReactionReverse
    from .ReactionSRI import ReactionSRI
    from .ReactionTROE import ReactionTROE

    tokenClasses = [
        Comments, Whitespace,
        EndSection, ElementSection, SpeciesSection, QssSpeciesSection, ThermoSection, TransSection, ReactionSection,
        ReactionDuplicate, ReactionLOW, ReactionSRI, ReactionTROE,
        ReactionUnitsA, ReactionUnitsE,
        ReactionHV, ReactionLT, ReactionRLT, ReactionReverse, ReactionFORD,
        ReactionEfficiencies, Reaction
        ]

    return tokenClasses
    

def parameterTokenClasses():

    from .Whitespace import Whitespace
    from .Comments import Comments
    from .ParameterList import ParameterList

    tokenClasses = [
        Comments, Whitespace, ParameterList
        ]

    return tokenClasses


# version
__id__ = "$Id$"

#
# End of file
