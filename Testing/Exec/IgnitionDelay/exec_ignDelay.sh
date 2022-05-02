#!/bin/bash

# Find first executable name available
execname=`find . -name "Pele*.ex" | head -1`

if [ -z "$execname" ]
then
    echo ERROR: No executable found, cannot compute ignition delay
    echo Compile PelePhysics first
    exit 1
else
    if [[ -f "PPreaction.txt" ]]; then
        rm PPreaction.txt
    fi
    if [[ -f "log" ]]; then
        rm log
    fi

    echo INFO: Using $execname to compute ignition delay
    # Rough estimate of ignition time
    $execname inputs/inputs.0d_firstpass
    python computeIgnitionDelay.py -v -est -f inputs/inputs.0d_firstpass
    
    # Refined estimate of ignition time
    $execname inputs/inputs.0d_refine
    python computeIgnitionDelay.py -ref -f inputs/inputs.0d_refine
fi


