#!/bin/bash

if [ -z "$1" ]
  then
    echo "Error: no coordinates file supplied."
    echo "   usage: plot_positions_in_Chimera.sh path_to_coordinates_file"
fi

if [ ! -f "$1" ]
  then
    echo "Error: coordinates file $1 not found"
else
  cp $1 .coordinates_file.tmp
  chimera --script ./ChimeraScripts/plot_positions.py && rm .coordinates_file.tmp
fi

