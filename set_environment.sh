#!/bin/bash
# set_environment.sh - preparing framework environment under GPL v2+
#-------------------------------------------------------------------

export PPSFRAMEWORKROOT=`pwd`
export WORKINGPPSDIR=$PPSFRAMEWORKROOT/working

#if [-d "$WORKINGPPSDIR" -a ! -h "$WORKINGPPSDIR"]
if [ -d $WORKINGPPSDIR ]
then 
   echo -e "\n$WORKINGPPSDIR is present. Nothing is needed to do.\n"
else
   echo -e "\nCreating a working dir.\n"
   mkdir $PPSFRAMEWORKROOT/working
   cp $PPSFRAMEWORKROOT/Skimmer/test/RunMissingMassSearches.py $PPSFRAMEWORKROOT/working/.
   cp $PPSFRAMEWORKROOT/Skimmer/data/*.root $PPSFRAMEWORKROOT/working/.
   cp -r $PPSFRAMEWORKROOT/Skimmer/test/json_dpg $PPSFRAMEWORKROOT/working/.
   echo -e "\nChecking: $WORKINGPPSDIR\n"
fi

alias PREPAREMCGENERATION="cp $PPSFRAMEWORKROOT/MCProduction/Execution/*.sh $WORKINGPPSDIR/." 
