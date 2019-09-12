#! /bin/bash
##############################################################
# Program: x.repack.sh
# Programmer: Robert Seigel
#             University of Miami
#             rseigel@rsmas.miami.edu
#             16 March 2013
# Purpose: This script repacks the HDF5 files while a distributive
#          memory RAMS run is occurring. Because Parallel HDF5 is
#          not yet supporting online compression, this script 
#          takes its place to conserve disk space.

##############################################################
# FUNCTION to find the new file and repack it according to 
# specified gzip level
repack ()
{
# Loop over .h5 files in requested directory
#   for f in $thisdir*'.h5'
files=$(find $mydir -maxdepth 1 -type f -name '*.h5' -size +2G)

   for f in $files 
   do
echo $f
         # Wait UNTIL header file exists to make sure HDF5 write is done
         #if [ $runtype == "analysis" ]; then
            hfile=${f:0:$((${#f}-5))}'head.txt'

            until [ -e $hfile ]; do sleep 1; done
            sleep 3
            echo 'hfile exists'
         #fi

         # Repack the file
         h5repack -f SHUF -f GZIP=$gziplevel $f $f.temp   #***UNCOMMENT
         mv $f.temp $f                                    #***UNCOMMENT

   done
}
######################################################################################
######################################################################################
# ONLYEDIT THESE

# Define some parameters
gziplevel=6        # gzip level
checkint=600         # in seconds
#cumulusdir='/Volumes/data/Cumulus/'$1'.' # Directory for file transfers

# END ONLYEDIT
######################################################################################
######################################################################################

##############################################################
# Execute repacking of HDF5 files
mydir=$1
      # Copy RAMSIN to the directory on cumulus
      #scp RAMSIN rseigel@cumulus.rsmas.miami.edu:$cumulusdir 

      # Continue to execute repack if (1) the model is running OR 
      #                               (2) the number of files have changed since
      #                                   last repack, i.e. new files were generated
#      postfilecount=`ls -1 $1*'.h5' | wc -l` 
#      echo $postfilecount
#      while [[ $2 -gt $postfilecount ]]
#      do   
        
         # Execute the repack function.
#         repack

         # Regrab file count. If this number is different than prefilecount, the 
         # repack function will execute again
#         postfilecount=`ls -1 $1*'.txt' | wc -l`         
#         echo $postfilecount
         # execute time controller
#         sleep $checkint
#      done

      repack

      # We are done
      echo $1" finished and x.repack.sh is terminating"

exit

###############################################################


