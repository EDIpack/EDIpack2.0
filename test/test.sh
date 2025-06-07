#Code to run all tests in test/bin/
#N.B. test should end by .x

set -e

if [ ! -d "bin" ]
then
    echo "\e[31m ERROR \e[0m"
    echo " There is no *bin* directory"
    echo " Try  'make all' before testing"
    return 1
fi

CHECK_LIB=$(pkg-config --libs edipack)
if [ -z ${CHECK_LIB} ]
then
    echo "\e[31m ERROR \e[0m"
    echo " EDIpack not loaded"
    return 1
fi

WITH_MPI=$(pkg-config --variable=mpi edipack)
BUILD_TYPE=$(pkg-config --variable=build_type edipack)

cd bin/
HERE=`pwd`
echo $HERE

#FOR FUTURE REFERENCE
# IFS='_'
# read -r -a RAV <<< $VAR;
# BATH=${RAV[0]^^}
# MODE=${RAV[1]^^}
# INEQ=${RAV[2]^^}


#Clean the list_dir input list to avoid repetitions.
awk '!seen[$0]++' list_dir > tmp
mv tmp list_dir



while read DIR; do   
    if [ -d $DIR ]; then
    	if [ -f $DIR/*.x ] ; then
    	    echo "TESTING $DIR"
    	    cd $DIR
    	    pwd
    	    if [ -f DONE.out ];then
    		echo "Test $DIR has already been passed. Skip"
    	    else		
    		for exe in *.x
    		do
    		    echo "Running $exe:"
    		    if [ "$BUILD_TYPE" = "DEBUG" ]
    		    then
    			echo "./$exe ED_VERBOSE=3 LOGFILE=6"
    			./$exe ED_VERBOSE=1 LOGFILE=6
    			echo ""
    			echo ""
    			echo "leaving $DIR..."
    			echo ""
    			echo ""
    			touch DONE.out			
    		    else
    			if [ -z ${WITH_MPI} ]
    			then
    			    echo "./$exe ED_VERBOSE=3 LOGFILE=6"
    			    ./$exe ED_VERBOSE=1 LOGFILE=6
    			    echo ""
    			    echo ""
    			    echo "leaving $DIR..."
    			    echo ""
    			    echo ""
    			    touch DONE.out
    			else
    			    echo "mpiexec -np 2 ./$exe ED_VERBOSE=3 LOGFILE=6"
    			    mpiexec -np 2 ./$exe ED_VERBOSE=1 LOGFILE=6 < /dev/null
    			    echo ""
    			    echo ""
    			    echo "leaving $DIR..."
    			    echo ""
    			    echo ""
    			    touch DONE.out
    			fi
    		    fi
    		done
    	    fi
    	    cd $HERE
    	fi
    fi
done<list_dir
    

