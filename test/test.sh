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

CHECK_LIB=$(pkg-config --libs edipack2)
if [ -z ${CHECK_LIB} ]
then
    echo "\e[31m ERROR \e[0m"
    echo " EDIpack2 not loaded"
    return 1
fi

WITH_MPI=$(pkg-config --variable=mpi edipack2)
BUILD_TYPE=$(pkg-config --variable=build_type edipack2)

cd bin/
HERE=`pwd`
echo $HERE

while read DIR; do
    if [ -d $DIR ]; then
	if [ -f $DIR/*.x ] ; then
	    echo "TESTING $DIR"
	    cd $DIR
	    pwd
	    for exe in *.x
	    do
		echo "Running $exe:"
		if [ "$BUILD_TYPE" = "DEBUG" ]
		then
		    echo "./$exe ED_VERBOSE=3 LOGFILE=6"
		    ./$exe ED_VERBOSE=1 LOGFILE=6
		    echo ""
		    echo ""
		    sleep 1
		else
		    if [ -z ${WITH_MPI} ]
		    then
			echo "./$exe ED_VERBOSE=3 LOGFILE=6"
			./$exe ED_VERBOSE=1 LOGFILE=6
			echo ""
			echo ""
			sleep 1
		    else
			echo "mpiexec -np 2 ./$exe ED_VERBOSE=3 LOGFILE=6"
			mpiexec -np 2 ./$exe ED_VERBOSE=1 LOGFILE=6 < /dev/null
			echo ""
			echo ""
			sleep 1
		    fi
		fi
	    done
	    cd $HERE
	fi
    fi
done<list_dir
    

