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
		./$exe ED_VERBOSE=1 LOGFILE=6 
		echo ""
		echo ""
		sleep 1
	    done
	    cd $HERE
	fi
    fi
done<list_dir
    

