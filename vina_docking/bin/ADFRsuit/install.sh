#!/bin/sh

# ADFRsuite installation script
pythonargs=" "
pyoptimize=0
TarDir=`pwd`
export ADS_ROOT=""
InstName=ADFRsuite
silent=0
usage() {
    echo "Usage: ./install.sh [-d InstDir] [-c optimization]"
    exit 
}
# Parse the command-line arguments
opts=`getopt "hlc:d:" "$@"`
if [ "$?" != 0 ]
then
    usage
fi

set -- $opts

while true; do
    case "$1" in 

    -c) pythonargs="$pythonargs -c"; pyoptimize="$2"; shift; shift ;;
    -d) export ADS_ROOT="$2"; shift ; shift ;;
    -l) pythonargs="$pythonargs -l"; silent=1; shift ;;
     
    -h) echo "Optional parameters:"
    echo "[-h]  help message;"
    echo "[ -d  InstDir] specifies installation directory (default-current directory)"
    echo "[ -c optimization] compile Python code with or without optimization:"
    echo "    0 - no optimization (generates .pyc files)"
    echo "    1 - with optimization (generates .pyo files);"
    exit ;;
    --) shift;;
    *)  if [ -z "$1" ] ; then break ; else echo "$1 is not a valid option" ; usage; fi ;; 
    esac
done


#echo "script options" python args "'$pythonargs'"  dest "'$ADS_ROOT'" pyoptimize $pyoptimize

if [ -e "./Version" ]; then
    vv=$(cat "./Version")
    InstName=${InstName}-$vv
    #echo "InstName" $InstName
fi
   
if [ "$ADS_ROOT" != "" ]; then
    # check if the user has write access to the installation directory
    if [ -e "$ADS_ROOT" ]; then
	if [ -d "$ADS_ROOT" ]; then
	    if [ ! -w  "$ADS_ROOT" ]; then 
		echo "Can not complete installation - specified directory $ADS_ROOT does not have write access."
		exit 1
	    fi
	    #check if the directory is empty
	    if [ "$(ls -A $ADS_ROOT)" ]; then
		# $ADS_ROOT directory is not empty, try to create ADFRsuite-#version ($InstName) to install our application in it:
		# 1 - check if the name of the $ADS_ROOT directory == $InstName (this is the case with the GUI installer)
		# 2 if not , check if the directory already contains $InstName
		#echo  "$(basename $ADS_ROOT)"
		createInstName=1
		if [ "$(basename $ADS_ROOT)" = $InstName ] ; then
		    createInstName=0
		    export ADS_ROOT="$ADS_ROOT"
		    #ADS_ROOT directory is the same as the $InstName
		    if [ $silent -eq 0 ]; then
			echo "$ADS_ROOT directory exists and is not empty."
			while true; do
			    read -p "Do you wish to overwrite it? (y/n)" yn
			    case $yn in
				[Yy]* ) rm -rf $ADS_ROOT/CCSBpckgs;  break;;
				[Nn]* ) exit;;
				* ) echo "Please answer y or n.";;
			    esac
			done
		    fi
		elif [ -e $ADS_ROOT/$InstName ]; then
		    createInstName=0
		    if [ $silent -eq 0 ]; then
			echo "$InstName exists in $ADS_ROOT directory."
			while true; do
			    read -p "Do you wish to overwrite it? (y/n)" yn
			    case $yn in
				[Yy]* ) rm -rf $ADS_ROOT/$InstName/CCSBpckgs; break;;
				[Nn]* ) exit;;
				* ) echo "Please answer y or n.";;
			    esac
			done
		    fi
		    export ADS_ROOT="$ADS_ROOT/$InstName"
		fi
		if [ $createInstName -eq 1 ] ; then
		    export ADS_ROOT="$ADS_ROOT/$InstName"
		    echo 
		    echo "Creating directory $ADS_ROOT"
		    mkdir -p "$ADS_ROOT"
		fi
	    fi
	else 
	    echo "$ADS_ROOT exists and is not a directory. "
	    exit 1
	fi
    else
	echo
	echo "Creating directory $ADS_ROOT"
	mkdir  -p "$ADS_ROOT"
    fi

else
    export ADS_ROOT="$(pwd)/$InstName"
    echo "Creating directoy $ADS_ROOT"
    mkdir  -p "$ADS_ROOT"
fi

echo
echo "Installing ADFRsuite to $ADS_ROOT;"
cd "$ADS_ROOT"
echo
echo "Installing Python Interpreter to $ADS_ROOT ... "
echo
tar xzf $TarDir/Python*.tar.gz

if [ "$?" != 0 ]; then
    echo "Error in Python installation"
    exit 1
fi
echo "Python installed, please wait for the rest of ADFRsuite to be installed."

cd $TarDir

## plaform we run on

export ADS_ARCHOSV=`$TarDir/Tools/archosv`

## add the path to the directory holding the python interpreter to your path

export PATH="$ADS_ROOT/bin:"$PATH

## use Python interpreter locally installed

PYTHON="$ADS_ROOT/bin/python2.7"
export PYTHONHOME="$ADS_ROOT"
if [ "`uname -s`" = "Linux" ] ; then
    export LD_LIBRARY_PATH="$ADS_ROOT/lib"
fi

## run python script - install.py - to install the packages and create scripts

if [ "$pyoptimize" -eq 1 ]; then
    echo "Running $PYTHON -O Tools/install.py $pythonargs"
    $PYTHON -O Tools/install.py $pythonargs
else
    echo "Running $PYTHON Tools/install.py $pythonargs"
    $PYTHON Tools/install.py  $pythonargs
fi

unset PYTHONHOME
