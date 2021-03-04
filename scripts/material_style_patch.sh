#!/bin/bash

# Make sure we're running this in the correct place
CWD=${PWD##*/}
if [ "$CWD" != "source" ]
then
    echo -e "\033[0;31mThis file should only be run from inside pyscf-doc/source\033[0m"
    exit
fi

# Apply/revert the patch
if [ "$1" == "apply" ] 
then
    echo "Applying patch and using material style for docs"
    cp conf.py conf.py.bkp
    cp index.rst index.rst.bkp
    cp ../patch_files/conf.py .
    cp ../patch_files/index.rst .
elif [ "$1" == "revert" ]
then
    echo "Reverting the patch to use material style for docs"
    cp conf.py ../patch_files/conf.py
    cp index.rst ../patch_files/index.rst
    mv conf.py.bkp conf.py
    mv index.rst.bkp index.rst
else
    echo "\033[0;31mIncorrect syntax! Try one of the following\033[0m"
    echo "1) material_style_path.sh apply"
    echo "2) material_style_path.sh revert"
fi

