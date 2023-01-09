#!/bin/bash

echo -e "\033[0;32mGenerating PySCF API Docs using sphinx-apidoc\033[0m"

# Show the Python where PySCF is installed
echo "Using $(which python)"

# Get PySCF path
PYSCF_PATH=$(python -c "import pyscf; print(pyscf.__path__[0])")

if [ -z $PYSCF_PATH ]
then
    echo -e "\033[0;31mPySCF not found in current environment. Please switch environments or install it.\033[0m"
else
    echo -e "Found path to PySCF: \n\t${PYSCF_PATH}"
fi

# Run sphinx
DESTINATION=source/pyscf_api_docs
LOGFILE=_api_docs.log

mkdir -p $DESTINATION
echo -e "Output directory for API docs: \n\t$(pwd)/${DESTINATION}"
echo -e "Log file:\n\t$(pwd)/${LOGFILE}"

sphinx-apidoc -o $DESTINATION $PYSCF_PATH > ${LOGFILE}
