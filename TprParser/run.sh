#!/bin/bash

# clear
function clean
{
    rm -rf dist TprParser.egg-info 
}

# build
function build 
{
    python setup.py sdist 
}

# help
function showHelp
{
    echo "Usage:"
    echo $0 build
    echo $0 test
    echo $0 pypi
    echo $0 clean
}

# check if input
opt=($*)
N=${#opt[@]}
[[ $N == 0 ]] && showHelp && exit

inp=$1
# for testpypi
if [[ "$inp" == "test" ]]; then
    clean && build
    twine upload --repository testpypi dist/*
# for pypi
elif [[ "$inp" == "pypi" ]]; then
    clean && build
    twine upload --repository pypi dist/*
# only build
elif [[ "$inp" == "build" ]]; then
    clean && build
# only clean old file
elif [[ "$inp" == "clean" ]]; then
    clean
# show help
else
    showHelp
fi

