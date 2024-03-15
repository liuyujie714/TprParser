#!/bin/bash

rm -rf dist TprParser.egg-info && python setup.py sdist 
twine upload dist/*
