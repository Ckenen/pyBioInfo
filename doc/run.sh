#!/bin/sh

sphinx-apidoc -o source ../src/pyBioInfo
sphinx-build -b html ./source ./build

