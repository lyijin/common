#!/usr/bin/env bash

# simple shell script to give python/shell/R scripts exec permissions.

chmod 755 ~/csiro/github/common/*.py
chmod 755 ~/csiro/github/common/*.sh
chmod 755 ~/csiro/github/common/*.R
chmod 755 *.py
chmod 755 *.sh
chmod 755 *.R

echo Granted 755 permissions to all \*.py, \*.sh and \*.R in ~/csiro/github/common!
echo Also, granted 755 permissions to *.py *.sh *.R in this directory!
