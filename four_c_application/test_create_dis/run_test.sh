#!/bin/sh
export LD_LIBRARY_PATH=/opt/4C/lib:$LD_LIBRARY_PATH
python pytest_create_dis.py

