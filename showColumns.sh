#!/bin/bash

# aim: go from 'column -t $file | less -S' to simply 'showColumns.sh $file'

file=$1
column -t $file | less -S
