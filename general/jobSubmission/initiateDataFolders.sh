#!/bin/bash
# create default data directories in given location

set -e

# one argument: directories location

cd "$1"

mkdir 0_metadata 1_raw 2_trimmed 3_aligned 4_catMemes


