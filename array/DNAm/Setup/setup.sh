#!/bin/bash

usage() {
cat << USAGE
============================================================================
$(basename "$0")
============================================================================
Purpose: Execute this script before running the QC pipeline. It will install
the necessary dependencies required automatically for you.
Contact: s.o.fletcher@exeter.ac.uk
============================================================================
USAGE
    exit 0
}

main() {
    if check_cmd "conda"; then
        install_conda
    fi
    setup_conda_environment
    install_renv
    install_r_libraries
    check_installation
}

if [[ $# -ne 0 ]]; then usage; fi
main
