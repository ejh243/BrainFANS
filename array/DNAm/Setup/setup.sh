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

check_cmd() {
    command -v "$1" > /dev/null
}

get_conda_path() {
    default_conda_path="${HOME}/miniconda3"
    echo "Please enter the path you want to install conda to (default is $default_conda_path): "
    read -r conda_path
    if [[ -z "${conda_path}" ]]; then conda_path="${default_conda_path}"; fi
}

install_conda() {
    get_conda_path
    mkdir -p "${conda_path}"
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O "${conda_path}/miniconda.sh"
    bash "${conda_path}/miniconda.sh" -b -u -p "${conda_path}"
    rm "${conda_path}/miniconda.sh"

    source "${conda_path}/bin/activate"
    conda init --all
}

main() {
    if check_cmd "conda"; then
cat << MESSAGE
WARNING:
Conda is not installed/is not on PATH!
Would you like to install conda? (y/n)
MESSAGE
    read -r continue_install
    if [[ "${continue_install}" == "n" ]]; then exit 0; fi
        install_conda
    fi
    setup_conda_environment
    install_renv
    install_r_libraries
    check_installation
}

if [[ $# -ne 0 ]]; then usage; fi
main
