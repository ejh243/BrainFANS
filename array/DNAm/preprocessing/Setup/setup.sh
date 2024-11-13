#!/bin/bash

usage() {
cat << USAGE
============================================================================
$(basename "$0")
============================================================================
Purpose: Execute this script before running the QC pipeline. It will install
the necessary dependencies required automatically for you.
Arguments: \$1 -> full/relative path to config.txt file
Contact: s.o.fletcher@exeter.ac.uk
============================================================================
USAGE
    exit 0
}

source_config_file() {
    config_file_path=$1
    source "${config_file_path}" || {
            echo "${RED}ERROR: config file not found at ${config_file_path}. Please check this path is correct.${NO_COLOUR}"
            exit 1
        }
}

check_cmd() {
    command -v "$1" > /dev/null
}

print_conda_missing_message() {
cat << MESSAGE
${RED}
WARNING:
Conda is not installed/is not on PATH!
Would you like to install conda? (y/n)
${NO_COLOUR}
MESSAGE
}

get_conda_path() {
    default_conda_path="${HOME}/miniconda3"
    echo "Please enter the path you want to install conda to (default is ${default_conda_path}, press enter for this path): "
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

get_conda_shell() {
cat << EOF
Could not find conda.sh at ${conda_shell_location}. Please provide the full
path to this file. If you are finding this difficult, it should be inside your
conda installation directory under:
.../conda-installation-folder/etc/profile.d/conda.sh
EOF
    read -r conda_shell_location
    conda_path=${conda_shell_location%/etc/profile.d/conda.sh}
}

find_conda_shell() {
    # Initial guess will work if this script installed conda for the user
    conda_shell_location="${conda_path}/etc/profile.d/conda.sh"
    while [[ ! -f "${conda_shell_location}" ]]; do
        get_conda_shell
    done
}

setup_conda_environment() {
    environment_name="BrainFANS_DNAm"
    if ! conda env list | grep -q "${environment_name}"; then
        conda config --add channels bioconda
        conda config --add channels conda-forge
        conda create \
            -y \
            --name "${environment_name}" \
            --file "$(dirname "$0")/requirements-${environment_name}.txt"
    fi
}

install_r_libraries() {
    # We need to ensure that renv sees the packages installed in the conda
    # environment to avoid certain compilation errors with system based 
    # libraries (e.g. systemfonts).
cat > "${SCRIPTSDIR}/array/DNAm/preprocessing/.Rprofile" << EOF
source("renv/activate.R")
.libPaths(c("$conda_path/envs/$environment_name/lib/R/library", .libPaths()))
EOF
    echo "Installing R libraries using renv, please follow on-screen instructions."
    cd "${SCRIPTSDIR}/array/DNAm/preprocessing/" || exit 1
    Rscript "installPackages.R"
}

print_installation_unsuccessful_message() {
cat << MESSAGE
${RED}
ERROR:
Installation of R libraries was not successful. 
Please create a bug report at: 
${BLUE}https://github.com/ejh243/BrainFANS/issues/new/choose
${NO_COLOUR}
MESSAGE
}

check_installation() {
    if Rscript "${SCRIPTSDIR}/array/DNAm/preprocessing/checkRPackages.R"; then
        echo "${GREEN}Installation successful${NO_COLOUR}"
    else
        print_installation_unsuccessful_message
        exit 1
    fi
}

add_to_config_file() {
    config_file_path=$1
cat >> "${config_file_path}" << EOF 

CONDA_SHELL=${conda_shell_location}
DNAM_CONDA_ENVIRONMENT=${environment_name}

EOF
}

main() {
    config_file_path=$1
    source_config_file "$config_file_path"

    if ! check_cmd "conda"; then
        print_conda_missing_message
        read -r install
        if [[ "${install}" == "n" ]]; then exit 1; fi
            install_conda
    else
        conda_path=${CONDA_EXE%/bin/conda}
    fi
    find_conda_shell
    source "${conda_shell_location}" || exit 1
    setup_conda_environment
    conda activate "${environment_name}"
    
    # This assists R in finding the correct pkg-config files so that curl
    # can be properly configured.
    export PKG_CONFIG_PATH="${conda_path}/envs/${environment_name}/lib/pkgconfig"
    install_r_libraries
    check_installation
    add_to_config_file "${config_file_path}"
}

if [[ $# -ne 1 ]]; then usage; fi
RED='[0;31m'
BLUE='[0;34m'
GREEN='0;32m'
NO_COLOUR='[0m'
main "$1"
