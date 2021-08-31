#!/usr/bin/env bash
#
# Copyright 2020 Axel Huebl
#
# License: BSD-3-Clause-LBNL

# search recursive inside a folder if a file contains tabs
#
# @result 0 if no files are found, else 1
#

set -eu -o pipefail

sudo apt-get update

sudo apt-get install -y --no-install-recommends\
    build-essential     \
    ca-certificates     \
    gnupg               \
    libopenmpi-dev      \
    openmpi-bin         \
    pkg-config          \
    wget

sudo wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/7fa2af80.pub
sudo apt-key add 7fa2af80.pub
echo "deb https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64 /" | sudo tee /etc/apt/sources.list.d/cuda.list
sudo apt-get update
CUDA_DASH_VERSION=11-2
CUDA_DOT_VERSION=11.2
sudo apt-get install -y cuda-command-line-tools-${CUDA_DASH_VERSION} cuda-compiler-${CUDA_DASH_VERSION} cuda-cupti-dev-${CUDA_DASH_VERSION} cuda-minimal-build-${CUDA_DASH_VERSION} cuda-nvml-dev-${CUDA_DASH_VERSION} cuda-nvtx-${CUDA_DASH_VERSION} libcurand-dev-${CUDA_DASH_VERSION} libcusparse-dev-${CUDA_DASH_VERSION} libcusolver-dev-${CUDA_DASH_VERSION} libcublas-dev-${CUDA_DASH_VERSION}
sudo ln -s cuda-${CUDA_DOT_VERSION} /usr/local/cuda
