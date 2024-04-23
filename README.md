Repository to benchmark various integer factorisation algorithms

The project is designed to be build on Linux. Debian native installation and Ubuntu on WSL were tested

For building gcc, g++, and make are used


The project has several dependencies: GMP, Pari/GP, Nvidia CUDA, Python uses matplotlib, scipy, sagemath

-GMP

It is a preinstalled package on most distributions, including Ubuntu and Debian

If it is not present then it can be installed in the following way:

Download the archive from: https://gmplib.org/

Unpack:
tar -xvzf filename

Configure:
./configure

Build:
make

Run Unit Tests:
make test

Install:
make install

-Pari/GP

Archive with the code can be downloaded from 

Unpack:
tar -xvzf filename

Configure:
./Configure --with-gmp --with-readline

Build:
make install

-Nvidia CUDA

The drivers and devkit can be obtained from nvidia site: https://developer.nvidia.com/cuda-downloads
By following the commands provided in the end

In some cases there might be bath issues with nvcc not added to path:

then ~/.bashrc should be modified by appending following lines in the end:

export PATH=/usr/local/cuda/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/usr/local/cuda/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
export LD_LIBRARY_PATH=/usr/local/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}


Python dependencies can be installed with apt with the following commands:

sudo apt install python3-matplotlib
sudo apt install python3-scipy
sudo apt install python3-sage
