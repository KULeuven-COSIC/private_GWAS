# private χ2 statistic
code accompanying Towards Practical Privacy-Preserving Genome-Wide Association Study

HE component

1. Getting started

Open the file HE_solution.cpp, go to the main part of the file at the bottom and change the necessary parameters (number of medical centers, maximum input value, parameters of the encoding, ...). Take care if you change the number of centers and the maximum input size, other parameters might need to be fine tuned. After saving the file you compile with g++ -std=c++11 -funroll-loops -Ofast -Wall -g -I nfllib/include HE_solution.cpp -o HE_solution -L nfllib/lib -lnfllib -lmpfr -lgmpxx -lgmp with the paths adapted to the correct folders. Then you can run the program with: ./HE_solution

2. Installation of the prerequisites

For ubuntu use the following commands:
apt-get install gcc
apt-get install libgmp-dev
apt-get install libmpfr-dev
apt-get install git
apt-get install cmake
git clone https://github.com/quarkslab/NFLlib.git
cd NFLlib/
mkdir _build
cd build/
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/nfllib
make
make test
make install
export LD_LIBARY_PATH=/../nfllib/lib/(path naar folder met libnfllib.so file)

3. Licensing

The license of the FV-NFLlib and the NFLlib is GPLv3.

MPC Component

1. Getting Started

Open the file MPC_GWAS.mpc, go to the main part of the file and change the necessary parameters (number of medical centers, threshold value for the χ2 significance, etc.). Save the file and compile it following the instructions detailed on the SPDZ-2 repository. 
Secret shared inputs can be introduced by following the instructions on the SPDZ-2 repository and Google group https://groups.google.com/forum/#!forum/spdz. 
Run the file and compile it following the instructions detailed on the SPDZ-2 repository. 

2. Prerequisites

Install SPDZ-2. Tested up to version 0.2. https://github.com/bristolcrypto/SPDZ-2. 
Configure SPDZ-2 to run in 256 bits inputs. 
Configure for the number of parties involved in the computation. 

3. Licensing

The license of the FV-NFLlib and the NFLlib is GPLv3.
