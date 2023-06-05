#!/bin/bash

# Installation script
# Code to install all the required packages to run the workflow

# Usage ./installation_code.sh path/to/folder

# Almost all required packages will be installed using conda within snakemake
# You just need to add the flag '--use-conda' to the snakemake command

# Install packages not present in the conda environment
WD=${PWD}

# Create the binaries folder
mkdir -p bin

# Create program packages folder
mkdir -p programs
cd programs

# FAMSA (multiple sequence alignments)
wget https://github.com/refresh-bio/FAMSA/releases/download/v1.2.1/famsa-1.2.1-linux -O "${WD}"/bin/famsa
chmod +x "${WD}"/bin/famsa
###############################################################################

# OD-seq (Outlier detection in multiple sequence alignments)
# In the mac systems we need to install gcc with brew (brew install gcc), then call specifically that gcc version (ex: /usr/local/Cellar/gcc/8.3.0/bin/g++-8)
wget http://www.bioinf.ucd.ie/download/od-seq.tar.gz
tar -zxf od-seq.tar.gz
cd OD-Seq/
g++ -fopenmp -o "${WD}"/bin/OD-seq AliReader.cpp Bootstrap.cpp DistCalc.cpp DistMatReader.cpp DistWriter.cpp FastaWriter.cpp IQR.cpp ODseq.cpp PairwiseAl.cpp Protein.cpp ResultWriter.cpp runtimeargs.cpp util.cpp
cd ..
###############################################################################

# HHblits
git clone https://github.com/soedinglab/hh-suite
cd hh-suite
git submodule init
git submodule update
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX="${WD}"/bin/hh-suite ..
make -j 8
make install
cd ../../

###############################################################################

# MMseqs2
git clone https://github.com/soedinglab/MMseqs2.git
cd MMseqs2
git checkout 9cc89aa594131293b8bc2e7a121e2ed412f0b931
mkdir build
cd build
cmake -DHAVE_MPI=1 -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX="${WD}" ..
make -j 8
make install
cd ../..
###############################################################################

# Kaiju
git clone https://github.com/bioinformatics-centre/kaiju.git
cd kaiju/src
make -j 8
cd ../
cp bin/* "${WD}"/bin/
cd ../
###############################################################################

# HMMER MPI-mode
wget http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz
tar zxf hmmer-3.3.tar.gz
cd hmmer-3.3
./configure --prefix="${WD}" --enable-mpi
make -j 8
make check
make install
cd ..
####################################################################################

# Eggnog-mapper
git clone https://github.com/eggnogdb/eggnog-mapper.git
cd eggnog-mapper
python setup.py install
pip install biopython
python download_eggnog_data.py -y
cd ../
###############################################################################

# Parasail
git clone https://github.com/jeffdaily/parasail
cd parasail
sed -i "s|int matches = parasail_result_get_matches(result);|int matches = parasail_result_get_similar(result);|" apps/parasail_aligner.cpp
autoreconf -fi
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX="${WD}" ..
make -j 8
make install
cd ../..
export LD_LIBRARY_PATH=/vol/cloud/agnostos-wf/lib:${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
###############################################################################

# Igraph C-library
wget https://github.com/igraph/igraph/releases/download/0.10.4/igraph-0.10.4.tar.gz
tar xvfz igraph-0.10.4.tar.gz
cd igraph-0.10.4
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX="${WD}"/bin/igraph ..
cmake --build .
#cmake --build . --target check
cmake --install .
cd ../../../

export LD_LIBRARY_PATH=/vol/cloud/agnostos-wf/bin/igraph/lib:${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

gcc workflow/scripts/is_connected.c -o workflow/scripts/is_connected -Ibin/igraph/include -Lbin/igraph/lib -ligraph -lm -lstdc++ -lgomp -lpthread
gcc workflow/scripts/filter_graph.c -o workflow/scripts/filter_graph -Ibin/igraph/include -Lbin/igraph/lib -ligraph -lm -lstdc++ -lgomp -lpthread
