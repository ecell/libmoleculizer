#!/bin/bash

export LD_LIBRARY_PATH=/home/naddy/usr/lib

rm -f Moleculizer.so Moleculizer.o


cd ~/Sources/libmoleculizer
make -j 3
make install -j 3
cd src/python

# Compile the wrapper class
echo "Recompiling..."
g++ ReactionNetworkGenerator.cpp -g -I.. -I.. -I/usr/include/libxml++-2.6 -I/usr/lib/libxml++-2.6/include -I/usr/include/libxml2 -I/usr/include/glibmm-2.4 -I/usr/lib/glibmm-2.4/include -I/usr/include/sigc++-2.0 -I/usr/lib/sigc++-2.0/include -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include  -c
g++ -fPIC -c -g -o Moleculizer.o MoleculizerPythonWrapper.cpp -I/usr/include/python2.5  -I.. -I/usr/include/libxml++-2.6 -I/usr/lib/libxml++-2.6/include -I/usr/include/libxml2 -I/usr/include/glibmm-2.4 -I/usr/lib/glibmm-2.4/include -I/usr/include/sigc++-2.0 -I/usr/lib/sigc++-2.0/include -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include  
g++ -shared -o Moleculizer.so ReactionNetworkGenerator.o Moleculizer.o -L/home/naddy/usr/lib -lboost_python -lxml++-2.6 -lxml2 -lglibmm-2.4 -lgobject-2.0 -lsigc-2.0 -lglib-2.0  -lmoleculizer_mzr -lmoleculizer_cpx -lmoleculizer_ftr -lmoleculizer_dimer -lmoleculizer_mol -lmoleculizer_nmr -lmoleculizer_utl -lmoleculizer_stoch


cp Moleculizer.so /home/naddy/Sources/FirstCollisionSimulator/src/Simulator/rulesparsing/