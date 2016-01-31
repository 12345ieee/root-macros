To compile a standalone with root libraries:  
`g++ -Wall $(root-config --cflags) -c "%f"`

To compile and link (generate) a standalone with root libraries:  
`g++ -Wextra -Og -g -Wall $(root-config --cflags --glibs) -o "%e" "%f"`

To add RooFit to the mix, add:  
`-L $ROOTSYS/lib -lRooFitCore -lRooFit` 

My ROOT distribution was configured using:  
`./configure linuxx8664gcc --prefix=/usr/local --all --enable-gsl-shared`

To run pyROOT (`import ROOT` hangs otherwise):  
`export PYTHONPATH="/usr/local/lib/root/"`
