To compile a standalone with root libraries:  
`g++ -Wextra -Og -g -Wall $(root-config --cflags --glibs) -o "%e" "%f"`

My ROOT distribution was configured using:  
`./configure linuxx8664gcc --prefix=/usr/local --all --enable-minuit2 --enable-gsl-shared`