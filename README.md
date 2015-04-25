To compile a standalone with root libraries:  
`g++ -Wextra -Og -g -Wall $(root-config --cflags --glibs) -o "%e" "%f"`