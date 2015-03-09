To compile a standalone wiyìth root libraries
g++ -Wextra -Og -g -Wall $(root-config --cflags --glibs) -o "%e" "%f"