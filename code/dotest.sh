gcc -g -Wall -c -lm permutationreduction.c
gcc -g -lm meschach/*.o permutationreduction.o randn.o qrmcp.o reduction.o search.o test2.c -o test2.o
./test2.o
