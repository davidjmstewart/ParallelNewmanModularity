# ParallelNewmanModularity

Development was done in WSL 2. VS Code launch configurations exist for the main application and all benchmark files. Running this (from the debug tab in VS Code) will rebuild and run the selected application. 

Running the `main` application will performance community detection on the Zachary Karate Network and output the communities to a text file. 

## Compiling

If you do not wish to use VS Code to build and run, the following compilation commands will be necessary:

Compile the libraries: `cd ./lib && gcc -o CDUtils.o -c CDUtils.c -fopenmp -O3 -march=native`

Compiling main: `gcc -fopenmp -g main.c -o main ./lib/CDUtils.o -lm -Wall -Werror -Wpedantic -Waggressive-loop-optimizations -O3 -march=native`

Run with `./main`

## Output

The program will place the communities file in `./benchmarking/matlab_testing_files/matlab-communities.txt`. The file should contain:

```
24 25 26 28 29 32 
9 10 15 16 19 21 23 27 30 31 33 34 
1 5 6 7 11 12 17 
2 3 4 8 13 14 18 20 22 
```

Each line is a community. 
