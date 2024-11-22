# Matrix Chain Essential

This project is a streamlined codebase that supports the generation and computation of the cost of all orders of computation for a matrix chain of any size.

The project includes the implementation of approximation algorithms for the matrix chain: Chandra's (A.K. Chandra, 1975), Chin's (F.Y. Chin, 1978), and two novel algorithms that always yield a better solution than the former.

## Requirements

* A compiler that supports C++17.
* CMake 3.15 or higher.

## Setup and compilation

Clone the repository using:

```bash
git clone git@github.com:HPAC/essential.git
```

The directory `src` contains the implementation of the library that generates and computes the cost of the different orders of computation. There are some other functionalities included, such as the generation of random instances and implementations of various approximation algorithms for the matrix chain problem. 

The directory `test` contains executables that leverage the library in `src`. The files in `test` are effectively experiments performed on the set of parenthesisations and/or the approximation algorithms for the matrix chain.

Use the following commands to compile it:

```bash
cd essential
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -S .. -B .
cmake --build .
```

## Executing experiments

After compiling, the directory `build/test` will contain some executables. These are and can be used as:

* `build/test/experiment` takes two arguments: 1) the length of the chain; 2) the number of instances to test upon. The program returns metrics (max, avg, freq-penalty) for each approximation algorithm on the set of randomly generated instances (whose sizes are in the range 1-1000). Example: `./experiment 7 100000`, where `7` is the length of the chain and `100000` is the number of instances to generate.

* `build/test/max_pen` takes two arguments: 1) the length of the chain; 2) the number of instances to test upon. The program prints metrics for Chin's algorithm and our improved version. It also prints the instance for which maximum penalty was found for each approximation algorithm.

* `build/test/single_instance` takes as arguments: 1) the length of the chain; and as many dimensions as needed to specify an instance for a chain of the given length. Example `./single_instance 5 100 77 94 42 212 44`, where `5` is the length of the chain and the other six numbers specify the input instance. The program returns the penalty for each approximation algorithm on the input instance.
