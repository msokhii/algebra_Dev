set -euo pipefail

g++ -std=c++17 main2.cpp fastVSolve.cpp helperF.cpp integerMath.cpp interpAlgo.cpp polyMath.cpp rfrWrapper.cpp -I../header_FILES -O3 -march=native -flto -DNDEBUG -o time
