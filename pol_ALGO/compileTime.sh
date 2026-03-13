set -euo pipefail

g++ main2.cpp fastVSolve.cpp helperF.cpp integerMath.cpp interpAlgo.cpp polyMath.cpp rfrWrapper.cpp -I../header_FILES -O3 -march=native -flto -DNDEBUG -o time
