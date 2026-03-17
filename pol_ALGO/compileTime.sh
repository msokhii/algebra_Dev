set -euo pipefail

g++ -std=c++17 main2.cpp helperF.cpp integerMath.cpp interpAlgo.cpp polyMath.cpp int128g.cpp -I../header_FILES -march=native -flto -O2 -finline-functions -funroll-loops -o time
