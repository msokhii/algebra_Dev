set -euo pipefail

g++ -std=c++17 main2.cpp helperF.cpp integerMath.cpp interpAlgo.cpp polyMath.cpp -I../header_FILES -O3 -march=native -flto -Ofast -fexceptions -fno-rtti -o time
