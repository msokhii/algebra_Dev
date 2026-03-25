set -euo pipefail 

g++ -O3 -std=c++17 -shared -fPIC -Iheader_FILES -o pol_ALGO/newton.so pol_ALGO/newtonWrapper.cpp pol_ALGO/polyMath.cpp pol_ALGO/integerMath.cpp pol_ALGO/helperF.cpp pol_ALGO/int128g.cpp pol_ALGO/interpAlgo.cpp
