set -euo pipefail 

g++ -O3 -std=c++17 -shared -fPIC -Iheader_FILES -o pol_ALGO/rfr.so pol_ALGO/rfrWrapper.cpp pol_ALGO/polyMath.cpp pol_ALGO/integerMath.cpp pol_ALGO/helperF.cpp
