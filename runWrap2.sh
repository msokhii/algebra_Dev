set -euo pipefail

g++ -O3 -std=c++17 -shared -fPIC -Iheader_FILES -o so_Files/fastVSolve.so pol_ALGO/vSolveWrapper.cpp pol_ALGO/fastVSolve.cpp pol_ALGO/polyMath.cpp pol_ALGO/integerMath.cpp pol_ALGO/helperF.cpp
