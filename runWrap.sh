set -euo pipefail 

g++ -O3 -shared -fPIC -Iheader_FILES -o pol_ALGO/wrapOBJ.so pol_ALGO/wrapperEx.cpp pol_ALGO/polyMath.cpp pol_ALGO/integerMath.cpp pol_ALGO/helperF.cpp
