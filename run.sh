set -euo pipefail

g++ -I . -Iheader_FILES curr_BUILD/main.cpp pol_ALGO/integerMath.cpp pol_ALGO/polyMath.cpp pol_ALGO/helperF.cpp -o output_FILES/ratTEST
