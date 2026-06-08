gcc main.c tridiagonal_system.c helper.c ../utils.c -o broyden -O3 -march=native -mavx -fopt-info-vec -Wall -Wextra -Wpedantic -Werror -flto
