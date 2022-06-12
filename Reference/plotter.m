function plotter()
N1 = 16;
N2 = 32;
N3 = 64;
N4 = 128;

LU_tridiagonal_solve(N1);
hold on;
LU_tridiagonal_solve(N2);
hold on;
LU_tridiagonal_solve(N3);
hold on;
LU_tridiagonal_solve(N4);
legend(strcat('N = ', num2str(N1)),strcat('N = ', num2str(N2)),strcat('N = ', num2str(N3)),strcat('N = ', num2str(N4)));
hold off;

end