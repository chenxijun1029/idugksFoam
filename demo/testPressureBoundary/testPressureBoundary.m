%test pressure boundary condition

%parameter
L = 1.0;
H = 1.0; 
vis  = 0.001
P_in  = 0.33334
P_out = 1/3
G = (P_in - P_out) / L;

%analytical solution
y_ect = (0:0.001:1);
U_ect = 0.5 * G * (y_ect .* H - y_ect .* y_ect) / vis;
%idugksFoam solution
row = 1;
col_Ux = 0;
col_Py = 7;
y = csvread('U_center.csv', row, col_Py, [row col_Py 101 col_Py]);
U = csvread('U_center.csv', row, col_Ux, [row col_Ux 101 col_Ux]);

%plot
figure('Name','testPressureBoundary');hold on;
plot(y_ect, U_ect, 'r');
plot(y, U, 'ro');
xlabel('y/H');
ylabel('Ux');
legend('analytical solution','idugksFoam solution')