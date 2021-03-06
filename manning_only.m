%% A matlab script for numerically solving the sediment transport of uniform flow
clear,clc

%% parameters
h_0 = 0.5;
S_b = 0.15;
p = 0.4;
L = 1.0;
dt = 0.1;
d_fifty = 0.001;
rho_s = 2600;
rho_w = 1000;
rho = rho_s*(1-p) + rho_w*p;
g = 9.81;
n = 0.035;
miu = 0.5;
A_g = 0.005;

%
tau = rho*g*h_0*S_b;
u_0 = sqrt(h_0*S_b/(n^2*h_0^(-1/3)));
du = 0.1*u_0;
hu_0 = h_0*(u_0+du);
q_b_star = A_g*u_0^3;
hC_0 = q_b_star/u_0;
C_0 = hC_0/h_0;
Froud = u_0/sqrt(g*h_0);

%
hC(1) = hC_0;
hu(1) = hu_0;
h(1) = 1.0*h_0;
C(1) = hC(1)/h(1);
u(1) = hu(1)/h(1);
%%
N = 6000;

%

for i = 1:N
% Fully implicit scheme for sediment transport
   q_b =  hu(i)*hC(i)/h(i);
   u = hu(i)/h(i);
   C = hC(i)/h(i);
   U_0 = [hC(i) h(i) u]';
   U_1 = [hC(i) h(i) u]';
   for j = 1:10
    J = myfunc_Jac(A_g, L, g, n, p, U_1(1), U_1(2), U_1(3), rho_s, rho_w);
    q_b_star = A_g*U_1(3)^3;
    S_1 = -(U_1(3)*U_1(1)-q_b_star)/L;
    S_2 = -(U_1(3)*U_1(1)-q_b_star)/L/(1-p);
    S_3 = -S_2*U_1(3)/U_1(2) + g*S_b - g*n^2*U_1(2)^(-4/3)*U_1(3)^2 -(rho_s - rho_w)/(rho_s - (rho_s -rho_w)*U_1(1)/U_1(2))/U_1(2)*U_1(3)*S_2*(1-p-U_1(1)/U_1(2));
    S = [S_1 S_2 S_3]';
    I = eye(3,3);
    U_2 = U_1 + (I-dt*J)^(-1)*(dt*S-(U_1-U_0));
    U_1 = U_2;
   end
   hC(i+1) = U_1(1);
   h(i+1) = U_1(2);
   hu(i+1) = U_1(2)*U_1(3);
end
%%
figure;
plot((1:N)*dt, h(1:N)/h_0,'k-','linewidth',1.0);
ylabel('$h/h_0$','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');
title('(a) $h$','Interpreter','latex');
set(gca,'fontsize',18);
%%
figure
plot((1:N)*dt, hC(1:N)./h(1:N)/(1-p),'k-','linewidth',1.0);
ylabel('$C/(1-p)$','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');
title('(b) $C$','Interpreter','latex');
set(gca,'fontsize',18);
%%
figure
plot((1:N)*dt, hu(1:N)./h(1:N)/u_0,'k-','linewidth',1.0);
ylabel('$u/u_0$','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');
title('(c) $u$','Interpreter','latex');
set(gca,'fontsize',18);
