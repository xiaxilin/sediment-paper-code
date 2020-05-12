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
g = 9.81;
n = 0.035;
mu = 0.5;
A_g = 0.005;


f_1 = g*S_b;
f_2 = A_g;
syms x;
eqn = f_1 - rho_w/(rho_w*h_0 + (rho_s-rho_w)*f_2*x^2)*g*n^2*(h_0 - f_2*x^2/(1-p))^(-1/3)*x^2 ...
            - (rho_s-rho_w)/(rho_w*h_0 + (rho_s-rho_w)*f_2*x^2)*g*f_2*x^2*mu;
u_0 = vpasolve(eqn,x);
u_0 = abs(double(u_0));
r = f_1 - rho_w/(rho_w*h_0 + (rho_s-rho_w)*f_2*u_0^2)*g*n^2*(h_0 - f_2*u_0^2/(1-p))^(-1/3)*u_0^2 ...
            - (rho_s-rho_w)/(rho_w*h_0 + (rho_s-rho_w)*f_2*u_0^2)*g*f_2*u_0^2*mu;


du = 0.1*u_0;
hu_0 = h_0*(u_0+du);
q_b_star = A_g*u_0^3;
hC_0 = q_b_star/u_0;
C_0 = hC_0/h_0;

h_f = h_0 - hC_0/(1-p);
rho = rho_w + (rho_s - rho_w)*C_0;
r_1 = g*S_b - g*rho_w/rho*n^2*h_f^(-1/3)/h_0*u_0^2 - g*(rho_s - rho_w)/rho/h_0*hC_0*mu;

%
hC(1) = hC_0;
hu(1) = hu_0;
h(1) = 1.0*h_0;
C(1) = hC(1)/h(1);
u(1) = hu(1)/h(1);

%%
N = 6000;
for i = 1:N
    % first step of RK2
    q_b =  hu(i)*hC(i)/h(i);
    u = hu(i)/h(i);
    C = hC(i)/h(i);
    rho = rho_s*C + rho_w*(1-C);
    q_b_star = A_g*u^3;
    dhdt = -1.0*(q_b - q_b_star)/L/(1-p);
    dhcdt = -1.0*(q_b - q_b_star)/L;
    h_f = h(i) - hC(i)/(1-p);
    S_f = rho_w/rho*n^2*h_f^(-1/3)/h(i)*u^2 + (rho_s - rho_w)/rho/h(i)*hC(i)*mu;
    dhudt = g*h(i)*(S_b-S_f) - (rho_s - rho_w)/rho*u*dhdt*(1-p-C);
    hu_1 = hu(i) + dt*dhudt;
    h_1 = h(i) + dt*dhdt;
    hC_1 = hC(i) + dt*dhcdt;
    % Second step of RK2
    q_b =  hu_1*hC_1/h_1;
    u = hu_1/h_1;
    C = hC_1/h_1;
    rho = rho_s*C + rho_w*(1-C);
    q_b_star = A_g*u^3;
    dhdt_1 = -1.0*(q_b - q_b_star)/L/(1-p);
    dhcdt_1 = -1.0*(q_b - q_b_star)/L;
    h_f = h_1 - hC_1/(1-p);
    S_f = rho_w/rho*n^2*h_f^(-1/3)/h_1*u^2 + (rho_s - rho_w)/rho/h_1*hC_1*mu;
    dhudt_1 = g*h_1*(S_b-S_f) - (rho_s - rho_w)/rho*u*dhdt_1*(1-p-C);
    hu(i+1) = hu(i) + 0.5*dt*(dhudt + dhudt_1);
    h(i+1) = h(i) + 0.5*dt*(dhdt + dhdt_1);
    hC(i+1) = hC(i) + 0.5*dt*(dhcdt + dhcdt_1);
end
%%
figure;
plot((1:N)*dt, h(1:N)/h_0,'k-','linewidth',1.0);
ylabel('$h/h_0$','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');
title('(a) $h$','Interpreter','latex');
axis([0 600 0.99 1.05]);
set(gca,'fontsize',18);
%%
figure
plot((1:N)*dt, hC(1:N)./h(1:N)/C_0,'k-','linewidth',1.0);
ylabel('$C/C_\infty$','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');
title('(b) $C$','Interpreter','latex');
axis([0 600 0.99 1.1]);
set(gca,'fontsize',18);
%%
figure
plot((1:N)*dt, hu(1:N)./h(1:N)/u_0,'k-','linewidth',1.0);
ylabel('$u/u_\infty$','Interpreter','latex');
xlabel('$t$ (s)','Interpreter','latex');
title('(c) $u$','Interpreter','latex');
axis([0 600 0.99 1.11]);
set(gca,'fontsize',18);