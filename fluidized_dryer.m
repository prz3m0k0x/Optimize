%% Cleaned & Structured MATLAB Script for Drying Process
% Author: University Assignment
% Purpose: Mass and energy balance of KCl drying with fluidized bed

clear; clc; close all;

%% ---------------- Parameters ----------------
% Dryer & Material
G_k = 12;             % ton/h, final product flowrate
mu_k = 0.004;         % final moisture content in KCl
mu_p = 0.05;          % initial moisture content
t_mat = 18;           % C, initial temperature KCl
t_1 = 650;            % C, inlet gas temperature
t_2 = 140;            % C, outlet gas temperature
q_loss_percent = 0.1; % % heat losses
rho_KCl = 1980;       % kg/m3
d_avg = 0.35e-3;      % m, average particle diameter
d_min = 0.15e-3;      % m, minimal particle size
h_bed = 0.4;          % m, bed height (fluidization stabilization zone)

% Gases
M_KCl = 39.10 + 35.45;
M_CO2 = 44.01; M_H2O = 18.02; M_O2 = 32.00; M_N2 = 28.01;
M_i = [M_CO2; M_H2O; M_O2; M_N2];
r_i = [0.05; 0.1; 0.13; 1 - 0.05 - 0.1 - 0.13];

R = 8.314; g = 9.81; % Gas constant & gravity
p = 1e5;              % Pa, process pressure

% Sutherland constants
T0 = 273.15; 
eta0 = [1.375;0.861;1.166;1.652]*1e-5;
C = [254;650;125;104];

%% ---------------- Mass Balance ----------------
G_p = G_k*(1-mu_k)/(1-mu_p)*1000/3600;  % kg/s, raw material + moisture
G_k = G_k*1000/3600;                     % kg/s, product
W = G_p - G_k;                           % kg/s, removed water
W_total = G_p*mu_p;                       % total initial water
W_remain = G_k*mu_k;                      % remaining water

%% ---------------- Energy Balance ----------------
% Heat capacity functions
fun_KCl = @(t) 35.42 + 70.03*t - 91.3823*t.^2 + 52.52*t.^3 + 0.1534./t.^2; 
DelH_KCl = integral(fun_KCl, 0.29115, 0.41315);
q_KCl = ((1-mu_p)*G_p*DelH_KCl*1000)/M_KCl; % kJ/s

fun_H2O_liquid = @(t) -203.6060 + 1523.290*t - 3196.413*t.^2 + 2474.455*t.^3 + 3.855326./t.^2;
DelH_H2O_liquid = integral(fun_H2O_liquid, 0.29115, 0.37315);
q_H2O_liquid = W_total*DelH_H2O_liquid*1000/M_H2O;

fun_H2O_steam = @(t) 30.092 + 6.832514*t + 6.793435*t.^2 - 2.53448*t.^3 + 0.082139./t.^2;
DelH_H2O_steam = integral(fun_H2O_steam, 0.37315, 0.41315);
q_H2O_steam = W*DelH_H2O_steam*1000/M_H2O;

DelH_H2O_vapor = 40.65; % kJ/mol, enthalpy of vaporization
q_H2O_vapor = W*DelH_H2O_vapor*1000/M_H2O;

q_H2O_total = q_H2O_vapor + q_H2O_steam + q_H2O_liquid; 
q_loss = q_loss_percent*(q_H2O_total + q_KCl); 

q_gas_total = q_H2O_total + q_KCl + q_loss;  % kJ/s
q_per_kg_water = q_gas_total / W;            % kJ/kg water

%% ---------------- Gas Properties ----------------
t_mean = (t_1 + t_2)/2; % mean gas temperature

% Heat capacities at 300 and 400 C for interpolation
T_points = [300;400];
Cp_CO2 = [0.949;0.983]; Cp_H2O = [1.919;1.948];
Cp_O2 = [0.95;0.965]; Cp_N2 = [1.057;1.066];
Cp_gases = [Cp_CO2 Cp_H2O Cp_O2 Cp_N2];
Cp_interpolated = interp1(T_points,Cp_gases,t_mean)'; % transpose to column


% Molar mass of gas mixture
M_z = sum(M_i.*r_i);
g_i = r_i.*M_i/M_z;  % mass fraction

Cp_gas_mix = sum(g_i .* Cp_interpolated);

% Gas flowrate based on heat balance
G_s = q_gas_total/(Cp_gas_mix*(t_1-t_2)); % kg/s

% Gas densities at inlet and outlet
rho_gas = M_z*p ./ ([t_2+T0; t_1+T0]*R*1000); % kg/m3

%% ---------------- Fluidization ----------------
% Gas viscosity
eta_i = eta0 .* (T0 + C)./(t_2 + T0 + C) .* ((t_2 + T0)/T0).^(3/2);
ln_eta_g = sum(log(eta_i).*g_i);
eta_g = exp(ln_eta_g);

% Archimedes number & fluidization velocity
Ar = (d_avg^3*rho_KCl*g)/eta_g^2;
Ly_kr = 1e-3; Ly_r = 1.5;
w_kr = (Ly_kr*eta_g*g*rho_KCl)^(1/3); 
K_r = (Ly_r/Ly_kr)^(1/3);
w = K_r*w_kr; % m/s, fluidization velocity
w_ds = w*(T0+t_1)/(T0+t_2); 
w_otw = w_ds/0.1;

% Distributor area
F_ds = G_s/(rho_gas(2)*w_ds);
D_ds = sqrt(4*F_ds/pi);

%% ---------------- Outlet & Separator ----------------
h_sep = 4*h_bed;
H_1 = h_sep + h_bed;
H_2 = 0.75*H_1;

G_o = G_s + G_p - G_k; % outlet gas mass flow
r_w = [0;1;0;0];       % water entrained
g_o_i = (G_s*g_i + r_w*W)/G_o;
M_o_z = sum(M_i.*g_o_i);
rho_o = M_o_z*p/(R*(t_2+T0)*1000);
w_rz = G_o/(rho_o*F_ds);

% Adjust separator area to avoid entrainment
F_sep = F_ds;
while w_rz > w
    F_sep = 1.05*F_sep;
    w_rz = G_o/(rho_o*F_sep);
end
D_wylot = sqrt(4*F_sep/pi);

%% ---------------- Plots ----------------
% Plot 1: Heat requirement vs water removed
figure; 
bar([q_H2O_liquid q_H2O_steam q_H2O_vapor]/1000);
set(gca,'XTickLabel',{'Heating','Steam','Vaporization'});
ylabel('Heat [kJ/s]');
title('Heat Contribution to Water Removal');

% Plot 2: Gas velocity profile
figure;
plot([t_1 t_2], [w_ds w_otw],'o-','LineWidth',2);
xlabel('Gas Temperature [°C]');
ylabel('Gas Velocity [m/s]');
title('Gas Velocity in Fluidized Bed');

% Plot 3: Gas flowrate vs Cp
figure;
bar(Cp_interpolated);
set(gca,'XTickLabel',{'CO2','H2O','O2','N2'});
ylabel('Cp [kJ/kg*K]');
title('Specific Heat of Gas Components');

%% ---------------- Display Key Results ----------------
fprintf('\n--- Key Results ---\n');
fprintf('Removed water W = %.4f kg/s\n', W);
fprintf('Heat per kg water q = %.2f kJ/kg\n', q_per_kg_water);
fprintf('Gas flowrate G_s = %.4f kg/s\n', G_s);
fprintf('Distributor diameter D_ds = %.4f m\n', D_ds);
fprintf('Separator diameter D_wylot = %.4f m\n', D_wylot);
