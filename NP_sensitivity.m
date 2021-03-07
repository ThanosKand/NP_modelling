% Data
L = 200; % [m]
param.D = 5*8.64; % [m^2/day]
param.u = 0.04*24; % [m/day]

param.t_range = 0:1000;
nopoints = 200;
param.dx = L / nopoints; %[m]
param.z = 0.5*param.dx:param.dx:(L-0.5*param.dx); 
param.c0 = [ones(1, nopoints)*100 ones(1, nopoints)*1e06]; % [mmol nutrient/m3] [cells/ m3] 

param.H_I = 20*86400; % [μmol photons m-2 day-1] half-saturation constant of light-limitead growth 
param.H_N = 0.02; % [mmol nutrient/m3]

param.I0 = 450*86400; %incident light intensity [μmol photons m-2 day-1]

param.kw = 0.045; % background turbidity [m-1]
param.kp =6e-12; % specific light attenuation of phytoplankton [m2/ cell]

param.mu_max = 0.04*24; %[1/day]
param.loss = 0.01*24; %[1/day]
param.alpha = 1*10^-9; % Nutrient content of phytoplankton [mmol nutrient/cell]
param.eps = 0.5; % Nutrient recycling coefficient []

%%
% [t, C] = NPZD(param);
% 
n = length(param.z);
% N = C(:, 1:n);
% P = C(:, n+1:2*n);

%%


u_new = [0, param.u, param.u*8];
% u_new = param.u*10;


maxvalue = zeros;
index = zeros;
I_light = zeros(length(u_new), n);
nutr_lm = zeros(length(u_new), n);
light_lm = zeros(length(u_new), n);
PP_last = zeros(length(u_new), n);


for i = 1:length(u_new)

param.u = u_new(i);

[t, C] = NPZD(param);
N = C(:, 1:n);
P = C(:, n+1:2*n);
% ADD function to plot stuff for each different value of u!!!

plotting_2(N,P,param,t)


P_last = P(end,:); 
PP_last(i,:) = P_last;
N_last = N(end,:);

[maxvalue(i), index(i)] = max(P_last);

I_light(i,:) = calclight(param.z,t, P_last, param.dx, param.kp, param.kw, param.I0);
nutr_lm(i,:) = (N_last./ (param.H_N + N_last));
light_lm(i,:) = (I_light(i,:) ./ (param.H_I + I_light(i,:)));
end

c = param.z(index);

%%
figure;

subplot(1,2,1)
plot(u_new, maxvalue,'-o', 'Linewidth', 1.5)
xlabel('Settling velocity [m/day]')
ylabel('Max [P] in steady state [cells/ m3] ')

subplot(1,2,2)
plot(u_new, c,'-ro', 'Linewidth', 1.5)
xlabel('Settling velocity [m/day]')
ylabel('Depth of max [P] in steady state [m]')
axis ij

%%

figure;

subplot(1,3,1)
line(PP_last(1,:), param.z, 'Linewidth', 1.5) % Half time
ax1 = gca; % current axes
axis ij
ylabel( 'Depth [m]')
xlabel('Phytoplankton conc. [cells/ m3]')
%ax1.XColor = '';
%ax1.YColor = 'b';
% ax1.XAxisLocation = 'top';

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top', 'YAxisLocation','right', 'Color','none');
line(nutr_lm(1,:), param.z, 'Parent', ax2, 'Linewidth', 1.2, 'Color','r', 'LineStyle', '--')
axis ij
line(light_lm(1,:), param.z, 'Parent', ax2, 'Linewidth', 1.2, 'Color','k', 'LineStyle', '-.')
axis ij
xlabel('Limitation by light or nutrients')

% legend('Nutrients','Light', 'Location', 'southeast')
title('v = 0 m/day')

subplot(1,3,2)
line(PP_last(2,:), param.z, 'Linewidth', 1.5) % Half time
ax1 = gca; % current axes
axis ij
ylabel( 'Depth [m]')
xlabel('Phytoplankton conc. [cells/ m3]')
%ax1.XColor = '';
%ax1.YColor = 'b';
% ax1.XAxisLocation = 'top';

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top', 'YAxisLocation','right', 'Color','none');
line(nutr_lm(2,:), param.z, 'Parent', ax2, 'Linewidth', 1.2, 'Color','r', 'LineStyle', '--')
axis ij
line(light_lm(2,:), param.z, 'Parent', ax2, 'Linewidth', 1.2, 'Color','k', 'LineStyle', '-.')
axis ij
xlabel('Limitation by light or nutrients')

% legend('Nutrients','Light', 'Location', 'southeast')
title('v = 0.96 m/day')

subplot(1,3,3)
line(PP_last(3,:), param.z, 'Linewidth', 1.5) % Half time
ax1 = gca; % current axes
axis ij
ylabel( 'Depth [m]')
xlabel('Phytoplankton conc. [cells/ m3]')
%ax1.XColor = '';
%ax1.YColor = 'b';
% ax1.XAxisLocation = 'top';

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top', 'YAxisLocation','right', 'Color','none');
line(nutr_lm(3,:), param.z, 'Parent', ax2, 'Linewidth', 1.2, 'Color','r', 'LineStyle', '--')
axis ij
line(light_lm(3,:), param.z, 'Parent', ax2, 'Linewidth', 1.2, 'Color','k', 'LineStyle', '-.')
axis ij
xlabel('Limitation by light or nutrients')

legend('Growth limiting factor with respect to nutrients','Growth limiting factor with respect to light', 'Location', 'southeast')
title('u = 7.7 m/day')


%%
figure
plot(I_light, param.z)
axis ij



%%
function [t, C] = NPZD(param)

n = length(param.z);

% opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
options1 = odeset('Refine',4);
opts = odeset(options1, 'NonNegative',4);

[t,C] = ode45(@derivative, param.t_range, param.c0, opts);

function dCdt = derivative(t,C)

N = C(1:n);
P = C(n+1:2*n);
    
ix = 2:n;
%--- Nutrient fluxes ---
% Diffusive fluxes:
Jdiff_N(ix) = -param.D*(N(ix)-N(ix-1))/param.dx;
Jdiff_N(1) = 0; % No flux at the surface...
Jdiff_N(n+1) = -param.D * (100 - N(end))/ param.dx;  % ...or the bottom

J_N = Jdiff_N;

%--- Phytoplankton fluxes ---
% Advective fluxes
Jadv_P(ix) = param.u*P(ix-1);
Jadv_P(1) = 0;  % No input from the surface
Jadv_P(n+1) = 0; % Closed bottom

% Diffusive fluxes:
Jdiff_P(ix) = -param.D*(P(ix)-P(ix-1))/param.dx;
Jdiff_P(1) = 0; % No flux at the surface...
Jdiff_P(n+1) = 0;  % ...or the bottom

% Rate-of-change due to advection and diffusion:
J_P = Jadv_P + Jdiff_P;    

% Light levels for each cell in each timestep
I = calclight(param.z,t,P,param.dx,param.kp,param.kw,param.I0);

% Specific growth rate in each cell as a function of light
ix_f = 1:n;

% mu = zeros;
% for i = 1:n
%     mu(i) = param.mu_max * min((N(i)/ (param.H_N + N(i))), (I(i)/ (param.H_I + I(i))));
% end

mu(ix_f) = param.mu_max*min((N(ix_f)'./(param.H_N+N(ix_f)')),(I(ix_f)./(param.H_I+I(ix_f))));

%--- Nutrients
uptake = param.alpha .* mu(ix_f) .*P(ix_f)';
recyc = param.eps .* param.alpha .* param.loss .* P(ix_f)';


% dNdt = -uptake + recyc - (Jdiff_N(2:(n+1))-Jdiff_N(1:n))/param.dx;
dNdt = -uptake + recyc - (J_N(2:(n+1))-J_N(1:n))/param.dx;


%--- Phtoplankton
dPdt = -(J_P(2:(n+1)) - J_P(1:n))/param.dx + mu .* P' - param.loss .* P';

dCdt = [dNdt  dPdt]';

end
end









