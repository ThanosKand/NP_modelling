% Data
L = 300; % [m]
param.D = 5*8.64; % [m^2/day]
param.u = 0.04*24; % [m/day]

param.t_range = 0:1000;
nopoints = 100;
param.dx = L / nopoints; %[m]
param.z = 0.5*param.dx:param.dx:(L-0.5*param.dx); 
param.c0 = [ones(1, nopoints)*100 ones(1, nopoints)*1e06]; % [mmol nutrient/m3] [cells/ m3] 

param.H_I = 20*86400; % [μmol photons m-2 day-1] half-saturation constant of light-limitead growth 
param.H_N = 0.04; % [mmol nutrient/m3]

param.I0 = 450*86400; %incident light intensity [μmol photons m-2 day-1]

param.kw = 0.045; % background turbidity [m-1]
param.kp =6e-12; % specific light attenuation of phytoplankton [m2/ cell]

param.mu_max = 0.04*24; %[1/day]
param.loss = 0.01*24; %[1/day]
param.alpha = 1*10^-9; % Nutrient content of phytoplankton [mmol nutrient/cell]
param.eps = 0.5; % Nutrient recycling coefficient []

%%
[t1, C] = NPZD(param);

n = length(param.z);
N = C(:, 1:n);
P = C(:, n+1:2*n);
%% Calculate light 
PP = P;

for i = 1:length(param.t_range)
    for j = 1:n
        
PP(i,:) = calclight(param.z,t1(i),P(i,:),param.dx,param.kp,param.kw,param.I0);

    end
end
%% Plot light in surface plot
figure;
hold on
set(gca,'Ydir','reverse')
surface(t1,param.z,PP')
shading interp
colorbar
xlabel('Time [days]');
ylabel('Depth [m]');
% xlim([0, max(param.t_range)])
title('Light intensity')

%% Depth profile in time

fig1 = figure;

subplot(3,1,1)
hold on
set(gca,'Ydir','reverse')
surface(t1,param.z,N')
shading interp
colorbar
xlabel('Time [days]');
ylabel('Depth [m]');
% xlim([0, max(param.t_range)])
title('Nutrients [mmol nutrient/m3]' )
 
subplot(3,1,3)
hold on
set(gca,'Ydir','reverse')
surface(t1,param.z,P')
shading interp
colorbar
xlabel('Time [days]');
ylabel('Depth [m]');
% xlim([0, max(param.t_range)])
title('Phytoplankton [cells/ m3] ')

subplot(3,1,2)
hold on
set(gca,'Ydir','reverse')
surface(t1,param.z,PP')
shading interp
colorbar
xlabel('Time [days]');
ylabel('Depth [m]');
% xlim([0, max(param.t_range)])
title('Light intensity [μmol photons m-2 day-1]')


%%

figure;
subplot(2,3,1)
hold on
plot(P(1,:), param.z, 'Linewidth', 1.5) % Half time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 0')
ylabel( 'Depth [m]')

%hold on
subplot(2,3,2)
hold on
plot(P(10,:), param.z, 'Linewidth', 1.5) % Half time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
title('Concentration of plytoplankton [cells/m3]')
xlabel(' t = 10 days')

%hold on
subplot(2,3,3)
hold on
plot(P(36,:), param.z, 'Linewidth', 1.5) % End time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 36 days')

%hold on
subplot(2,3,4)
hold on
plot(P(50,:), param.z, 'Linewidth', 1.5) % End time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 50 days')

%hold on
subplot(2,3,5)
plot(P(100,:), param.z, 'Linewidth', 1.5) % End time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 100 days')

%hold on
subplot(2,3,6)
hold on
plot(P(end,:), param.z, 'Linewidth', 1.5) % End time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 300 days')
%%

figure;
subplot(2,3,1)
hold on
plot(N(1,:), param.z, 'Linewidth', 1.5) % Half time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 0')
ylabel( 'Depth [m]')

%hold on
subplot(2,3,2)
hold on
plot(N(10,:), param.z, 'Linewidth', 1.5) % Half time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
title('Concentration of plytoplankton [cells/m3]')
xlabel(' t = 10 days')
title('Concentration of nutrients [mmol nutrient/m3]')

%hold on
subplot(2,3,3)
hold on
plot(N(33,:), param.z, 'Linewidth', 1.5) % End time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 20 days')

%hold on
subplot(2,3,4)
hold on
plot(N(50,:), param.z, 'Linewidth', 1.5) % End time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 50 days')

%hold on
subplot(2,3,5)
plot(N(100,:), param.z, 'Linewidth', 1.5) % End time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 100 days')

%hold on
subplot(2,3,6)
hold on
plot(N(end,:), param.z, 'Linewidth', 1.5) % End time
axis ij
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('t = 200 days')

%% Plot light
figure
hold on
plot(calclight(param.z,t1(1),P(1,:),param.dx,param.kp,param.kw,param.I0), param.z, 'Linewidth', 1.5)
plot(calclight(param.z,t1(36),P(36,:),param.dx,param.kp,param.kw,param.I0), param.z, 'Linewidth', 1.5)
plot(calclight(param.z,t1(end),P(end,:),param.dx,param.kp,param.kw,param.I0), param.z, 'Linewidth', 1.5)

xlabel('Light intensity [μmol photons m-2 day-1]')
ylabel( 'Depth [m]')
axis ij
legend('t = 0 days', 't = 36 days', 't = end')
title('Light')


%%
PPP = PP;
nutr_lm = zeros(length(param.t_range),n);
light_lm = zeros(length(param.t_range),n);
for i = 1:length(param.t_range)
    for j = 1:n
        
nutr_lm(i,j) = (N(i,j)/ (param.H_N + N(i,j)));
light_lm(i,j) = (PP(i,j) / (param.H_I + PP(i,j)));

if nutr_lm(i,j) < light_lm(i,j)
    PPP(i,j) = 0;
else
    PPP(i,j) = 1;
end

    end
end

figure;
hold on
set(gca,'Ydir','reverse')
surface(t1,param.z,PPP')
shading interp
colorbar
xlabel('Time [days]');
ylabel('Depth [m]');
% xlim([0, max(param.t_range)])
title('Limiting factor')

%%

target_time = 596;

figure;
line(P(target_time,:), param.z, 'Linewidth', 1.5) % Half time
ax1 = gca; % current axes
axis ij
ylabel( 'Depth [m]')
xlabel('Phytoplankton concentration [cells/ m3]')
%ax1.XColor = '';
%ax1.YColor = 'b';
% ax1.XAxisLocation = 'top';

ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos, 'XAxisLocation','top', 'YAxisLocation','right', 'Color','none');
line(nutr_lm(target_time,:), param.z, 'Parent', ax2, 'Linewidth', 1.2, 'Color','r', 'LineStyle', '--')
axis ij
line(light_lm(target_time,:), param.z, 'Parent', ax2, 'Linewidth', 1.2, 'Color','k', 'LineStyle', '-.')
axis ij
xlabel('Limitation by light or nutrients')


legend('Nutrients','Light', 'Location', 'southeast')
title('t = 50 days')
% ylabel( 'Depth [m]')

%%
function [t, C] = NPZD(param)

n = length(param.z);

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[t,C] = ode45(@derivative, param.t_range, param.c0,opts);

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

mu = zeros;
for i = 1:n
    mu(i) = param.mu_max * min((N(i)/ (param.H_N + N(i))), (I(i)/ (param.H_I + I(i))));
end

% mu(ix_f) = param.mu_max*min((N(ix_f)'./(param.H_N+N(ix_f)')),(I(ix_f)./(param.H_I+I(ix_f))));


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






