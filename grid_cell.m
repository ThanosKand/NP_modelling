% Data
L = 300; % [m]
% D = 1; % [m^2/day]
D = 5*8.64;
% u = 0.2; % 
u = 0.04*24; % [m/day]
% t_range = [0,50]; % [days]
t_range = 0:1500;

% ASK!!! How do we decide the best # of points/ delta
% delta affects the diffusion parameter -> if many points, more focused/ if
% less points, more spread
nopoints = 50;
dx = L / nopoints; %[m]
% dx = 0.3; % [m]

% Grid cells
% z = linspace(0.5*dx, L-0.5*dx, nopoints);
z = 0.5*dx:dx:(L-0.5*dx); % ASK !!! Why do we start from 0.5? Are we sure that we use the middle of each cell?

% Initial conditions
% c0 = zeros(1,500); 
%c0 = zeros(1,nopoints); 
% c0(z<L/5) = 20;
c0 = ones(1, nopoints)*20;
%%
% We want to solve: dΦ/dt = - d(Ja + Jd)/ dz
[t, C] = NPZD(dx,D,u,z,c0,t_range);

%% Depth profile in time

fig1 = figure;
set(gca,'Ydir','reverse')
surface(t,z,C')
shading interp

h = colorbar;
ylabel(h, 'Nutrient concentration [mg/L]')
xlabel('Time [days]');
ylabel('Depth [m]');
% xlim([0, max(tspan1)])
title('dz = 0.333 m')
 
%% Depth against concentration in different times (start / end)

fig2 = figure;
hold on
plot(C(1,:), z, 'Linewidth', 1.5) % Start time
plot(C(t(round(length(t)/2)),:), z, 'Linewidth', 1.5) 
plot(C(end,:), z, 'Linewidth', 1.5) % End time
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
xlabel('Population density [cells/ m-3]')
ylabel( 'Depth [m]')
axis ij
xlim([0,50])
legend('t = 0', 't = t/2', 't = end', 'Location','southwest')
 
%% Concentration of each grid cell against time

fig3 = figure;
plot(t,C(:,:), 'Linewidth', 1.5) 
ylabel('Nutrient concentration in each grid cell [mg/L]')
xlabel('Time [days]')
% ylim([0,20])
title('dz = 0.333 m')

%% Check steps

figure

plot(C(end-50:end,:), z, '-', 'Linewidth', 1.1)




%%
figure
hold on 
for nopoints_n = [500, 300, 200, 100, 50, 20]

dx = L / nopoints_n; %[m]
z = 0.5*dx:dx:(L-0.5*dx); % ASK !!! Why do we start from 0.5? Are we sure that we use the middle of each cell?
c0 = ones(1, nopoints_n)*20 ;
    
[t, C] = NPZD(dx,D,u,z,c0,t_range);

plot(C(end,:), z, '-', 'Linewidth', 1.1)
axis ij
end
legend()
xlabel('Concentration (X/m^3)')
ylabel('Depth (m)')

%%
figure
hold on 
i = 1;
for nopoints_n = [900, 700, 500, 300, 150, 50]

dx = L / nopoints_n; %[m]
z = 0.5*dx:dx:(L-0.5*dx); % ASK !!! Why do we start from 0.5? Are we sure that we use the middle of each cell?
c0 = ones(1, nopoints_n)*20 ;
    
[t, C] = NPZD(dx,D,u,z,c0,t_range);
hold on
subplot(3,2,i)
plot(t, C(:,:), '-', 'Linewidth', 1.1)
axis ij
i = i + 1;
title(num2str(dx))
end
% legend()
xlabel('Concentration (X/m^3)')
ylabel('Depth (m)')






%%

function [t, C] = NPZD(dx,D,u,z,c0,t_range)

n = length(z);

[t,C] = ode45(@derivative, t_range, c0);

function dCdt = derivative(~,C)
   
Ja = zeros(n+1,1);
Jd = zeros(n+1,1);

% n is the left boundary of the last cell
% n+1 is the right boundary of the last cell

for i = 2:n
    Ja(i) = u*C(i-1);
    Jd(i) = - D*(C(i) - C(i-1))/dx;
end

% Boundaries
Ja(1) = 0;
Jd(1) = 0;

% ASK!!! Do we need to set Ja(50) and Jd(50) equals to 0 for "closed system"??
% J(50) is what comes in from advection and diffsusin fluxes?
Ja(n+1) = 0;
Jd(n+1) = 0;

J = Ja + Jd;

% Initialize
dCdt = zeros(n,1);

for j = 1:n
 
  dCdt(j) = -(J(j+1) - J(j))/dx ; 
    
end

end
end



