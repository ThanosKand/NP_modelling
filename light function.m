lat1 = 85;
lat2 = 5;
tt = 0:365;
syms t

II = 450*86400;

g1 = II * (1-0.8.*sin((pi.*lat1)/180).*cos(2.*pi.*(t./365)));
g2 = II * (1-0.8.*sin((pi.*lat2)/180).*cos(2.*pi.*(t./365)));

gg1 = subs(g1,t, tt)
gg2 = subs(g2,t, tt)

figure
hold on
plot(tt, gg1,'-', 'Linewidth', 1.5);
plot(tt, gg2,'-', 'Linewidth', 1.5);
yline(II, '--')
xlabel('Time [days]')
ylabel('Light intensity (surface) [μmol photons m-2 day-1]')
legend('Latitude = 85 degrees', 'Latitude = 5 degrees', 'Average incident light (I0)')



tick = ['-', '--', ':', '-.'];