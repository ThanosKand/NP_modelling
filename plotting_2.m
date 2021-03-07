%% Depth profile in time
function plotting_2(N,P, param,t)
figure

subplot(2,1,1)
hold on
set(gca,'Ydir','reverse')
surface(t,param.z,N')
shading interp
colorbar
xlabel('Time [days]');
ylabel('Depth [m]');
% xlim([0, max(param.t_range)])
title('Nutrients')

subplot(2,1,2)
hold on
set(gca,'Ydir','reverse')
surface(t,param.z,P')
shading interp
colorbar
xlabel('Time [days]');
ylabel('Depth [m]');
% xlim([0, max(param.t_range)])
title('Phytoplankton')

end
