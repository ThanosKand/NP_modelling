function I = calclight(z,t,C,dx,kp,kw,I0)
int=zeros;

for i = 2:length(z)
    int(i) = int(i-1) + kp*C(i)*dx;
   
end

% lat = 85;
% g1 = (1-0.8.*sin((pi.*lat)/180).*cos(2.*pi.*(t./365)));
% g1 = (1-0.8.*sin((pi.*lat)/180).*cos(2.*pi.*(t./365)));
% g1 = 1+cos(2.*pi.*(t./365));

% I0 = I0 * g1;

I = I0.*(exp(-kw*z) .* exp(-int));
I(1) = I0.*(exp(-kw.*z(1)) .* exp(-kp.*C(1).*dx));

% I = I0*exp(-kw*z - kp*cumsum(C')*dx);

end



% zzz = 1:50;
% 
% cc = 1000 * (exp(-2) * exp(-0.2 * zzz));
% % cc = 1000 * (exp(-0.2 * zzz))
% figure
% plot(zzz, cc)
