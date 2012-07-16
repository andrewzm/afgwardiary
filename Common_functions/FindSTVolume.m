function Volume = FindSTVolume(s1,s2,t,GaussKernel,mu)

% Find the volume (numerically) of a 2D GRBF using the vectors s1,s2

Volume = 0;
dt = mean(diff(t));
ds1 = mean(diff(s1));
ds2 = mean(diff(s2));

[S1,S2] = meshgrid(s1,s2);

for i = 1:length(t)
    Volume = Volume + trapz(s2,trapz(s1,GaussKernel(S1,S2,t(i),mu(1),mu(2),mu(3))))*dt;
end