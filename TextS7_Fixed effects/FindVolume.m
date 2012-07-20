function Volume = FindVolume(s1,s2,GaussKernel)

% Find the volume of a 2-D GRBF in spatial grid defined by (s1,s2) numerically
Volume = 0;
ds1 = mean(diff(s1));
ds2 = mean(diff(s2));

[S1,S2] = meshgrid(s1,s2);
Volume = Volume + trapz(s2,trapz(s1,GaussKernel(S1,S2,0,0)));
