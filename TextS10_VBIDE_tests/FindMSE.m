function Error = FindMSE(s,intensity,field)

J = length(s);
s_sub = linspace(s(1),s(end),25);
[S1_sub,S2_sub] = meshgrid(s_sub);
[s1,s2] = meshgrid(s);
bias = -2;

N = size(field,3);
Error = zeros(N,1);

for i = 1:N
    Field_original = @(s1,s2) interp2(S1_sub,S2_sub,field(:,:,i),s1,s2);
    Error(i) = sum(sum(((exp(bias + Field_original(s1,s2)) -intensity(:,:,i)).^2)))/J^2;
end