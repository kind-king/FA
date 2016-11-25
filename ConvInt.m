function T=ConvInt(P,p1,p2)
T=0;
P=conv(P(1:end),P(1:end));
for k=0:length(P)-1
    T=P(length(P)-k).*intPoly(p1,p2,k+1)+T;
    %U=P(length(P)-n).*x.^n+U;
end
end