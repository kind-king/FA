function U=basisPoly(x,P)
    U=0;
    for n=0:length(P)-1
        U=P(length(P)-n).*x.^n+U;
    end
end