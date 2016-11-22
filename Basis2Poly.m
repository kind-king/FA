function U=Basis2Poly(x,P1,P2)
    U=basisPoly(x,P1(1:end)).*basisPoly(x,P2(1:end));
end