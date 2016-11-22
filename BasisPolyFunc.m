function U=BasisPolyFunc(x,P)
    U=basisPoly(x,P(1:end)).*func(x);
end