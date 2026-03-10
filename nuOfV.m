function nu = nuOfV(v,p)
    if(v==0)
        nu=0;
        return;
    end
    v0 = p(4);
    gamma = p(5);
    beta = p(3);
    alpha = p(2);
    x = v./v0;
    numerator = p(1)*(alpha+beta);
    denominator = beta*(x.^-alpha)+alpha*x.^beta;
    nu = numerator./denominator.*v.^gamma;
end