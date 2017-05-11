"""Calculate the R integral from Popel divided by x, (x = q*L0, L0 = zero-temp London depth).

Parameters
----------
a : float

b : float

x : float

Returns
-------
r : float
    The R integral from Popel, divided by x.

Note
----
See R. Pöpel (1989), doi: 10.1063/1.343622 for more details."""
function intR(a, b, x)

    z2 = a^2+b^2

    #for small x
    if x < 0.01*sqrt(z2)
        r = b/(3.0*z2) #This is really r/x

    #for large x
    elseif x > 100*sqrt(z2)
        r = (π*(1+(b^2-a^2)/x^2)/4-b/x)/x #This is really r/x

    #in between x
    else
        #calculate all the terms of r
        r = (1/x^2)*(-0.5*b*x+0.25*a*b*log((z2+x^2+2*a*x)/(z2+x^2-2*a*x))+
        0.25*(x^2+b^2-a^2)*atan2(2*b*x,(z2-x^2)))/x #This is really r/x
    end

    return r
end

"""Calculate the R integral from Popel divided by x, (x = q*L0, L0 = zero-temp London depth).

Parameters
----------
a : float

b : float

x : float

Returns
-------
s : float
    The S integral from Popel, divided by x.

Note
----
See R. Pöpel (1989), doi: 10.1063/1.343622 for more details."""
function intS(a, b, x)

    z2 = a^2+b^2

    #for small x
    if x < 0.01*sqrt(z2)
        s = a/(3.0*z2) #This is really s/x

    #for large x
    elseif x > 100*sqrt(z2)
        s = (a/x - a*b*π/(2*x^2))/x #This is really s/x

    #in between x
    else
        #calculate all the terms of s
        s = (1/x^2)*(0.5*(a*x)+
            0.125*(x^2+b^2-a^2)*log((z2+x^2+2*a*x)/(z2+x^2-2*a*x))-
            0.5*b*a*atan2(2*b*x,(z2-x^2)))/x #This is really s/x
    end

    return s
end

"""Similar to intS, but divided by a. This is specifically useful for when
a=0, as lim(a->0) intS/a should not blow up!"""
function intS_over_a(b, x)

    z2 = b^2

    #for small x
    if x < 0.01*sqrt(z2)
        s = 1.0/(3.0*z2) #This is really s/(x*a)

    #for large x
    elseif x > 100*sqrt(z2)
        s = (1.0/x - b*π/(2*x^2))/x #This is really s/(x*a)

    #in between x
    else
        #calculate all the terms of s
        s = (1/x^2)*(x-0.5*b*atan2(2*b*x,(z2-x^2)))/x #This is really s/(x*a)
    end

    return s
end
