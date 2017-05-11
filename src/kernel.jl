include("volume_integrals.jl")
include("fermi_functions.jl")

"""Energy integral of real part of K, limits -dr to 0."""
function reKlint1_a1(u, args)
    #u-substitution, u^2 = en+dr
    #from -dr to 0, only fires if f >= 2*dr
    #ulims are: 0, sqrt(dr)

    x, a0, b0, tr, fr, dr, bcs = args

    en = u^2-dr
    e1 = u*sqrt(2*dr-u^2)
    e2 = sqrt((en+fr)^2-dr^2)
    a1 = a0*e1
    a2 = a0*e2

    return 2*(1-2*fermi_fun(en+fr, tr/bcs))*(((en^2+dr^2+en*fr)/(sqrt(2*dr-u^2)*e2))*intR(a2, a1+b0, x)+u*intS(a2, a1+b0, x))
end

"""Energy integral of real part of K, limits dr-fr to dr-0.5*fr."""
function reKlint1_a2(u, args)
    #u-substitution, u^2 = (en+fr)-dr
    #from dr-fr to dr-0.5*fr
    #ulims are: 0, sqrt(0.5*fr)

    x, a0, b0, tr, fr, dr, bcs = args

    en = u^2-fr+dr
    e1 = sqrt(abs(en^2-dr^2))
    e2 = u*sqrt(u^2+2*dr)
    a1 = a0*e1
    a2 = a0*e2

    return 2*(1-2*fermi_fun(en+fr, tr/bcs))*(((en^2+dr^2+en*fr)/(e1*sqrt(u^2+2*dr)))*intR(a2, a1+b0, x)+u*intS(a2, a1+b0, x))
end

"""Energy integral of real part of K, limits are max(0, dr-0.5*fr) to dr."""
function reKlint1_b(u, args)
    #u-substitution, u^2 = dr-en
    #from 0 to dr OR dr-0.5*fr to dr, depending on whether fr > 2*dr
    #ulims are: 0 to sqrt(dr) OR 0, sqrt(0.5*fr)

    x, a0, b0, tr, fr, dr, bcs = args

    en = dr-u^2
    e1 = u*sqrt(2*dr-u^2)
    e2 = sqrt((en+fr)^2-dr^2)
    a1 = a0*e1
    a2 = a0*e2

    return 2*(1-2*fermi_fun(en+fr, tr/bcs))*(((en^2+dr^2+en*fr)/(sqrt(2*dr-u^2)*e2))*intR(a2, a1+b0, x)+u*intS(a2, a1+b0, x))
end

"""Energy integral of real part of K, limits are dr-fr to -fr/2."""
function reKlint2_a(u, args)
    #u-substitution, u^2 = en+fr-dr
    #from dr-fr to -fr/2, only fires if fr > 2*dr
    #ulims are: 0, sqrt(0.5*fr-dr)

    x, a0, b0, tr, fr, dr, bcs = args

    en = u^2-fr+dr

    e1 = sqrt(en^2-dr^2)
    e2 = sqrt((en+fr)^2-dr^2)
    gfun = (en^2+dr^2+en*fr)/(e1*sqrt(u^2+2*dr))

    return (1-2*fermi_fun(en+fr, tr/bcs))*((gfun+u)*intS(a0*(e2-e1), b0, x)
                                                       -(gfun-u)*intS(a0*(e2+e1), b0, x))
end

"""Energy integral of real part of K, limits are -fr/2 to -dr."""
function reKlint2_b(u, args)
    #u-substitution, u^2 = -en-dr
    #from -fr/2 to -dr, only fires if fr > 2*dr
    #ulims are: 0, sqrt(0.5*fr-dr)

    x, a0, b0, tr, fr, dr, bcs = args

    en = -u^2-dr

    e1 = sqrt(en^2-dr^2)
    e2 = sqrt((en+fr)^2-dr^2)
    gfun = (en^2+dr^2+en*fr)/(sqrt(u^2+2*dr)*e2)

    return (1-2*fermi_fun(en+fr, tr/bcs))*((gfun+u)*intS(a0*(e2-e1), b0, x)
                                                       -(gfun-u)*intS(a0*(e2+e1), b0, x))
end

"""Energy integral of real part of K, limits are dr to infinity."""
function reKlint3(u, args)
    #u-substitution, u^2 = en-dr
    #from dr to infinity
    #ulims are: 0, 31.6 <-- this is big enough, rounding errors happen at big numbers

    x, a0, b0, tr, fr, dr, bcs = args

    en = u^2+dr

    #Fixes the divide by zero error nicely
    #This works because intS goes to zero like 'a' as 'a'->0.
    #This leading behavior cancels out the part that makes the integrand blow up at the lower limit.
    if (fr == 0) && (u == 0)
        tplus = 1-fermi_fun(en, tr/bcs)-fermi_fun(en+fr, tr/bcs)

        return -2*tplus*a0*2*intS_over_a(b0, x)*(en^2+dr^2+en*fr)/sqrt(2*dr)

    else

        e1 = sqrt((en+dr)*(en-dr))
        e2 = sqrt(((en+fr)+dr)*((en+fr)-dr))
        gfun = (en^2+dr^2+en*fr)/(sqrt(u^2+2*dr)*e2)
        tplus = 1-fermi_fun(en, tr/bcs)-fermi_fun(en+fr, tr/bcs)
        tminus = (1-2*fermi_fun(fr, tr/bcs))*(1-(1-2*fermi_fun(en, tr/bcs))*(1-2*fermi_fun(en+fr, tr/bcs)))

        return -2*tplus*(gfun-u)*intS(a0*(e2+e1), b0, x) + tminus*(gfun+u)*intS(a0*(e2-e1), b0, x)
    end
end


"""Energy integral of imaginary part of K, limits are dr-fr to -dr."""
function imKlint1_a(u, args)
    #from dr-fr to -dr
    #ulims are: 0, sqrt(0.5*fr-dr)

    x, a0, b0, tr, fr, dr, bcs = args

    en = u^2-fr+dr

    e1 = sqrt(en^2-dr^2)
    e2 = sqrt((en+fr)^2-dr^2)
    gfun = (en^2+dr^2+en*fr)/(e1*sqrt(u^2+2*dr))

    return -(1-2*fermi_fun(en+fr, tr/bcs))*((gfun+u)*intR(a0*(e2-e1), b0, x)+(gfun-u)*intR(a0*(e2+e1), b0, x))
end

"""Energy integral of imaginary part of K, limits are dr-fr to -dr."""
function imKlint1_b(u, args)
    #from dr-fr to -dr
    #ulims are: 0, sqrt(0.5*fr-dr)

    x, a0, b0, tr, fr, dr, bcs = args

    en = -u^2-dr

    e1 = sqrt(en^2-dr^2)
    e2 = sqrt((en+fr)^2-dr^2)
    gfun = (en^2+dr^2+en*fr)/(sqrt(u^2+2*dr)*e2)

    return -(1-2*fermi_fun(en+fr, tr/bcs))*((gfun+u)*intR(a0*(e2-e1), b0, x)+(gfun-u)*intR(a0*(e2+e1), b0, x))
end

"""Energy integral of imaginary part of K, limits are dr to infinity."""
function imKlint2(u, args)
    #u-substitution, u^2 = en-dr
    #from dr to infinity
    #ulims are 0, infinity

    x, a0, b0, tr, fr, dr, bcs = args

    en = u^2+dr

    e1 = sqrt(en^2-dr^2)
    e2 = sqrt((en+fr)^2-dr^2)
    gfun = (en^2+dr^2+en*fr)/(sqrt(u^2+2*dr)*e2)
    tminus = (1-2*fermi_fun(fr, tr/bcs))*(1-(1-2*fermi_fun(en, tr/bcs))*(1-2*fermi_fun(en+fr, tr/bcs)))

    return tminus*((gfun+u)*intR(a0*(e2-e1), b0, x)+(gfun-u)*intR(a0*(e2+e1), b0, x))
end


"""Calculate the Kl integral over energy. Kl = K*mfp^2 where K is the Mattis-Bardeen Kernel.

Parameters
----------
tr : float
    Reduced temperature tr = T/Tc.

fr : float
    Reduced frequency fr = h*f/delta0

x : float
    Reduced momentum x = q*london0 where london0 is the zero-temperature
    London penetration depth and q is the coordinate resulting from a
    fourier transform of position.

xk : float
    BCS coherence length (ksi0) divided by mean free path (mfp). ksi0 =
    hbar*vf/(π*delta0), where hbar is Planck's constant divided by 2*π, vf
    is the Fermi velocity, and delta0 is the zero-temprature superconducting
    energy gap.

xm : float
    Mean free path (mfp) divided by the zero-temperature London penetration
    depth (london0).

dr : float
    Reduced BCS energy gap dr = delta(T)/delta(T=0) where delta is the
    superconducting energy gap.

Keyword Arguments
-----------------
verbose : int
    Whether to print out some debugging information. Very much a work in
    progress. functionault is 0. Currently does nothing.

Returns
-------
Kl : complex float
    Kl = K*london0^2 where K is the complex Mattis-Bardeen Kernel as a
    function of x = q*london0, london0 is the zero-temperature London
    penetration depth and q is the coordinate resulting from a fourier
    transform of position.

Note
----
See Mattis and Bardeen (1958), Pöpel (1989), or Gao (2008)."""
function cmplx_kernel(tr, fr, x, xk, xm, dr, bcs, verbose=false)

    #a0 and b0
    a0 = 1.0/(xk*π)

    b0 = 1.0/xm

    #If temperature is >= Tc, then skip all the integration, and
    #use the normal Kernel
    if tr >= 1

        #Calculate the actual number out front (1/x dependance included in R and S)
        prefactor = 3*fr*a0

        reKl = intS(fr*a0, b0, x)
        imKl = intR(fr*a0, b0, x)

    else
        #arguments to pass the integrand functions
        iargs = (x, a0, b0, tr, fr, dr, bcs,)

        #Calculate the actual number out front (1/x dependance included in R and S)
        prefactor = 3*a0

        #Now run the integrals

        if (fr == 0)
            #Only three of the integrals actually
            #contribute. reKl1_a2 and reKl1_b have an analytic expression,
            #and reKl3 has a limit built into it.

            #reKl1_a2 = reKl1_b = Fermi*intR*0.25*π
            reKl1 = 2*(1-2*fermi_fun(dr, tr/bcs))*intR(0, b0, x)*0.5*π*dr
            reKl1err = 0
            reKl2 = 0
            reKl2err = 0
            imKl1 = 0
            imKl1err = 0
            imKl2 = 0
            imKl2err = 0

        elseif (fr/dr < 2)
            #from (dr-fr) to dr (inside the gap), limits look wierd due to u-substitution
            reKl1a, reKl1aerr = quadgk(u->reKlint1_a2(u, iargs), 0, sqrt(0.5*fr))
            reKl1b, reKl1berr = quadgk(u->reKlint1_b(u, iargs), 0, sqrt(0.5*fr))

            reKl1 = reKl1a + reKl1b
            reKl1err = sqrt(reKl1aerr^2 + reKl1berr^2)

            #from dr-fr to -dr (below the gap), limits look wierd due to u-substitution
            #These ones are zero if the photons aren't big enough to break pairs
            reKl2 = 0
            reKl2err = 0
            imKl1 = 0
            imKl1err = 0
        else
            #from -dr to dr (inside the gap), limits look wierd due to u-substitution

            reKl1a, reKl1aerr = quadgk(u->reKlint1_a1(u, iargs), 0, sqrt(dr))
            reKl1b, reKl1berr = quadgk(u->reKlint1_b(u, iargs), 0, sqrt(dr))

            reKl1 = reKl1a + reKl1b
            reKl1err = sqrt(reKl1aerr^2 + reKl1berr^2)

            #from dr-fr to -dr (below the gap), limits look wierd due to u-substitution
            reKl2a, reKl2aerr = quadgk(u->reKlint2_a(u, iargs), 0, sqrt(0.5*fr-dr))
            reKl2b, reKl2berr = quadgk(u->reKlint2_b(u, iargs), 0, sqrt(0.5*fr-dr))

            reKl2 = reKl2a+reKl2b
            reKl2err = sqrt(reKl2aerr^2+reKl2berr^2)

            imKl1a, imKl1aerr = quadgk(u->imKlint1_a(u, iargs), 0, sqrt(0.5*fr-dr))
            imKl1b, imKl1berr = quadgk(u->imKlint1_b(u, iargs), 0, sqrt(0.5*fr-dr))

            imKl1 = imKl1a+imKl1b
            imKl1err = sqrt(imKl1aerr^2+imKl1berr^2)
        end

        #from dr to infinity (above the gap), limits look wierd due to u-substitution
        #Wierd things happen numerically at infinity. 31.6 is far enough. equals about 1000 gaps.
        reKl3, reKl3err = quadgk(u->reKlint3(u, iargs), 0, 31.6)
        imKl2, imKl2err = quadgk(u->imKlint2(u, iargs), 0, Inf)

        reKl = reKl1+reKl2+reKl3
        imKl = imKl1+imKl2
    end

#         if verbose
#             print 'x = '+str(x)+'\n'
#             print 'prefactor = '+str(prefactor)+'\n'
#             print 'First reK integral = '+str(reKl1)+' +/- '+str(reKl1err) + '\n'
#             print 'Second reK integral = '+str(reKl2)+' +/- '+str(reKl2err) + '\n'
#             print 'Third reK integral = '+str(reKl3)+' +/- '+str(reKl3err) + '\n'
#             print 'First imK integral = '+str(imKl1)+' +/- '+str(imKl1err) + '\n'
#             print 'Second imK integral = '+str(imKl2)+' +/- '+str(imKl2err) + '\n'

    return prefactor*(reKl+1im*imKl)
end
