include("physical_constants.jl")

"Calculate delta0 from Tc and a custom value of the BCS constant."
function get_delta0(tc, bcs=1.76)
    return bcs*kb_eV*tc
end

"Calculate Tc from delta0 and a custom value of the BCS constant."
function get_tc(delta0, bcs=1.76)
    return delta0/(bcs*kb_eV)
end

"Calculate the BCS constant given a Tc and a delta0."
function get_bcs(tc, delta0)
    return delta0/(kb_eV*tc)
end

"Calculate the BCS coherence length from the Fermi velocity and delta0."
function get_ksi0(vf, delta0)
    return vf*hbar_eV/(π*delta0)
end

"Calculate the Fermi velocity from the BCS coherence length and delta0."
function get_vf(ksi0, delta0)
    return ksi0*π*delta0/hbar_eV
end

"""Calculate some unitless quantities needed for the computation from
physical data.

Parameters
----------
tc : float
    Critical temperature in Kelvin.

vf : float
    Fermi velocity in m/s.

mfp : float
    Mean free path in m.

london0 : float
    London penetration depth at zero temperature in m.

bcs : float (optional)
    The constant that gives the zero-temperature superconducting energy gap
    delta0 according to the equation delta0 = bcs*kB*Tc, where kB is
    Boltzmann's constant and Tc is the superconducting critical temperature.
    functionault value is the Bardeen-Cooper-Schrieffer value of 1.76.

ksi0 : float (optional)
    The BCS coherence length in m. ksi0 = hbar*vf/(pi*delta0), where hbar is
    Planck's constant divided by 2*pi, vf is the Fermi velocity, and delta0
    is the zero-temprature superconducting energy gap.

delta0 : float (optional)
    The zero-temperature superconducting energy gap in electron volts.
    delta0 = bcs*kB*Tc, where kB is Boltzmann's constant and Tc is the
    superconducting critical temperature. functionault value is the
    Bardeen-Cooper-Schrieffer value of 1.76.

Keyword Arguments
-----------------
verbose : int
    Whether to print out some debugging information. Very much a work in
    progress. functionault is 0.

    * 0: Don't print out any debugging information
    * 1: Print out minimal debugging information
    * 2: Print out full output from quad routines

Returns
-------
output_dict : dict
    output_dict contains the following keys: bcs, vf, x0, x1. These are all
    necessary inputs to impedance.cmplx_impedance."""
function init_from_physical_data(tc, vf, london0, mfp, bcs=1.76, ksi0=nothing, delta0=nothing, verbose=false)


    if delta0 == nothing
        delta0 = bcs*kb_eV*tc
    end

    if ksi0 == nothing
        ksi0 = vf*hbar_eV/(π*delta0)
        ksi0_calc = nothing
    else
        ksi0_calc = vf*hbar_eV/(π*delta0)
    end

    #Reduced lengths
    xk = ksi0/london0
    xm = mfp/london0

    # if verbose > 1
    #     if ksi0_calc != nothing
    #         print "calculated ksi0 = " + str(ksi0_calc*1e9)+" nm"
    #         print "supplied ksi0 = " + str(ksi0*1e9)+" nm"
    #     else
    #         print "calculated ksi0 = " + str(ksi0*1e9)+" nm"
    #     end
    #
    #     print "calculated delta0 = " + str(delta0*1e6) + " ueV"
    #     print "xk = ksi0/london0 = " + str(xk)
    #     print "xm = mfp/london0 = " + str(xm)
    #     print "xk/xm = ksi0/mfp = " + str(xk/xm)
    # end

    output_dict = Dict(:bc=>bcs,
                    :vf=>vf,
                    :xk=>xk,
                    :xm=>xm)

    return output_dict
end
