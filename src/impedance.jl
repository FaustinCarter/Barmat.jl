include("physical_constants.jl")
include("tools.jl")
include("gap_functions.jl")
include("kernel.jl")

"""Return Z at several temperatures or frequencies.

Parameters
----------
input_vector : list-like
    List of either reduced temperatures (T/Tc) or reduced frequencies
    (h*f/delta0), where Tc is the superconducting critical temperature, h is
    Planck"s constant and delta0 is the zero-temperature superconducting
    energy gap.

tc : float
    Critical temperature in Kelvin.

vf : float
    Fermi velocity in m/s.

london0 : float
    London penetration depth at zero temperature in m.

axis : "string" (optional)
    Acceptable values are ``"temperature" or "frequency" or "mean free
    path"``. Specifies what the units of the ``input_vector`` parameter are.
    Default is ``"temperature"``.

Keyword Arguments
-----------------
mfp, tr, fr : float (required)
    If ``axis == "temperature"`` must specify reduced frequency value fr =
    h*f/delta0 and mean-free-path in meters. If ``axis == "frequency"`` must specify a reduced
    temperature value tr = T/Tc and mean-free-path in meters. If ``axis == "mean free path"`` must
    specify both fr and tr.

bcs : float (optional)
    The constant that gives the zero-temperature superconducting energy gap
    delta0 according to the equation delta0 = bcs*kB*Tc, where kB is
    Boltzmann"s constant and Tc is the superconducting critical temperature.
    Default value is the Bardeen-Cooper-Schrieffer value of 1.76.

gap : python function or string (optional)
    Python function that gives a value for the reduced superconducting
    energy gap deltar(T) = delta(T)/delta0, where delta is the
    superconducting energy gap and delta0 = delta(T=0). Function signature
    is float(float) with return value between zero and one. Default is
    tabulated values from Muhlschlegel (1959) via the deltar_bcs function.
    Optionally, one may pass the string ``"cos"`` to use the built-in cosine
    approximation of the gap.


output_depths : bool (optional)
    Sets the output units. False returns complex impedance in Ohms. True
    converts the complex impedance to a skin-depth (real part) and a
    superconducting penetration depth (imaginary part), both in meters.
    Default is False.

boundary : string (optional)
    Options are ``"diffuse"``/``"d"`` or ``"specular"``/``"s"``. Determines
    whether the impedance calculation assumes diffuse or specular scattering
    at the boundaries of the superconductor. Default is ``"diffuse"``.

verbose : int
    Whether to print out some debugging information. Very much a work in
    progress. Default is 0.

    * 0: Don"t print out any debugging information
    * 1: Print out minimal debugging information
    * 2: Print out full output from quad routines

Returns
-------
zs : numpy array
    The complex impedance calculated at each value of input_vector.

Note
----
fr = 0 may return slightly wrong number for specular boundary conditions.
Not sure why."""
function get_Zvec(input_vector, tc, vf, london0, axis="temperature"; kwargs...)

    allowed_axes = ["temperature", "t", "frequency", "f", "mean free path", "mfp", "l"]

    @assert axis in allowed_axes "Invalid axis."

    kwargs = Dict(kwargs)

    #This is the value that specifies delta0 = bcs*kB*Tc
    bcs = pop!(kwargs, :bcs, 1.76)
    delta0 = bcs*kb_eV*tc

    #Kwargs to pass to the cmplx_impedance function
    zs_kwargs = Dict()

    #Allow for passing in a custom gap function
    gap = pop!(kwargs, :gap, nothing)
    if gap != nothing
        if gap == "cos"
            zs_kwargs[:gap] = deltar_cos
        else
            zs_kwargs[:gap] = gap
        end
    end

    verbose = pop!(kwargs, :verbose, 0)

    zs_kwargs[:verbose] = verbose
    #Optionally can output penetration/skin depth in meters instead of Ohms
    output_depths = pop!(kwargs, :output_depths, false)
    zs_kwargs[:output_depths] = output_depths

    #See if specular reflection is requested
    boundary = pop!(kwargs, :boundary, "diffuse")
    @assert boundary in ["diffuse", "d", "specular", "s"] "Invalid boundary type."
    zs_kwargs[:boundary] = boundary

    if axis in ["mfp", "mean free path", "l"]
        @assert haskey(kwargs, :fr) "Must supply reduced frequency"
        fr = kwargs[:fr]

        @assert haskey(kwargs, :tr) "Must supply reduced temperature"
        tr = kwargs[:tr]

        mfps = input_vector

        zs = Array{Complex128}(length(input_vector))

        for (ix, mfp) in enumerate(mfps)
            #Convert physical data to params
            params_dict = init_from_physical_data(tc, vf, london0, mfp, bcs)

            #Add the parameters from physical data into the kwargs dict
            merge!(zs_kwargs, params_dict)

            xk = pop!(zs_kwargs, :xk, nothing)
            xm = pop!(zs_kwargs, :xm, nothing)
            vf = pop!(zs_kwargs, :vf, nothing)

            #Calculate the next impedance
            zs[ix] = cmplx_impedance(tr, fr, tc, xk, xm, vf; zs_kwargs...)
        end

    else
        @assert haskey(kwargs, :mfp) "Must supply mean free path"
        mfp = kwargs[:mfp]

        #Convert physical data to params
        params_dict = init_from_physical_data(tc, vf, london0, mfp, bcs)

        #Add the parameters from physical data into the kwargs dict
        merge!(zs_kwargs, params_dict)

        xk = pop!(zs_kwargs, :xk, nothing)
        xm = pop!(zs_kwargs, :xm, nothing)
        vf = pop!(zs_kwargs, :vf, nothing)

        if axis in ["temperature", "t"]
            @assert haskey(kwargs, :fr) "Must supply reduced frequency"
            fr = kwargs[:fr]

            trs = input_vector
            zs = Array{Complex128}(length(input_vector))

            for (ix, tr) in enumerate(trs)
              zs[ix] = cmplx_impedance(tr, fr, tc, xk, xm, vf; zs_kwargs...)
            end

        elseif axis in ["frequency", "f"]
            @assert haskey(kwargs, :tr) "Must supply reduced temperature"
            tr = kwargs[:tr]

            frs = input_vector
            zs = Array{Complex128}(length(input_vector))

            for (ix, fr) in enumerate(frs)
              zs[ix] = cmplx_impedance(tr, fr, tc, xk, xm, vf; zs_kwargs...)
            end
        end
    end

    return zs
end

"""Calculate the complex surface impedance (Z) of a superconductor at a
given temperature and frequency.

Parameters
----------
tr : float
    Reduced temperature tr = T/Tc, where T is temperature and Tc is the
    superconducting critical temperature.

fr : float
    Reduced frequency fr = h*f/delta0 where h is Planck"s constant, f is the
    frequency in Hz, and delta0 is the zero-temperature superconducting
    energy gap.

tc : float
    Superconducting critical temperature in Kelvin.

xk : float
    BCS coherence length (ksi0) divided by mean free path (mfp). ksi0 =
    hbar*vf/(pi*delta0), where hbar is Planck"s constant divided by 2*pi, vf
    is the Fermi velocity, and delta0 is the zero-temprature superconducting
    energy gap.

xm : float
    Mean free path (mfp) divided by the zero-temperature London penetration
    depth (london0).

vf : float
    Fermi velocity in m/s.


Keyword Arguments
-----------------
bcs : float (optional)
    The constant that gives the zero-temperature superconducting energy gap
    delta0 according to the equation delta0 = bcs*kB*Tc, where kB is
    Boltzmann"s constant and Tc is the superconducting critical temperature.
    Default value is the Bardeen-Cooper-Schrieffer value of 1.76.

boundary : string (optional)
    Options are ``"diffuse"``/``"d"`` or ``"specular"``/``"s"``. Determines
    whether the impedance calculation assumes diffuse or specular scattering
    at the boundaries of the superconductor. Default is ``"diffuse"``.

gap : python function or string (optional)
    Python function that gives a value for the reduced superconducting
    energy gap deltar(T) = delta(T)/delta0, where delta is the
    superconducting energy gap and delta0 = delta(T=0). Function signature
    is float(float) with return value between zero and one. Default is
    tabulated values from Muhlschlegel (1959) via the deltar_bcs function.
    Optionally, one may pass the string ``"cos"`` to use the built-in cosine
    approximation of the gap.

output_depths : bool (optional)
    Sets the output units. False returns complex impedance in Ohms. True
    converts the complex impedance to a skin-depth (real part) and a
    superconducting penetration depth (imaginary part), both in meters.
    Default is False.

verbose : int
    Whether to print out some debugging information. Very much a work in
    progress. Default is 0.

    * 0: Don"t print out any debugging information
    * 1: Print out minimal debugging information
    * 2: Print out full output from quad routines

Returns
-------
Z: The complex surface impedance in Ohms.

Note
----
fr = 0 may return slightly wrong number for specular boundary conditions.
Not sure why."""
function cmplx_impedance(tr, fr, tc, xk, xm, vf; kwargs...)

    kwargs = Dict(kwargs)
    #Optionally can output penetration/skin depth in meters instead of Ohms
    output_depths = pop!(kwargs, :output_depths, false)

    units = "Ohms"
    if output_depths
        units = "meters"
    end

    verbose = pop!(kwargs, :verbose, false)
    boundary = pop!(kwargs, :boundary, "diffuse")

    gap = pop!(kwargs, :gap, nothing)
    if gap == nothing
        gap = deltar_bcs
    elseif gap == "cos"
        gap = deltar_cos
    end

    dr = gap(tr)
    bcs = pop!(kwargs, :bcs, 1.76)
    delta0 = bcs*kb_eV*tc

    #Calculate the prefactor. Mostly this doesn"t matter since we care about ratios.
    if output_depths
        #Units here are meters
        prefactor = vf*hbar_eV/(xk*delta0)
    else
        #Units here are Ohms
        prefactor = fr*mu_0*vf/xk
    end


    k(x) = cmplx_kernel(tr, fr, x, xk, xm, dr, bcs, verbose)


    if (boundary == "diffuse") || (boundary == "d")
        #Need two integrals to help with singularity removal

        invZint_a(x) = log(x^2+k(x))-2*log(x)
        invZint_b(x) = log(x^2+k(x))

        invZ_a, invZerr = quadgk(invZint_a, 1, Inf)
        invZ_b, invZerr = quadgk(invZint_b, 0, 1)
        #quadgk(x->-2*log(x), 0, 1) = 2

        invZ = invZ_a + invZ_b + 2

        Z = 1.0/invZ

    elseif (boundary == "specular") || (boundary == "s")
        zint(x) = 1.0/(x^2+k(x))

        Z, Zerr = quadgk(zint, -Inf, Inf)

        Z = Z/π^2
    end

    Z *= prefactor*1im

    return Z
end
