"""Calculate the Fermi Function given some energy and some temperature.

Parameters
----------
en : float
    Energy relative to the fermi energy (E-Ef).

kbt : Temperature in same units as Energy (using kB).

Returns
-------
ffun : float
    The Fermi Function at en and temp."""
function fermi_fun(en, kbt)
    if en == 0
        ffun = 0.5
    elseif kbt == 0
        if en < 0
            ffun = 1.0
        elseif en > 0
            ffun = 0.0
        else
            throw(DomainError("Energy must be a real number."))
        end
    elseif kbt > 0
        #Using a tanh here instead of the more traditional formulation because
        #the function is better behaved for very large or small values of en/kbt
        ffun = 0.5*(1-tanh(0.5*en/kbt))

    elseif kbt < 0
        throw(DomainError("Temperature must be >= 0."))
    end

    return ffun
end
