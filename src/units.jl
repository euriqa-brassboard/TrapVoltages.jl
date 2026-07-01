#

module Units

struct TrapUnits
    l_unit::Float64
    V_unit::Float64
    l_unit_um::Float64
    V_unit_uV::Float64

    TrapUnits(l_unit, V_unit) = new(l_unit, V_unit, l_unit * 1e6, V_unit * 1e6)
end

const q_e = 1.60217663e-19 # C
const ε_0 = 8.8541878128e-12

# Unit such that electric potential that creates 1MHz trapping frequency
# for a particular mass has the form X^2/2,
# and the electric potential between two ions is 1/r.
function MHz_unit(mass, q=1)
    MHz = 1e6
    q = q * q_e
    A = mass * (2π * MHz)^2 / q
    B = q / (4π * ε_0)

    return TrapUnits(cbrt(B / A), cbrt(A * B^2))
end

const N_A = 6.02214076e23

public SI, Yb171, Yb172, Yb174, Ba133, Ba137, Ba138, euriqa

const SI = TrapUnits(1, 1)
const Yb171 = MHz_unit(170.9363315e-3 / N_A)
const Yb172 = MHz_unit(171.936386654e-3 / N_A)
const Yb174 = MHz_unit(173.938867546e-3 / N_A)
const Ba133 = MHz_unit(132.9060074e-3 / N_A)
const Ba137 = MHz_unit(136.90582721e-3 / N_A)
const Ba138 = MHz_unit(137.90524706e-3 / N_A)

const euriqa = Yb171

end
