```@meta
CurrentModule = QuadPol
```

# QuadPol.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/JuliaHCI/QuadPol.jl)
[![Build Status](https://github.com/JuliaHCI/QuadPol.jl/workflows/CI/badge.svg)](https://github.com/JuliaHCI/QuadPol.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaHCI/QuadPol.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaHCI/QuadPol.jl)
[![License](https://img.shields.io/github/license/JuliaHCI/QuadPol.jl?color=yellow)](https://github.com/JuliaHCI/QuadPol.jl/blob/main/LICENSE)

**Primary author:** Miles Lucas <mdlucas@hawaii.edu>

Quadrant polarization statistics based on [Schmid et al., 2021](https://www.aanda.org/articles/aa/abs/2021/11/aa40405-21/aa40405-21.html).

## Installation


## Getting Started

```julia-repl
julia> Q, U = # load your Stokes data as arrays
```

position angle is defined as the near-side minor axis in degrees clockwise from North

```julia-repl
julia> pa = -118 # degrees
-118

julia> using QuadPol # needed for quadpol function

julia> stats = quadpol(Q, U; pa=pa);
Relative quadrant polarization statistics
Q000/Qphi =     -0.006138270764788408   U045/Qphi =     -0.03825566505358086
Q090/Qphi =     0.47013684691873975     U135/Qphi =     0.24783192367196974
Q180/Qphi =     -0.020413603056487466   U225/Qphi =     -0.14450180863634055
Q270/Qphi =     0.3894757166335373      U315/Qphi =     -0.08027138677188023

Quadrant sum statistics
ΣQxxx/Qphi =    0.8330606897310011      ΣUxxx/Qphi =    -0.015196936789831882
Σ|Qxxx|/Qphi =  0.886164437373553       Σ|Uxxx|/Qphi =  0.5108607841337713

Left-right asymmetry deviations
Δ090/270 =      9.38342850084427%       Δ045/315 =      -35.448212936377%
                                Δ135/225 =      26.337300753540255%

Back-front parameter ratios
Λ000/180 =      0.3006951172609217      Λ045/135 =      0.1543613287859398
                                Λ315/225 =      0.5555043741625038

Special back-front parameter ratios
Λa = (|Q090|+|Q270|)/(2|Q180|)  =               21.054895629487888
Λb = (|Q090|+|Q270|)/(|U135|+|U225|) =          2.1910238472096557
```


the stats are just a `NamedTuple` and can be exported to a table easily

```julia-repl
julia> using DataFrames

julia> table = DataFrame([stats])
1×22 DataFrame
 Row │ Q000      Q090     Q180     Q270     U045      U135     U225      U315      Qd       Ud        Qd_abs   Ud_abs   Qphi     Uphi      delta_090_270  delta_045_ ⋯
     │ Float64   Float64  Float64  Float64  Float64   Float64  Float64   Float64   Float64  Float64   Float64  Float64  Float64  Float64   Float64        Float64    ⋯
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ -212.889  16305.4  -707.99  13507.9  -1326.79  8595.37  -5011.65  -2783.99  28892.4  -527.064  30734.2  17717.8  34682.3  -1405.09      0.0938343      -0.354 ⋯
                                                                                                                                                     7 columns omitted
```

### With Errors

```julia-repl
julia> Qerr, Uerr = # optionally load Qerr and Uerr arrays

julia> stats = quadpol(Q, Qerr, U, Uerr; pa=pa);
Relative quadrant polarization statistics
Q000/Qphi =     -0.0061 ± 0.004 U045/Qphi =     -0.0383 ± 0.0034
Q090/Qphi =     0.4701 ± 0.0047 U135/Qphi =     0.2478 ± 0.0045
Q180/Qphi =     -0.0204 ± 0.0046        U225/Qphi =     -0.1445 ± 0.0039
Q270/Qphi =     0.3895 ± 0.004  U315/Qphi =     -0.0803 ± 0.0036

Quadrant sum statistics
ΣQxxx/Qphi =    0.8331 ± 0.0096 ΣUxxx/Qphi =    -0.0152 ± 0.0074
Σ|Qxxx|/Qphi =  0.8862 ± 0.0099 Σ|Uxxx|/Qphi =  0.5109 ± 0.0083

Left-right asymmetry deviations
Δ090/270 =      9.38 ± 0.48%    Δ045/315 =      -35.4 ± 4.4%
                                Δ135/225 =      26.3 ± 1.4%

Back-front parameter ratios
Λ000/180 =      0.3 ± 0.21      Λ045/135 =      0.154 ± 0.014
                                Λ315/225 =      0.556 ± 0.028

Special back-front parameter ratios
Λa = (|Q090|+|Q270|)/(2|Q180|)  =               21.1 ± 4.7
Λb = (|Q090|+|Q270|)/(|U135|+|U225|) =          2.191 ± 0.033
```

```@index
```

```@autodocs
Modules = [QuadPol]
```

## Citation

Cite Schmid+2021

## Contributing

TODO