# FEM Benchmark Julia test run script.

# Copyright 2013-2022 Precise Simulation, Ltd.


workspace()

include( "src_julia/fem_poisson.jl" )


# n0    = 512
# nlev  = 1
# nruns = 10
params = readdlm( "testrun_param.txt", Int64 )
n0    = params[1]
nlev  = params[2]
nruns = params[3]


fem_poisson( n0 )
function run_level( n0, ilev, nruns )

    n = n0*2^(ilev-1)
    timings = zeros( Float64, 8 )

    for i = 1:max(3,nruns)
        _, timings_i = fem_poisson( n )
        if i != 1 && i != nruns
            timings = timings + timings_i
        end
    end
    timings = timings/(nruns - 2)

    return timings
end


timings = zeros( Float64, nlev, 9 )
for ilev = 1:nlev
    timings[ilev,1] = n0*2^(ilev-1)
    timings[ilev,2:end] = run_level( n0, ilev, nruns )
end


writedlm( "output/output_julia.txt", timings, " " )


# Profile.init(delay=0.01)
# @profile fem_poisson(256)
# Profile.print()
