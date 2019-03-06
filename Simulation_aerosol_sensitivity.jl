include("JlBoxModule.jl")
include("Configure_aerosol.jl")
using .Compute:run_simulation_aerosol_DDM,run_simulation_aerosol_adjoint
using DataFrames
using CSV

@time dSOA_mass_drate=run_simulation_aerosol_adjoint(linsolver=:Dense)
df=DataFrame(dSOA_mass_drate)
CSV.write("/data/jlbox_sensitivity_dSOAdrate_results.csv",df)
df
