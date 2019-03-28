include("JlBoxModule.jl")
include("Configure_aerosol.jl")
using .Compute:run_simulation_aerosol_adjoint
using DataFrames
using CSV

@time dSOA_mass_drate,dSOA_mass_percentk=run_simulation_aerosol_adjoint(linsolver=:Dense)
df=DataFrame(dSOA_mass_drate)
df2=DataFrame(dSOA_mass_percentk)
CSV.write("/data/jlbox_sensitivity_dSOAdrate_results.csv",df)
CSV.write("/data/jlbox_sensitivity_dSOApercentk_results.csv",df2)
df2
