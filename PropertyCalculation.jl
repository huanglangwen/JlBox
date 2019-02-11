using PyCall
unshift!(PyVector(pyimport("sys")["path"]),"../UManSysProp_public")
@pyimport umansysprop.boiling_points as boiling_points
@pyimport umansysprop.vapour_pressures as vapour_pressures
@pyimport umansysprop.critical_properties as critical_properties
@pyimport umansysprop.liquid_densities as liquid_densities
@pyimport umansysprop.partition_models as partition_models
@pyimport umansysprop.groups as groups
@pyimport umansysprop.activity_coefficient_models_dev as aiomfac
@pyimport umansysprop.forms as forms #need forms.CoreAbundanceField (class)
CoreAbundanceField=forms.CoreAbundanceField

