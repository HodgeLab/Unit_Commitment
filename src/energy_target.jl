struct EnergyTargetFF <: PSI.AbstractAffectFeedForward
    variable_source_problem::Symbol
    affected_variables::Vector{Symbol}
    target_period::Int
    penalty_cost::Float64
    cache::Union{Nothing, Type{<:PSI.AbstractCache}}
    function EnergyTargetFF(
        variable_source_problem::AbstractString,
        affected_variables::Vector{<:AbstractString},
        target_period::Int,
        penalty_cost::Float64,
        cache::Union{Nothing, Type{<:PSI.AbstractCache}},
    )
        new(
            Symbol(variable_source_problem),
            Symbol.(affected_variables),
            target_period,
            penalty_cost,
            cache,
        )
    end
end

function EnergyTargetFF(;
    variable_source_problem,
    affected_variables,
    target_period,
    penalty_cost,
)
    return EnergyTargetFF(
        variable_source_problem,
        affected_variables,
        target_period,
        penalty_cost,
        nothing,
    )
end

PSI.get_variable_source_problem(p::EnergyTargetFF) = p.variable_source_problem

function energy_target_ff(
    optimization_container::PSI.OptimizationContainer,
    cons_name::Symbol,
    param_reference::PSI.UpdateRef,
    var_name::Tuple{Symbol, Symbol},
    target_period::Int,
    penalty_cost::Float64,
)
    time_steps = PSI.model_time_steps(optimization_container)
    variable = PSI.get_variable(optimization_container, var_name[1])
    varslack = PSI.get_variable(optimization_container, var_name[2])
    axes = JuMP.axes(variable)
    set_name = axes[1]

    @assert axes[2] == time_steps
    container_ub =
        PSI.add_param_container!(optimization_container, param_reference, set_name)
    param_ub = PSI.get_parameter_array(container_ub)
    multiplier_ub = PSI.get_multiplier_array(container_ub)
    con_ub = PSI.add_cons_container!(optimization_container, cons_name, set_name)

    for name in axes[1]
        value = JuMP.upper_bound(variable[name, 1])
        param_ub[name] = PSI.add_parameter(optimization_container.JuMPmodel, value)
        # default set to 1.0, as this implementation doesn't use multiplier
        con_ub[name] = JuMP.@constraint(
            optimization_container.JuMPmodel,
            variable[name, target_period] == param_ub[name] # - varslack[name, target_period]
        )
        # PSI.linear_gen_cost!(
        #      optimization_container,
        #      var_name[2],
        #      name,
        #      penalty_cost,
        #      target_period,
        # )
    end
end

function PSI.feedforward!(
    optimization_container::PSI.OptimizationContainer,
    devices::IS.FlattenIteratorWrapper{T},
    ::PSI.DeviceModel{T, D},
    ff_model::EnergyTargetFF,
) where {T <: PSY.StaticInjection, D <: PSI.AbstractDeviceFormulation}
    PSI.add_variables!(optimization_container, PSI.EnergyShortageVariable, devices, D())
    for prefix in PSI.get_affected_variables(ff_model)
        var_name = PSI.make_variable_name(prefix, T)
        varslack_name = PSI.make_variable_name(PSI.ENERGY_SHORTAGE, T)
        parameter_ref = PSI.UpdateRef{JuMP.VariableRef}(var_name)
        energy_target_ff(
            optimization_container,
            PSI.make_constraint_name("FF_energy_target", T),
            parameter_ref,
            (var_name, varslack_name),
            ff_model.target_period,
            ff_model.penalty_cost,
        )
    end
end
