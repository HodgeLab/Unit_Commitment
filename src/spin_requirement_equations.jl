
function apply_spin_requirements!(
    problem::PSI.OperationsProblem{T},
    spin_device_names::Vector{String},
    required_spin::Vector{Float64},
    storage_reserve_names::Vector{String},
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        StochasticUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
        HourAheadUnitCommitmentCC
    },
}
    use_solar_spin = problem.ext["use_solar_spin"]
    use_storage_reserves = problem.ext["use_storage_reserves"]
    use_wind_reserves = problem.ext["use_wind_reserves"]

    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    obj_dict = jump_model.obj_dict
    use_supp =
        :total_supp in keys(obj_dict) ||
        :total_supp in keys(optimization_container.expressions)
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)

    # Needs registration into the optimization container
    spin = obj_dict[:spin]
    if use_solar_spin
        spin_S = obj_dict[:spin_S]
    end
    if use_wind_reserves
        spin_W = obj_dict[:spin_W]
    end
    if use_slack
        slack_spin = obj_dict[:slack_spin]
    end

    if use_supp
        if problem.ext["supp_type"] == "generic"
            total_supp = obj_dict[:total_supp]
        else
            total_supp = optimization_container.expressions[:total_supp]
        end
    end

    if use_supp
        scenarios = 1:size(total_supp)[1]
    elseif (use_solar_spin && ndims(spin_S) == 2)
        scenarios = 1:size(spin_S)[1]
    else
        scenarios = nothing
    end

    # 1D version
    if isnothing(scenarios)
        # Eq (19) Total spin
        spin_constraints = JuMP.@constraint(
            jump_model,
            [t in time_steps],
            sum(
                spin[g, t] for g in (
                    use_storage_reserves ?
                    union(spin_device_names, storage_reserve_names) : spin_device_names
                )
            ) +
            (use_solar_spin ? spin_S[t] : 0) +
            (use_wind_reserves ? spin_W[t] : 0) >=
            required_spin[t] - (use_slack ? slack_spin[t] : 0)
        )
        renewable_spin_constraint = JuMP.@constraint(
            jump_model,
            [t in time_steps],
            (use_storage_reserves ? sum(spin[g, t] for g in storage_reserve_names) : 0) +
            (use_solar_spin ? spin_S[t] : 0) +
            (use_wind_reserves ? spin_W[t] : 0) <=
            required_spin[t] .* problem.ext["renewable_spin_prop"]
        )
    else # 2D version
        # Eq (19) Total spin
        spin_constraints = JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            sum(
                spin[g, t] for g in (
                    use_storage_reserves ?
                    union(spin_device_names, storage_reserve_names) : spin_device_names
                )
            ) +
            (use_solar_spin ? spin_S[j, t] : 0) +
            (use_wind_reserves ? spin_W[t] : 0) >=
            required_spin[t] - (use_supp ? total_supp[j, t] : 0) -
            (use_slack ? slack_spin[t] : 0)
        )
        renewable_spin_constraint = JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            (use_storage_reserves ? sum(spin[g, t] for g in storage_reserve_names) : 0) +
            (use_solar_spin ? spin_S[j, t] : 0) +
            (use_wind_reserves ? spin_W[t] : 0) <=
            required_spin[t] .* problem.ext["renewable_spin_prop"]
        )
    end

    return
end
