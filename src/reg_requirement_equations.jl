
function apply_reg_requirements!(
    problem::PSI.OperationsProblem{T},
    reg⁺_device_names::Vector{String},
    reg⁻_device_names::Vector{String},
    required_reg⁺::Vector{Float64},
    required_reg⁻::Vector{Float64},
    storage_reserve_names::Vector{String},
) where {T <: Union{
    CVaRReserveUnitCommitmentCC,
    StochasticUnitCommitmentCC,
    BasecaseUnitCommitmentCC,
}}
    use_solar_reg = problem.ext["use_solar_reg"]
    use_wind_reserves = problem.ext["use_wind_reserves"]
    use_storage_reserves = false # Hack to overwrite storage reg reserves for now

    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    obj_dict = jump_model.obj_dict
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)

    # These also need registration in the optimization Containers

    reg⁺ = obj_dict[:reg⁺]
    reg⁻ = obj_dict[:reg⁻]
    if use_solar_reg
        reg⁺_S = obj_dict[:reg⁺_S]
        reg⁻_S = obj_dict[:reg⁻_S]
    end
    if use_wind_reserves
        reg⁺_W = obj_dict[:reg⁺_W]
        reg⁻_W = obj_dict[:reg⁻_W]
    end
    if use_slack
        slack_reg⁺ = obj_dict[:slack_reg⁺]
        slack_reg⁻ = obj_dict[:slack_reg⁻]
    end

    if (use_solar_reg && ndims(reg⁺_S) == 2)
        scenarios = 1:size(reg⁺_S)[1]
    else
        scenarios = nothing
    end

    # 1D version
    if isnothing(scenarios)
        # Eq (17) Total reg up
        reg⁺_constraints = JuMP.@constraint(
            jump_model,
            [t in time_steps],
            sum(
                reg⁺[g, t] for g in (
                    use_storage_reserves ?
                    union(reg⁺_device_names, storage_reserve_names) : reg⁺_device_names
                )
            ) +
            (use_solar_reg ? reg⁺_S[t] : 0) +
            (use_wind_reserves ? reg⁺_W[t] : 0) >=
            required_reg⁺[t] - (use_slack ? slack_reg⁺[t] : 0)
        )
        renewable_reg⁺_constraint = JuMP.@constraint(
            jump_model,
            [t in time_steps],
            (use_storage_reserves ? sum(reg⁺[g, t] for g in storage_reserve_names) : 0) +
            (use_solar_reg ? reg⁺_S[t] : 0) +
            (use_wind_reserves ? reg⁺_W[t] : 0) <=
            required_reg⁺[t] .* problem.ext["renewable_reg_prop"]
        )
        # Eq (18) Total reg down
        reg⁻_constraints = JuMP.@constraint(
            jump_model,
            [t in time_steps],
            sum(
                reg⁻[g, t] for g in (
                    use_storage_reserves ?
                    union(reg⁻_device_names, storage_reserve_names) : reg⁻_device_names
                )
            ) +
            (use_solar_reg ? reg⁻_S[t] : 0) +
            (use_wind_reserves ? reg⁻_W[t] : 0) >=
            required_reg⁻[t] - (use_slack ? slack_reg⁻[t] : 0)
        )
        renewable_reg⁻_constraint = JuMP.@constraint(
            jump_model,
            [t in time_steps],
            (use_storage_reserves ? sum(reg⁻[g, t] for g in storage_reserve_names) : 0) +
            (use_solar_reg ? reg⁻_S[t] : 0) +
            (use_wind_reserves ? reg⁻_W[t] : 0) <=
            required_reg⁻[t] .* problem.ext["renewable_reg_prop"]
        )
    else # 2D version
        # Eq (17) Total reg up
        reg⁺_constraints = JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            sum(
                reg⁺[g, t] for g in (
                    use_storage_reserves ?
                    union(reg⁺_device_names, storage_reserve_names) : reg⁺_device_names
                )
            ) +
            (use_solar_reg ? reg⁺_S[j, t] : 0) +
            (use_wind_reserves ? reg⁺_W[t] : 0) >=
            required_reg⁺[t] - (use_slack ? slack_reg⁺[t] : 0)
        )
        renewable_reg⁺_constraint = JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            (use_storage_reserves ? sum(reg⁺[g, t] for g in storage_reserve_names) : 0) +
            (use_solar_reg ? reg⁺_S[j, t] : 0) +
            (use_wind_reserves ? reg⁺_W[t] : 0) <=
            required_reg⁺[t] .* problem.ext["renewable_reg_prop"]
        )
        # Eq (18) Total reg down
        reg⁻_constraints = JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            sum(
                reg⁻[g, t] for g in (
                    use_storage_reserves ?
                    union(reg⁻_device_names, storage_reserve_names) : reg⁻_device_names
                )
            ) +
            (use_solar_reg ? reg⁻_S[j, t] : 0) +
            (use_wind_reserves ? reg⁻_W[t] : 0) >=
            required_reg⁻[t] - (use_slack ? slack_reg⁻[t] : 0)
        )
        renewable_reg⁻_constraint = JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            (use_storage_reserves ? sum(reg⁻[g, t] for g in storage_reserve_names) : 0) +
            (use_solar_reg ? reg⁻_S[j, t] : 0) +
            (use_wind_reserves ? reg⁻_W[t] : 0) <=
            required_reg⁻[t] .* problem.ext["renewable_reg_prop"]
        )
    end

    return
end
