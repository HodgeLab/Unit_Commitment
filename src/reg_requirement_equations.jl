
function apply_reg_requirements!(problem::PSI.OperationsProblem{T},
    reg⁺_device_names::Vector{String},
    reg⁻_device_names::Vector{String},
    storage_reserve_names::Vector{String}
    ) where T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC, BasecaseUnitCommitmentCC}
    use_solar_reg = problem.ext["use_solar_reg"]
    use_storage_reserves = false # Hack to overwrite storage reg reserves for now

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    case_initial_time = PSI.get_initial_time(problem)
    obj_dict = jump_model.obj_dict
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)

    reg_reserve_up = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "REG_UP")
    reg_reserve_dn =
        PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, "REG_DN")
    required_reg⁺ = get_time_series_values(
        Deterministic,
        reg_reserve_up,
        "requirement";
        start_time = case_initial_time,
    )
    required_reg⁻ = get_time_series_values(
        Deterministic,
        reg_reserve_dn,
        "requirement";
        start_time = case_initial_time,
    )

    reg⁺ = obj_dict[:reg⁺]
    reg⁻ = obj_dict[:reg⁻]
    if use_solar_reg
        reg⁺_S = obj_dict[:reg⁺_S]
        reg⁻_S = obj_dict[:reg⁻_S]
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
                    use_storage_reserves ? union(reg⁺_device_names, storage_reserve_names) :
                    reg⁺_device_names
                )
            ) + (use_solar_reg ? reg⁺_S[t] : 0) >=
            required_reg⁺[t] - (use_slack ? slack_reg⁺[t] : 0)
        )
        # Eq (18) Total reg down
        reg⁻_constraints = JuMP.@constraint(
            jump_model,
            [t in time_steps],
            sum(
                reg⁻[g, t] for g in (
                    use_storage_reserves ? union(reg⁻_device_names, storage_reserve_names) :
                    reg⁻_device_names
                )
            ) + (use_solar_reg ? reg⁻_S[t] : 0) >=
            required_reg⁻[t] - (use_slack ? slack_reg⁻[t] : 0)
        )
    else # 2D version
        # Eq (17) Total reg up
        reg⁺_constraints = JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            sum(
                reg⁺[g, t] for g in (
                    use_storage_reserves ? union(reg⁺_device_names, storage_reserve_names) :
                    reg⁺_device_names
                )
            ) + (use_solar_reg ? reg⁺_S[j, t] : 0)  >=
            required_reg⁺[t] - (use_slack ? slack_reg⁺[t] : 0)
        )
        # Eq (18) Total reg down
        reg⁻_constraints = JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            sum(
                reg⁻[g, t] for g in (
                    use_storage_reserves ? union(reg⁻_device_names, storage_reserve_names) :
                    reg⁻_device_names
                )
            ) + (use_solar_reg ? reg⁻_S[j, t] : 0) >=
            required_reg⁻[t] - (use_slack ? slack_reg⁻[t] : 0)
        )
    end

    return
end
