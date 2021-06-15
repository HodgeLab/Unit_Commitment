
function apply_spin_requirements!(problem::PSI.OperationsProblem{T},
    spin_device_names::Vector{String},
    storage_reserve_names::Vector{String}
    ) where T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC, BasecaseUnitCommitmentCC}
    use_solar_spin = problem.ext["use_solar_spin"]
    use_storage_reserves = problem.ext["use_storage_reserves"]

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    case_initial_time = PSI.get_initial_time(problem)
    obj_dict = jump_model.obj_dict
    use_supp = :supp in keys(obj_dict)
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)

    spin_reserve = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "SPIN")
    required_spin = get_time_series_values(
        Deterministic,
        spin_reserve,
        "requirement";
        start_time = case_initial_time,
    )

    spin = obj_dict[:spin]
    if use_solar_spin
        spin_S = obj_dict[:spin_S]
    end
    if use_slack
        slack_spin = obj_dict[:slack_spin]
    end

    if use_supp
        supp = obj_dict[:supp]
    end

    if use_supp || (use_solar_spin && ndims(spin_S) == 2)
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
                    use_storage_reserves ? union(spin_device_names, storage_reserve_names) :
                    spin_device_names
                )
            ) + (use_solar_spin ? spin_S[t] : 0) >=
            required_spin[t] - (use_slack ? slack_spin[t] : 0)
        )
    else # 2D version
        # Eq (19) Total spin
        spin_constraints = JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            sum(
                spin[g, t] for g in (
                    use_storage_reserves ? union(spin_device_names, storage_reserve_names) :
                    spin_device_names
                )
            ) + (use_solar_spin ? spin_S[j, t] : 0) >=
            required_spin[t] - (use_supp ? supp[j, t] : 0) - (use_slack ? slack_spin[t] : 0)
        )
    end

    return
end
