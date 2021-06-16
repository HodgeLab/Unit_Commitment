
function apply_solar!(
    problem::PSI.OperationsProblem{T},
    required_reg⁺::Vector{Float64},
    required_reg⁻::Vector{Float64},
    required_spin::Vector{Float64}
    ) where T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC}
    use_solar_reg = problem.ext["use_solar_reg"]
    use_solar_spin = problem.ext["use_solar_spin"]
    use_reg = problem.ext["use_reg"]
    solar_scale = problem.ext["solar_scale"]

    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    system = PSI.get_system(problem)
    case_initial_time = PSI.get_initial_time(problem)

    # Populate solar scenarios
    # Scenarios are for entire ERCOT area, but are hacked in under the FarWest name
    area = PSY.get_component(Area, system, "FarWest")
    area_solar_forecast_scenarios = permutedims(
        PSY.get_time_series_values(
            Scenarios,
            area,
            "solar_power";
            start_time = case_initial_time,
        ) .* solar_scale ./ 100,
    )
    scenarios = 1:size(area_solar_forecast_scenarios)[1]

    pS = JuMP.@variable(jump_model, pS[j in scenarios, t in time_steps] >= 0)
    if use_reg && use_solar_reg
        reg⁺_S = JuMP.@variable(
            jump_model,
            0 <= reg⁺_S[j in scenarios, t in time_steps] <=
            required_reg⁺[t] .* problem.ext["solar_reg_prop"]
            )
        reg⁻_S = JuMP.@variable(
            jump_model,
            0 <= reg⁻_S[j in scenarios, t in time_steps] <=
            required_reg⁻[t] .* problem.ext["solar_reg_prop"]
            )
    end
    if use_spin && use_solar_spin
        spin_S = JuMP.@variable(
            jump_model,
            0 <= spin_S[j in scenarios, t in time_steps] <=
            required_spin[t] .* problem.ext["solar_spin_prop"]
        )
    end

    # Eq (23) Solar scenarios
    solar_constraints = JuMP.@constraint(
        jump_model,
        [j in scenarios, t in time_steps],
        pS[j, t] <= area_solar_forecast_scenarios[j, t]
    )

    # Solar reserve holding
    if use_reg && (use_solar_reg || use_solar_spin)
        if use_solar_reg
            solar_reserve_dn_constraint =
            JuMP.@constraint(
                jump_model,
                [j in scenarios, t in time_steps],
                reg⁻_S[j, t] <= pS[j, t]
            )
        end
        solar_reserve_up_constraints =
        JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            (use_solar_reg ? reg⁺_S[j, t] : 0) +
            (use_solar_spin ? spin_S[j, t] : 0) <=
            area_solar_forecast_scenarios[j, t] - pS[j, t]
        )
    end

    return
end

function apply_solar!(
    problem::PSI.OperationsProblem{T},
    required_reg⁺::Vector{Float64},
    required_reg⁻::Vector{Float64},
    required_spin::Vector{Float64}
    ) where T <: BasecaseUnitCommitmentCC
    use_solar_reg = problem.ext["use_solar_reg"]
    use_reg = problem.ext["use_reg"]
    solar_scale = problem.ext["solar_scale"]

    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    system = PSI.get_system(problem)
    case_initial_time = PSI.get_initial_time(problem)

    pS = JuMP.@variable(jump_model, pS[t in time_steps] >= 0)
    if use_reg && use_solar_reg
        reg⁺_S = JuMP.@variable(
            jump_model,
            0 <= reg⁺_S[t in time_steps] <= required_reg⁺[t] .* problem.ext["solar_reg_prop"]
            )
        reg⁻_S = JuMP.@variable(
            jump_model,
            0 <= reg⁻_S[t in time_steps] <= required_reg⁻[t] .* problem.ext["solar_reg_prop"]
            )
    end
    if use_spin && use_solar_spin
        spin_S = JuMP.@variable(
            jump_model,
            0 <= spin_S[t in time_steps] <= required_spin[t] .* problem.ext["solar_spin_prop"]
        )
    end

    # Deterministic solar forecast
    total_solar = get_area_total_time_series(
        problem,
        RenewableGen;
        filter = x -> get_prime_mover(x) == PrimeMovers.PVe && get_available(x),
    ) .* solar_scale

    # Solar power
    solar_constraints = JuMP.@constraint(
        jump_model,
        [t in time_steps],
        pS[t] <= total_solar[t]
    )

    # Solar reserve holding
    if use_reg && use_solar_reg
        solar_reserve_dn_constraint =
        JuMP.@constraint(
            jump_model,
            [t in time_steps],
            reg⁻_S[t] <= pS[t]
        )
    end
    if (use_reg || use_solar_reg) && (use_spin || use_solar_spin)
        solar_reserve_up_constraints =
        JuMP.@constraint(
            jump_model,
            [t in time_steps],
            (use_reg && use_solar_reg ? reg⁺_S[t] : 0) +
            (use_spin && use_solar_spin ? spin_S[t] : 0) <=
            total_solar[t] - pS[t]
        )
    end

    return
end