
function apply_solar!(problem::PSI.OperationsProblem{T}
    ) where T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC}
    use_solar_reserves = problem.ext["use_solar_reserves"]
    use_reg = problem.ext["use_reg"]

    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    system = PSI.get_system(problem)
    case_initial_time = PSI.get_initial_time(problem)

    # Populate solar scenarios
    area = PSY.get_component(Area, system, "1")
    area_solar_forecast_scenarios = permutedims(
        PSY.get_time_series_values(
            Scenarios,
            area,
            "solar_power";
            start_time = case_initial_time,
        ) ./ 100,
    )
    scenarios = 1:size(area_solar_forecast_scenarios)[1]

    pS = JuMP.@variable(jump_model, pS[j in scenarios, t in time_steps] >= 0)
    if use_reg && use_solar_reserves
        reg⁺_S = JuMP.@variable(jump_model, reg⁺_S[j in scenarios, t in time_steps] >= 0)
        reg⁻_S = JuMP.@variable(jump_model, reg⁻_S[j in scenarios, t in time_steps] >= 0)
    end

    # Eq (23) Solar scenarios
    solar_constraints = JuMP.@constraint(
        jump_model,
        [j in scenarios, t in time_steps],
        pS[j, t] <= area_solar_forecast_scenarios[j, t]
    )

    # Solar reserve holding
    if use_reg && use_solar_reserves
        solar_reserve_dn_constraint =
        JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            reg⁻_S[j, t] <= pS[j, t]
        )
        solar_reserve_up_constraints =
        JuMP.@constraint(
            jump_model,
            [j in scenarios, t in time_steps],
            reg⁺_S[j, t] <=
            area_solar_forecast_scenarios[j, t] - pS[j, t]
        )
    end

    return
end

function apply_solar!(problem::PSI.OperationsProblem{T}
    ) where T <: BasecaseUnitCommitmentCC
    use_solar_reserves = problem.ext["use_solar_reserves"]
    use_reg = problem.ext["use_reg"]

    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    system = PSI.get_system(problem)

    pS = JuMP.@variable(jump_model, pS[t in time_steps] >= 0)
    if use_reg && use_solar_reserves
        reg⁺_S = JuMP.@variable(jump_model, reg⁺_S[t in time_steps] >= 0)
        reg⁻_S = JuMP.@variable(jump_model, reg⁻_S[t in time_steps] >= 0)
    end

    # Deterministic solar forecast
    total_solar = get_area_total_time_series(
        problem,
        RenewableGen;
        filter = x -> get_prime_mover(x) == PrimeMovers.PVe && get_available(x),
    )

    # Eq (23) Solar power
    solar_constraints = JuMP.@constraint(
        jump_model,
        [t in time_steps],
        pS[t] <= total_solar[t]
    )

    # Solar reserve holding
    if use_reg && use_solar_reserves
        solar_reserve_dn_constraint =
        JuMP.@constraint(
            jump_model,
            [t in time_steps],
            reg⁻_S[t] <= pS[t]
        )
        solar_reserve_up_constraints =
        JuMP.@constraint(
            jump_model,
            [t in time_steps],
            reg⁺_S[t] <=
            total_solar[t] - pS[t]
        )
    end

    return
end