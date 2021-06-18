
function apply_wind!(
    problem::PSI.OperationsProblem{T},
    required_reg⁺::Vector{Float64},
    required_reg⁻::Vector{Float64},
    required_spin::Vector{Float64},
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        StochasticUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
    },
}
    use_wind_reserves = problem.ext["use_wind_reserves"]
    use_reg = problem.ext["use_reg"]
    use_spin = problem.ext["use_spin"]

    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    obj_dict = jump_model.obj_dict

    pW = obj_dict[:pW]
    if use_reg && use_wind_reserves
        reg⁺_W = JuMP.@variable(
            jump_model,
            0 <=
            reg⁺_W[t in time_steps] <=
            required_reg⁺[t] .* problem.ext["wind_reg_prop"]
        )
        reg⁻_W = JuMP.@variable(
            jump_model,
            0 <=
            reg⁻_W[t in time_steps] <=
            required_reg⁻[t] .* problem.ext["wind_reg_prop"]
        )
    end
    if use_spin && use_wind_reserves
        spin_W = JuMP.@variable(
            jump_model,
            0 <=
            spin_W[t in time_steps] <=
            required_spin[t] .* problem.ext["wind_spin_prop"]
        )
    end

    # Perfect wind forecast
    total_wind = get_area_total_time_series(
        problem,
        RenewableGen;
        filter = x -> get_prime_mover(x) != PrimeMovers.PVe,
    )

    # Wind constraint
    wind_constraint =
        JuMP.@constraint(jump_model, [t in time_steps], pW[t] <= total_wind[t])

    # Wind reserve holding
    if use_reg && use_wind_reserves
        wind_reserve_dn_constraint =
            JuMP.@constraint(jump_model, [t in time_steps], reg⁻_W[t] <= pW[t])
    end
    if (use_reg || use_spin) && use_wind_reserves
        wind_reserve_up_constraints = JuMP.@constraint(
            jump_model,
            [t in time_steps],
            (use_reg ? reg⁺_W[t] : 0) + (use_spin ? spin_W[t] : 0) <= total_wind[t] - pW[t]
        )
    end

    return
end
