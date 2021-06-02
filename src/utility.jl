
function get_area_total_time_series(
    problem::Union{PSI.OperationsProblem{CVaRPowerUnitCommitmentCC}, PSI.OperationsProblem{CVaRReserveUnitCommitmentCC}},
    type::DataType; 
    filter = nothing
    )
    system = PSI.get_system(problem)
    case_initial_time = PSI.get_initial_time(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)

    total = zeros(length(time_steps))

    if isnothing(filter)
        iter = get_components(type, system)
    else
        iter = get_components(type, system, filter)
    end

    for l in iter
        total .+= get_time_series_values(
            Deterministic,
            l,
            "max_active_power";
            start_time = case_initial_time,
        )
    end

    return total
end

# Unconventional route. To be cleaned later.
function PSI.write_to_CSV(
    problem::Union{PSI.OperationsProblem{CVaRPowerUnitCommitmentCC}, PSI.OperationsProblem{CVaRReserveUnitCommitmentCC}},
    output_path::String;
    time = nothing
)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    exclusions = [:λ, :β] # PWL chunks, expensive to export and useless
    for (k, v) in jump_model.obj_dict
        if !(k in exclusions)
            df = PSI.axis_array_to_dataframe(v, [k])
            file_name = joinpath(output_path, string(k) * ".csv")
            CSV.write(file_name, df)
        end
    end

    _write_summary_stats(problem, output_path, time)
end

function _write_summary_stats(
    problem::Union{PSI.OperationsProblem{CVaRPowerUnitCommitmentCC}, PSI.OperationsProblem{CVaRReserveUnitCommitmentCC}},
    output_path::String,
    solvetime::Union{Nothing, Float64}
)
    optimization_container = PSI.get_optimization_container(problem)
    system = PSI.get_system(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)
    use_storage = problem.ext["use_storage"]
    use_reg = problem.ext["use_reg"]
    use_spin = problem.ext["use_spin"]
    C_RR = problem.ext["C_RR"]
    α = problem.ext["α"]
    C_res_penalty = problem.ext["C_res_penalty"]
    C_ener_penalty = problem.ext["C_ener_penalty"]

    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    no_load_cost = Dict(
        g => get_no_load(get_operation_cost(get_component(ThermalMultiStart, system, g))) for g in thermal_gen_names
    )
    shutdown_cost = Dict(
        g => get_shut_down(get_operation_cost(get_component(ThermalMultiStart, system, g))) for g in thermal_gen_names
    )
    startup_cost = Dict(
        g => get_start_up(get_operation_cost(get_component(ThermalMultiStart, system, g))) for g in thermal_gen_names
    )

    obj_dict = PSI.get_jump_model(optimization_container).obj_dict
    ug = PSI.axis_array_to_dataframe(obj_dict[:ug], [:ug])
    wg = PSI.axis_array_to_dataframe(obj_dict[:wg], [:wg])
    z = PSI.axis_array_to_dataframe(obj_dict[:z], [:z])
    δ_sg = obj_dict[:δ_sg]
    Cg = PSI.axis_array_to_dataframe(JuMP.value.(optimization_container.expressions[:Cg]))

    output = Dict(
        "Solve time (s)" => solvetime,
        "L_SUPP" => problem.ext["L_SUPP"],
        "C_RR" => C_RR,
        "C_res_penalty" => C_res_penalty,
        "C_ener_penalty" => C_ener_penalty,
        "alpha" => α,
        "Hot start cost" => JuMP.value(sum(startup_cost[g][:hot] * δ_sg[g, :hot, t] for g in thermal_gen_names, t in time_steps)),
        "Warm start cost" => JuMP.value(sum(startup_cost[g][:warm] * δ_sg[g, :warm, t] for g in thermal_gen_names, t in time_steps)),
        "Cold start cost" => JuMP.value(sum(startup_cost[g][:cold] * δ_sg[g, :cold, t] for g in thermal_gen_names, t in time_steps)),
        "No-load cost" => sum(sum(ug[!, n] .* no_load_cost[n]) for n in names(ug)),
        "Variable cost" => sum(sum(eachcol(Cg))),
        "Shut-down cost" => sum(sum(wg[!, n] .* shutdown_cost[n]) for n in names(wg)),
        "CVaR cost" =>
        C_RR * (PSI._jump_value(obj_dict[:β]) + 1 / (nrow(z) * (1 - α)) * sum(z[!, :z]))
    )
    output["Start-up cost"] =
        output["Hot start cost"] +
        output["Warm start cost"] +
        output["Cold start cost"]
    output["Total cost"] =
        output["No-load cost"] +
        output["Variable cost"] +
        output["Start-up cost"] +
        output["Shut-down cost"] +
        output["CVaR cost"]

    if use_slack
        if use_reg
            slack_reg⁺ = PSI.axis_array_to_dataframe(obj_dict[:slack_reg⁺], [:slack_reg⁺])
            slack_reg⁻ = PSI.axis_array_to_dataframe(obj_dict[:slack_reg⁻], [:slack_reg⁻])
        end
        if use_spin
            slack_spin = PSI.axis_array_to_dataframe(obj_dict[:slack_spin], [:slack_spin])
        end
        slack_energy⁺ = PSI.axis_array_to_dataframe(obj_dict[:slack_energy⁺], [:slack_energy⁺])
        slack_energy⁻ = PSI.axis_array_to_dataframe(obj_dict[:slack_energy⁻], [:slack_energy⁻])
    end

    output["Penalty cost unserved reg up"] = (use_slack && use_reg) ? C_res_penalty * sum(slack_reg⁺[!, :slack_reg⁺]) : nothing
    output["Penalty cost unserved reg down"] = (use_slack && use_reg) ? C_res_penalty * sum(slack_reg⁻[!, :slack_reg⁻]) : nothing
    output["Penalty cost unserved spin"] = (use_slack && use_spin) ? C_res_penalty * sum(slack_spin[!, :slack_spin]) : nothing
    output["Penalty cost unserved load"] = use_slack ? C_ener_penalty * sum(slack_energy⁺[!, :slack_energy⁺]) : nothing
    output["Penalty cost overgeneration"] = use_slack ? C_ener_penalty * sum(slack_energy⁻[!, :slack_energy⁻]) : nothing

    if use_storage
        pb_in = PSI.axis_array_to_dataframe(obj_dict[:pb_in], [:pb_in])
        pb_out = PSI.axis_array_to_dataframe(obj_dict[:pb_out], [:pb_out])
    end

    output["Total cost with penalties"] = output["Total cost"] +
        ((use_slack && use_reg) ? output["Penalty cost unserved reg up"] + output["Penalty cost unserved reg down"] : 0) +
        ((use_slack && use_spin) ? output["Penalty cost unserved spin"] : 0) +
        (use_slack ? output["Penalty cost unserved load"] : 0) +
        (use_slack ? output["Penalty cost overgeneration"] : 0)

    CSV.write(joinpath(output_path, "Summary_stats.csv"), output)
end
