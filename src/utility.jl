
function get_area_total_time_series(
    problem::PSI.OperationsProblem{T},
    type::DataType;
    filter = nothing,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
        StochasticUnitCommitmentCC,
    },
}
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

function get_thermal_generator_power_dataframe(
    problem::PSI.OperationsProblem{T},
    time_steps::UnitRange{Int64},
    scenario,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        StochasticUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
    },
}
    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)

    # Get power above min
    pg = jump_model.obj_dict[:pg]
    if length(size(pg)) == 2
        pg_dataframe =
            PSI.axis_array_to_dataframe(jump_model.obj_dict[:pg], [:pg])[time_steps, :]
    elseif length(size(pg)) == 3
        pg_dataframe = _scenario_in_3D_array_to_dataframe(pg, scenario, time_steps)
    else
        throw(ErrorException("Only support pg of 2 or 3 dimensions"))
    end
    # Add the min level
    Pg = PSI.axis_array_to_dataframe(jump_model.obj_dict[:ug], [:ug])[time_steps, :]
    for n in names(Pg)
        Pg[!, n] .*=
            get_active_power_limits(get_component(ThermalMultiStart, system, n)).min
    end
    pg_dataframe .+= Pg
    return pg_dataframe
end

# Unconventional route. To be cleaned later.
function PSI.write_to_CSV(
    problem::PSI.OperationsProblem{T},
    data_path::String,
    output_path::String;
    time = nothing,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
        StochasticUnitCommitmentCC,
    },
}
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    exclusions = [:λ, :β] # PWL chunks, expensive to export and useless
    for (k, v) in jump_model.obj_dict
        print("writing $k\n")
        if !(k in exclusions)
            df = PSI.axis_array_to_dataframe(v, [k])
            file_name = joinpath(output_path, string(k) * ".csv")
            CSV.write(file_name, df)
        end
    end

    _write_summary_stats(problem, output_path, time)
end

function _write_summary_stats(
    problem::PSI.OperationsProblem{T},
    output_path::String,
    solvetime::Union{Nothing, Float64},
) where {T <: CVaRReserveUnitCommitmentCC}
    optimization_container = PSI.get_optimization_container(problem)
    system = PSI.get_system(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)
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
    Cg = JuMP.value.(optimization_container.expressions[:Cg]).data

    output = Dict(
        "Solve time (s)" => solvetime,
        "C_RR" => C_RR,
        "C_res_penalty" => C_res_penalty,
        "C_ener_penalty" => C_ener_penalty,
        "alpha" => α,
        "Hot start cost" => JuMP.value(
            sum(
                startup_cost[g][:hot] * δ_sg[g, :hot, t] for g in thermal_gen_names,
                t in time_steps
            ),
        ),
        "Warm start cost" => JuMP.value(
            sum(
                startup_cost[g][:warm] * δ_sg[g, :warm, t] for
                g in thermal_gen_names, t in time_steps
            ),
        ),
        "Cold start cost" => JuMP.value(
            sum(
                startup_cost[g][:cold] * δ_sg[g, :cold, t] for
                g in thermal_gen_names, t in time_steps
            ),
        ),
        "No-load cost" => sum(sum(ug[!, n] .* no_load_cost[n]) for n in names(ug)),
        "Variable cost" => sum(Cg) / (ndims(Cg) == 3 ? size(Cg)[2] : 1),
        "Shut-down cost" => sum(sum(wg[!, n] .* shutdown_cost[n]) for n in names(wg)),
        "CVaR cost" =>
            C_RR *
            (PSI._jump_value(obj_dict[:β]) + 1 / (nrow(z) * (1 - α)) * sum(z[!, :z])),
    )
    output["Start-up cost"] =
        output["Hot start cost"] + output["Warm start cost"] + output["Cold start cost"]
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
        slack_energy⁺ =
            PSI.axis_array_to_dataframe(obj_dict[:slack_energy⁺], [:slack_energy⁺])
        slack_energy⁻ =
            PSI.axis_array_to_dataframe(obj_dict[:slack_energy⁻], [:slack_energy⁻])
    end

    if use_slack
        if use_reg
            output["Penalty cost unserved reg up"] =
                C_res_penalty * sum(slack_reg⁺[!, :slack_reg⁺])
            output["Penalty cost unserved reg down"] =
                C_res_penalty * sum(slack_reg⁻[!, :slack_reg⁻])
        end
        if use_spin
            output["Penalty cost unserved spin"] =
                C_res_penalty * sum(slack_spin[!, :slack_spin])
        end
        output["Penalty cost unserved load"] =
            C_ener_penalty * sum(slack_energy⁺[!, :slack_energy⁺])
        output["Penalty cost overgeneration"] =
            C_ener_penalty * sum(slack_energy⁻[!, :slack_energy⁻])

        output["Total cost with penalties"] =
            output["Total cost"] +
            (
                use_reg ?
                output["Penalty cost unserved reg up"] +
                output["Penalty cost unserved reg down"] : 0
            ) +
            (use_spin ? output["Penalty cost unserved spin"] : 0) +
            output["Penalty cost unserved load"] +
            output["Penalty cost overgeneration"]
    end

    CSV.write(joinpath(output_path, "Summary_stats.csv"), output)
end

function _write_summary_stats(
    problem::PSI.OperationsProblem{T},
    output_path::String,
    solvetime::Union{Nothing, Float64},
) where {T <: Union{BasecaseUnitCommitmentCC, StochasticUnitCommitmentCC}}
    optimization_container = PSI.get_optimization_container(problem)
    system = PSI.get_system(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)
    use_reg = problem.ext["use_reg"]
    use_spin = problem.ext["use_spin"]
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
    δ_sg = obj_dict[:δ_sg]
    Cg = JuMP.value.(optimization_container.expressions[:Cg]).data

    output = Dict(
        "Solve time (s)" => solvetime,
        "C_res_penalty" => C_res_penalty,
        "C_ener_penalty" => C_ener_penalty,
        "Hot start cost" => JuMP.value(
            sum(
                startup_cost[g][:hot] * δ_sg[g, :hot, t] for g in thermal_gen_names,
                t in time_steps
            ),
        ),
        "Warm start cost" => JuMP.value(
            sum(
                startup_cost[g][:warm] * δ_sg[g, :warm, t] for
                g in thermal_gen_names, t in time_steps
            ),
        ),
        "Cold start cost" => JuMP.value(
            sum(
                startup_cost[g][:cold] * δ_sg[g, :cold, t] for
                g in thermal_gen_names, t in time_steps
            ),
        ),
        "No-load cost" => sum(sum(ug[!, n] .* no_load_cost[n]) for n in names(ug)),
        "Variable cost" => sum(Cg) / (ndims(Cg) == 3 ? size(Cg)[2] : 1),
        "Shut-down cost" => sum(sum(wg[!, n] .* shutdown_cost[n]) for n in names(wg)),
    )
    output["Start-up cost"] =
        output["Hot start cost"] + output["Warm start cost"] + output["Cold start cost"]
    output["Total cost"] =
        output["No-load cost"] +
        output["Variable cost"] +
        output["Start-up cost"] +
        output["Shut-down cost"]

    if use_slack
        if use_reg
            slack_reg⁺ = PSI.axis_array_to_dataframe(obj_dict[:slack_reg⁺], [:slack_reg⁺])
            slack_reg⁻ = PSI.axis_array_to_dataframe(obj_dict[:slack_reg⁻], [:slack_reg⁻])
        end
        if use_spin
            slack_spin = PSI.axis_array_to_dataframe(obj_dict[:slack_spin], [:slack_spin])
        end
        slack_energy⁺ =
            PSI.axis_array_to_dataframe(obj_dict[:slack_energy⁺], [:slack_energy⁺])
        slack_energy⁻ =
            PSI.axis_array_to_dataframe(obj_dict[:slack_energy⁻], [:slack_energy⁻])
    end

    if use_slack
        if use_reg
            output["Penalty cost unserved reg up"] =
                C_res_penalty * sum(slack_reg⁺[!, :slack_reg⁺])
            output["Penalty cost unserved reg down"] =
                C_res_penalty * sum(slack_reg⁻[!, :slack_reg⁻])
        end
        if use_spin
            output["Penalty cost unserved spin"] =
                C_res_penalty * sum(slack_spin[!, :slack_spin])
        end
        output["Penalty cost unserved load"] =
            C_ener_penalty * sum(slack_energy⁺[!, :slack_energy⁺])
        output["Penalty cost overgeneration"] =
            C_ener_penalty * sum(slack_energy⁻[!, :slack_energy⁻])

        output["Total cost with penalties"] =
            output["Total cost"] +
            (
                use_reg ?
                output["Penalty cost unserved reg up"] +
                output["Penalty cost unserved reg down"] : 0
            ) +
            (use_spin ? output["Penalty cost unserved spin"] : 0) +
            output["Penalty cost unserved load"] +
            output["Penalty cost overgeneration"]
    end

    CSV.write(joinpath(output_path, "Summary_stats.csv"), output)
end

function write_missing_power(
    problem::PSI.OperationsProblem{T},
    system_ha::PSY.System,
    data_path::String,
    output_path::String,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        StochasticUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
    },
}
    case_initial_time = PSI.get_initial_time(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    time_steps = PSI.model_time_steps(optimization_container)

    initial_time = DateTime("2018-01-01")
    d = Day(case_initial_time - initial_time).value + 1
    total_solar_GW = h5open(joinpath(data_path, "ERCOT132.h5"), "r") do file
        return read(file, "Power")[(288 * (d - 1) + 1):(288 * (d - 1) + length(
            time_steps,
        ) * 12)] ./ 1000
    end

    if length(size(jump_model.obj_dict[:pg])) == 2
        total_thermal =
            sum.(
                eachrow(
                    get_thermal_generator_power_dataframe(problem, time_steps, nothing),
                ),
            )
    else
        nscenarios = 1:size(pg)[2]
        total_thermal = zeros(length(time_steps))
        for scenario in scenarios
            total_thermal .+=
                sum.(
                    eachrow(
                        get_thermal_generator_power_dataframe(
                            problem,
                            time_steps,
                            scenario,
                        ),
                    ),
                )
        end
        total_thermal ./= length(scenarios)
    end

    total_thermal_online = PSI.axis_array_to_dataframe(jump_model.obj_dict[:ug], [:ug])
    for n in names(total_thermal_online)
        total_thermal_online[!, n] .*=
            get_active_power_limits(get_component(ThermalMultiStart, system, n)).max
    end
    total_thermal_online = sum.(eachrow(total_thermal_online))

    # Enforce battery charging as committed
    pb_in =
        sum.(eachrow(PSI.axis_array_to_dataframe(jump_model.obj_dict[:pb_in], [:pb_in])))
    pb_out =
        sum.(eachrow(PSI.axis_array_to_dataframe(jump_model.obj_dict[:pb_out], [:pb_out])))

    _missing_power_by_hourly_setpoint(
        problem,
        total_solar_GW,
        pb_in,
        pb_out,
        total_thermal,
        total_thermal_online,
        data_path,
        output_path,
    )

    _missing_power_by_5_min(
        problem,
        system_ha,
        total_solar_GW,
        pb_in,
        pb_out,
        total_thermal,
        total_thermal_online,
        data_path,
        output_path,
    )
end

_missing_power_by_5_min = function (
    problem::PSI.OperationsProblem{T},
    system_ha::PSY.System,
    total_solar_GW::Vector{Float32},
    pb_in::Vector{Float64},
    pb_out::Vector{Float64},
    total_thermal::Vector{Float64},
    total_thermal_online::Vector{Float64},
    data_path::String,
    output_path::String,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        StochasticUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
    },
}
    case_initial_time = PSI.get_initial_time(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    time_steps_per_hour = 12

    # HACK WITH get max active power TO ADJUST HA TIME SERIES UNTIL JOSE UPDATES SYSTEMS
    total_load = zeros(length(time_steps) * time_steps_per_hour)
    for i in time_steps
        for l in get_components(PowerLoad, system_ha)
            total_load[((i - 1) * time_steps_per_hour + 1):(i * time_steps_per_hour)] .+=
                get_time_series_values(
                    Deterministic,
                    l,
                    "max_active_power";
                    start_time = case_initial_time + Hour(i - 1),
                    len = time_steps_per_hour,
                ) .* get_max_active_power(l) .* problem.ext["load_scale"]
        end
    end
    total_hydro = zeros(length(time_steps) * time_steps_per_hour)
    for i in time_steps
        for l in get_components(HydroGen, system_ha)
            total_hydro[((i - 1) * time_steps_per_hour + 1):(i * time_steps_per_hour)] .+=
                get_time_series_values(
                    Deterministic,
                    l,
                    "max_active_power";
                    start_time = case_initial_time + Hour(i - 1),
                    len = time_steps_per_hour,
                ) .* get_max_active_power(l)
        end
    end
    total_wind = zeros(length(time_steps) * time_steps_per_hour)
    for i in time_steps
        for l in get_components(
            RenewableGen,
            system_ha,
            x -> get_prime_mover(x) != PrimeMovers.PVe,
        )
            total_wind[((i - 1) * time_steps_per_hour + 1):(i * time_steps_per_hour)] .+=
                get_time_series_values(
                    Deterministic,
                    l,
                    "max_active_power";
                    start_time = case_initial_time + Hour(i - 1),
                    len = time_steps_per_hour,
                ) .* get_max_active_power(l)
        end
    end

    # Comparison by set-point
    missing_power =
        (
            total_load - total_hydro - total_wind +
            repeat(pb_in - pb_out - total_thermal, inner = 12)
        ) .* 100 ./ 1000 - total_solar_GW
    missing_power =
        DataFrame(transpose(reshape(missing_power, (12, length(time_steps)))), :auto)
    CSV.write(joinpath(output_path, "missing_power__5min_setpoint.csv"), missing_power)

    missing_capacity =
        (
            total_load - total_hydro - total_wind +
            repeat(pb_in - pb_out - total_thermal_online, inner = 12)
        ) .* 100 ./ 1000 - total_solar_GW
    missing_capacity =
        DataFrame(transpose(reshape(missing_capacity, (12, length(time_steps)))), :auto)
    CSV.write(
        joinpath(output_path, "missing_capacity__5min_setpoint.csv"),
        missing_capacity,
    )
end

_missing_power_by_hourly_setpoint = function (
    problem::PSI.OperationsProblem{T},
    total_solar_GW::Vector{Float32},
    pb_in::Vector{Float64},
    pb_out::Vector{Float64},
    total_thermal::Vector{Float64},
    total_thermal_online::Vector{Float64},
    data_path::String,
    output_path::String,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        StochasticUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
    },
}
    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    time_steps = PSI.model_time_steps(optimization_container)

    total_load =
        get_area_total_time_series(problem, PowerLoad) .* problem.ext["load_scale"]
    total_hydro = get_area_total_time_series(problem, HydroGen)
    pW = PSI.axis_array_to_dataframe(jump_model.obj_dict[:pW], [:pW])

    net_GW =
        (total_load + pb_in - pb_out - total_thermal - total_hydro - pW[:, 1]) .* 100 ./
        1000

    # Comparison by set-point
    missing_power = repeat(net_GW, inner = 12) - total_solar_GW
    missing_power =
        DataFrame(transpose(reshape(missing_power, (12, length(time_steps)))), :auto)
    CSV.write(joinpath(output_path, "missing_power__hourly_setpoint.csv"), missing_power)

    # Comparison by max capacity online
    total_wind = get_area_total_time_series(
        problem,
        RenewableGen;
        filter = x -> get_prime_mover(x) != PrimeMovers.PVe,
    )

    net_GW_online =
        (total_load + pb_in - pb_out - total_thermal_online - total_hydro - total_wind) .* 100 ./ 1000

    missing_capacity = repeat(net_GW_online, inner = 12) - total_solar_GW
    missing_capacity =
        DataFrame(transpose(reshape(missing_capacity, (12, length(time_steps)))), :auto)
    CSV.write(
        joinpath(output_path, "missing_capacity__hourly_setpoint.csv"),
        missing_capacity,
    )
end

function save_as_initial_condition(problem::PSI.OperationsProblem{T}, fname, hour) where {T}
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    ug = PSI.axis_array_to_dataframe(jump_model.obj_dict[:ug], [:ug])[[hour], :]
    CSV.write(fname, ug)
end
