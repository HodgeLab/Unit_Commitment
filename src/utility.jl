
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

function save_as_initial_condition(problem::PSI.OperationsProblem{T}, fname, hour) where {T}
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    ug = PSI.axis_array_to_dataframe(jump_model.obj_dict[:ug], [:ug])[[hour], :]
    CSV.write(fname, ug)
end

function write_reserve_summary(
    problem::PSI.OperationsProblem{T},
    output_path::String;
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
        StochasticUnitCommitmentCC,
    },
}
    use_storage_reserves = problem.ext["use_storage_reserves"]
    use_solar_reg = problem.ext["use_solar_reg"]
    use_solar_spin = problem.ext["use_solar_spin"]

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)
    obj_dict = PSI.get_jump_model(optimization_container).obj_dict
    case_initial_time = PSI.get_initial_time(problem)
    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    storage_reserve_names = problem.ext["storage_reserve_names"]
    time_steps = 1:24

    reg_reserve_up = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "REG_UP")
    reg_reserve_dn =
        PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, "REG_DN")
    spin_reserve = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "SPIN")
    total_required_reg⁺ = sum(get_time_series_values(
        Deterministic,
        reg_reserve_up,
        "requirement";
        start_time = case_initial_time,
        len = length(time_steps)
    ))
    total_required_reg⁻ = sum(get_time_series_values(
        Deterministic,
        reg_reserve_dn,
        "requirement";
        start_time = case_initial_time,
        len = length(time_steps)
    ))
    total_required_spin = sum(get_time_series_values(
        Deterministic,
        spin_reserve,
        "requirement";
        start_time = case_initial_time,
        len = length(time_steps)
    ))

    output = Dict{String, Float64}()
    output["Required REG_UP [GWh]"] = total_required_reg⁺
    output["Required REG_DN [GWh]"] = total_required_reg⁻
    output["Required SPIN [GWh]"] = total_required_spin

    reg⁺ = PSI.axis_array_to_dataframe(obj_dict[:reg⁺], [:reg⁺])[time_steps, :]
    reg⁻ = PSI.axis_array_to_dataframe(obj_dict[:reg⁻], [:reg⁻])[time_steps, :]
    spin = PSI.axis_array_to_dataframe(obj_dict[:spin], [:spin])[time_steps, :]

    output["Thermal REG_UP [GWh]"] = sum(sum.(eachrow(
        (reg⁺[:, intersect(names(reg⁺), thermal_gen_names)]))))
    output["Thermal REG_DN [GWh]"] = sum(sum.(eachrow(
        (reg⁻[:, intersect(names(reg⁻), thermal_gen_names)]))))
    output["Thermal SPIN [GWh]"] = sum(sum.(eachrow(
        (spin[:, intersect(names(spin), thermal_gen_names)]))))

    if use_slack
        output["Unserved REG_UP [GWh]"] = sum.eachcol(
            PSI.axis_array_to_dataframe(obj_dict[:slack_reg⁺], [:slack_reg⁺])[time_steps, :])
        output["Unserved REG_DN [GWh]"] = sum.eachcol(
            PSI.axis_array_to_dataframe(obj_dict[:slack_reg⁻], [:slack_reg⁻])[time_steps, :])
        output["Unserved SPIN [GWh]"] = sum.eachcol(
            PSI.axis_array_to_dataframe(obj_dict[:slack_spin], [:slack_spin])[time_steps, :])
    end
    if use_solar_spin
        output["Solar SPIN [GWh]"] = sum.eachcol(
            PSI.axis_array_to_dataframe(obj_dict[:spin_S], [:spin_S])[time_steps, :])
    end
    if use_solar_reg
        output["Solar REG_UP [GWh]"] = sum.eachcol(
            PSI.axis_array_to_dataframe(obj_dict[:reg⁺_S], [:reg⁺_S])[time_steps, :])
        output["Solar REG_DN [GWh]"] = sum.eachcol(
            PSI.axis_array_to_dataframe(obj_dict[:reg⁻_S], [:reg⁻_S])[time_steps, :])
    end
    if use_storage_reserves
        output["Storage SPIN [GWh]"] = sum(sum.(eachrow(
        (spin[:, intersect(names(spin), storage_reserve_names)]))))
    end

    # Scale from 100 MW to GW
    for v in keys(output)
        output[v] *= get_base_power(system) / 1000
    end

    output = _calc_reserve_percentages!(output,
        use_slack,
        use_solar_spin,
        use_solar_reg,
        use_storage_reserves
        )

    CSV.write(
        joinpath(output_path, "reserve_summary.csv"),
        output,
    )
end

function write_reserve_summary(
    res::PSI.SimulationProblemResults,
    system::PSY.System,
    output_path::String,
    use_slack::Bool,
    use_solar_reg::Bool,
    use_solar_spin::Bool,
    use_storage_reserves::Bool;
)
    Δt = 1 / 12

    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    storage_names = PSY.get_name.(get_components(PSY.GenericBattery, system))
    solar_names = PSY.get_name.(get_components(RenewableGen, system, x -> get_prime_mover(x) == PrimeMovers.PVe))

    output = Dict{String, Float64}()

    output["Required REG_UP [GWh]"] = _get_total_hauc_reserve_requirement(
        res,
        PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "REG_UP")
        )
    output["Required REG_DN [GWh]"] = _get_total_hauc_reserve_requirement(
        res,
        PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, "REG_DN")
        )
    output["Required SPIN [GWh]"] = _get_total_hauc_reserve_requirement(
        res,
        PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "SPIN")
        )

    output["Thermal REG_UP [GWh]"] = _2D_total(
        res,
        :REG_UP__VariableReserve_ReserveUp;
        filter = thermal_gen_names
    )
    output["Thermal REG_DN [GWh]"] = _2D_total(
        res,
        :REG_DN__VariableReserve_ReserveDown;
        filter = thermal_gen_names
    )
    output["Thermal SPIN [GWh]"] = _2D_total(
        res,
        :SPIN__VariableReserve_ReserveUp;
        filter = thermal_gen_names
    )

    if use_slack
        output["Unserved REG_UP [GWh]"] = _1D_total(res, :γ⁺__REG_UP)
        output["Unserved REG_DN [GWh]"] = _1D_total(res, :γ⁺__REG_DN)
        output["Unserved SPIN [GWh]"] = _1D_total(res, :γ⁺__SPIN)
    end
    if use_solar_spin
        output["Solar SPIN [GWh]"] = _2D_total(
            res,
            :SPIN__VariableReserve_ReserveUp;
            filter = solar_names
        )
    end
    if use_solar_reg
        output["Solar REG_UP [GWh]"] = _2D_total(
            res,
            :REG_UP__VariableReserve_ReserveUp;
            filter = solar_names
        )
        output["Solar REG_DN [GWh]"] = _2D_total(
            res,
            :REG_DN__VariableReserve_ReserveDown;
            filter = solar_names
        )
    end
    if use_storage_reserves
        output["Storage SPIN [GWh]"] = _2D_total(
            res,
            :SPIN__VariableReserve_ReserveUp;
            filter = storage_names
        )
    end

    for v in keys(output)
        output[v] *= Δt * get_base_power(system) / 1000
    end

    output = _calc_reserve_percentages!(output,
        use_slack,
        use_solar_spin,
        use_solar_reg,
        use_storage_reserves
        )

    CSV.write(
        joinpath(output_path, "reserve_summary.csv"),
        output,
    )
end

function _2D_total(res::PSI.SimulationProblemResults, sym::Symbol; filter = nothing)
    temp = read_realized_variables(res, names = [sym])[sym]
    if isnothing(filter)
        temp = temp[:, setdiff(names(temp), ["DateTime"])]
    else
        temp = temp[:, intersect(names(temp), filter)]
    end
    return sum(sum.(eachrow(temp)))
end

function _1D_total(res::PSI.SimulationProblemResults, sym::Symbol)
    temp = read_realized_variables(res, names = [sym])[sym]
    return sum(sum.(eachcol(temp[:, setdiff(names(temp), ["DateTime"])])))
end

function _get_total_hauc_reserve_requirement(
    res::PSI.SimulationProblemResults,
    reserve_component
    )

    timestamps = get_realized_timestamps(res)
    hour_timestamps = collect(first(timestamps):Hour(1):last(timestamps))

    required_reserve = zeros(length(timestamps))
    for i in 1:length(hour_timestamps)
        required_reserve[((i - 1) * 12 + 1):(i * 12)] =
            get_time_series_values(
                Deterministic,
                reserve_component,
                "requirement";
                start_time = hour_timestamps[i],
                len = 12,
            )
    end

    return sum(required_reserve)
end

function _calc_reserve_percentages!(output::Dict,
    use_slack::Bool,
    use_solar_spin::Bool,
    use_solar_reg::Bool,
    use_storage_reserves::Bool
    )

    total_spin = output["Thermal SPIN [GWh]"] + (use_slack ? output["Unserved SPIN [GWh]"] : 0.0) +
        (use_solar_spin ? output["Solar SPIN [GWh]"] : 0.0) +
        (use_storage_reserves ? output["Storage SPIN [GWh]"] : 0.0)
    total_reg_up = output["Thermal REG_UP [GWh]"] + (use_slack ? output["Unserved REG_UP [GWh]"] : 0.0) +
        (use_solar_reg ? output["Solar REG_UP [GWh]"] : 0.0)
    total_reg_dn = output["Thermal REG_DN [GWh]"] + (use_slack ? output["Unserved REG_DN [GWh]"] : 0.0) +
        (use_solar_reg ? output["Solar REG_DN [GWh]"] : 0.0)

    if use_slack
        output["Unserved SPIN [%]"] = output["Unserved SPIN [GWh]"] / total_spin * 100
        output["Unserved REG_UP [%]"] = output["Unserved REG_UP [GWh]"] / total_reg_up * 100
        output["Unserved REG_DN [%]"] = output["Unserved REG_DN [GWh]"] / total_reg_dn * 100
    end

    output["Thermal SPIN [%]"] = output["Thermal SPIN [GWh]"] / total_spin * 100
    output["Thermal REG_UP [%]"] = output["Thermal REG_UP [GWh]"] / total_reg_up * 100
    output["Thermal REG_DN [%]"] = output["Thermal REG_DN [GWh]"] / total_reg_dn * 100

    if use_solar_spin
        output["Solar SPIN [%]"] = output["Solar SPIN [GWh]"] / total_spin * 100
    end
    if use_solar_reg
        output["Solar REG_UP [%]"] = output["Solar REG_UP [GWh]"] / total_reg_up * 100
        output["Solar REG_DN [%]"] = output["Solar REG_DN [GWh]"] / total_reg_dn * 100
    end
    if use_storage_reserves
        output["Storage SPIN [%]"] = output["Storage SPIN [GWh]"] / total_spin * 100
    end

    return output
end