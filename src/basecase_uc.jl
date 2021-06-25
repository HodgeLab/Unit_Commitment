struct BasecaseUnitCommitmentCC <: PSI.PowerSimulationsOperationsProblem end
struct HourAheadUnitCommitmentCC <: PSI.PowerSimulationsOperationsProblem end

function PSI.problem_build!(problem::PSI.OperationsProblem{BasecaseUnitCommitmentCC};)
    _build_basecase_internal!(problem)
end

function PSI.problem_build!(problem::PSI.OperationsProblem{HourAheadUnitCommitmentCC};)
    _build_basecase_internal!(problem)

    _enforce_ha_commitments!(problem)
end

function _build_basecase_internal!(problem::PSI.OperationsProblem{T}
    ) where T <: Union{BasecaseUnitCommitmentCC, HourAheadUnitCommitmentCC}
    use_storage = problem.ext["use_storage"]
    use_storage_reserves = problem.ext["use_storage_reserves"]
    storage_reserve_names = problem.ext["storage_reserve_names"]
    use_reg = problem.ext["use_reg"]
    use_spin = problem.ext["use_spin"]
    C_res_penalty = problem.ext["C_res_penalty"]
    C_ener_penalty = problem.ext["C_ener_penalty"]

    if use_storage_reserves && !use_storage
        throw(ArgumentError("Can only add storage to reserves if use_storage is true"))
    end
    if use_storage_reserves && !(use_reg || use_spin)
        throw(ArgumentError("Can only use storage reserves if a reserve category is on"))
    end
    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    # Hack to reset the internal JuMP model during development. Do not remove.
    PSI.optimization_container_init!(
        optimization_container,
        PSI.CopperPlatePowerModel,
        system,
    )
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)
    case_initial_time = PSI.get_initial_time(problem)

    ug_t0 = problem.ext["init_conditions"][:ug_t0]
    Pg_t0 = problem.ext["init_conditions"][:Pg_t0]
    time_up_t0 = problem.ext["init_conditions"][:time_up_t0]
    time_down_t0 = problem.ext["init_conditions"][:time_down_t0]

    # -------------------------------------------------------------
    # Collect definitions from PSY model
    # -------------------------------------------------------------
    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    storage_names = PSY.get_name.(get_components(PSY.GenericBattery, system))
    pg_lim = Dict(
        g => get_active_power_limits(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    startup_categories = (:hot, :warm, :cold)
    no_load_cost = Dict(
        g => get_no_load(get_operation_cost(get_component(ThermalMultiStart, system, g))) for g in thermal_gen_names
    )
    shutdown_cost = Dict(
        g => get_shut_down(get_operation_cost(get_component(ThermalMultiStart, system, g))) for g in thermal_gen_names
    )
    startup_cost = Dict(
        g => get_start_up(get_operation_cost(get_component(ThermalMultiStart, system, g))) for g in thermal_gen_names
    )

    reg_reserve_up = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "REG_UP")
    reg_reserve_dn =
        PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, "REG_DN")
    spin_reserve = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "SPIN")
    # reg⁺_device_names = get_name.(get_contributing_devices(system, reg_reserve_up))
    # reg⁻_device_names = get_name.(get_contributing_devices(system, reg_reserve_dn))
    # spin_device_names = get_name.(get_contributing_devices(system, spin_reserve))
    reg⁺_device_names =
        get_name.(get_components(ThermalMultiStart, system, x -> !PSY.get_must_run(x)))
    reg⁻_device_names =
        get_name.(get_components(ThermalMultiStart, system, x -> !PSY.get_must_run(x)))
    spin_device_names =
        get_name.(get_components(ThermalMultiStart, system, x -> !PSY.get_must_run(x)))

    # -------------------------------------------------------------
    # Time-series data
    # -------------------------------------------------------------

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
    required_spin = get_time_series_values(
        Deterministic,
        spin_reserve,
        "requirement";
        start_time = case_initial_time,
    )
    total_load = get_area_total_time_series(problem, PowerLoad) .* problem.ext["load_scale"]
    total_hydro = get_area_total_time_series(problem, HydroGen)

    # Begin with solar equations
    apply_solar!(problem, required_reg⁺, required_reg⁻, required_spin)
    pS = jump_model.obj_dict[:pS]

    # -------------------------------------------------------------
    # Variables
    # -------------------------------------------------------------
    ug = JuMP.@variable(
        jump_model,
        ug[g in thermal_gen_names, t in time_steps],
        binary = true
    )

    optimization_container.variables[:On__ThermalMultiStart] = ug
    vg = JuMP.@variable(
        jump_model,
        vg[g in thermal_gen_names, t in time_steps],
        binary = true
    )
    wg = JuMP.@variable(
        jump_model,
        wg[g in thermal_gen_names, t in time_steps],
        binary = true
    )
    δ_sg = JuMP.@variable(
        jump_model,
        δ_sg[g in thermal_gen_names, s in startup_categories, t in time_steps],
        binary = true
    )
    pg = JuMP.@variable(jump_model, pg[g in thermal_gen_names, t in time_steps] >= 0) # power ABOVE MINIMUM
    pW = JuMP.@variable(jump_model, pW[t in time_steps] >= 0)

    if use_reg
        reg⁺ = JuMP.@variable(jump_model, reg⁺[
            # g in (use_storage_reserves ? union(reg⁺_device_names, storage_reserve_names) :
            #         reg⁺_device_names),
            g in reg⁺_device_names,
            t in time_steps,
        ] >= 0)
        reg⁻ = JuMP.@variable(jump_model, reg⁻[
            # g in (use_storage_reserves ? union(reg⁻_device_names, storage_reserve_names) :
            #         reg⁻_device_names),
            g in reg⁻_device_names,
            t in time_steps,
        ] >= 0)
    end
    if use_spin
        spin = JuMP.@variable(
            jump_model,
            spin[
                g in (use_storage_reserves ?
                      union(spin_device_names, storage_reserve_names) :
                      spin_device_names),
                t in time_steps,
            ] >= 0
        )
    end
    # Register expression of total reserves for each generator
    total_reserve⁺ = JuMP.@expression(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        (g in reg⁺_device_names ? reg⁺[g, t] : 0) +
        (g in spin_device_names ? spin[g, t] : 0)
    )
    total_reserve⁻ = JuMP.@expression(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        (g in reg⁻_device_names ? reg⁻[g, t] : 0)
    )
    optimization_container.expressions[:total_reserve⁺] = total_reserve⁺
    optimization_container.expressions[:total_reserve⁻] = total_reserve⁻

    if use_slack
        slack_reg⁺ = JuMP.@variable(jump_model, slack_reg⁺[t in time_steps] >= 0)
        slack_reg⁻ = JuMP.@variable(jump_model, slack_reg⁻[t in time_steps] >= 0)
        slack_spin = JuMP.@variable(jump_model, slack_spin[t in time_steps] >= 0)
        slack_energy⁺ = JuMP.@variable(jump_model, slack_energy⁺[t in time_steps] >= 0)
        slack_energy⁻ = JuMP.@variable(jump_model, slack_energy⁻[t in time_steps] >= 0)
    end

    # -------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------

    apply_wind!(problem, required_reg⁺, required_reg⁻, required_spin)

    apply_thermal_constraints!(problem,
        spin_device_names,
        ug_t0,
        Pg_t0,
        time_up_t0,
        time_down_t0)
    Cg = optimization_container.expressions[:Cg]

    if use_storage
        apply_storage!(problem, storage_reserve_names)
        pb_in = jump_model.obj_dict[:pb_in]
        pb_out = jump_model.obj_dict[:pb_out]
    end

    if use_reg
        apply_reg_requirements!(
            problem,
            reg⁺_device_names,
            reg⁻_device_names,
            required_reg⁺,
            required_reg⁻,
            storage_reserve_names,
        )
    end

    if use_spin
        apply_spin_requirements!(
            problem,
            spin_device_names,
            required_spin,
            storage_reserve_names,
        )
    end

    # Power balance constraint, hydro included
    power_balance_constraint = JuMP.@constraint(
        jump_model,
        [t in time_steps],
        sum(pg[g, t] + pg_lim[g].min * ug[g, t] for g in thermal_gen_names) +
        pS[t] +
        pW[t] +
        total_hydro[t] +
        (use_slack ? slack_energy⁺[t] : 0) +
        (use_storage ? sum(pb_out[b, t] - pb_in[b, t] for b in storage_names) : 0) ==
        total_load[t] + (use_slack ? slack_energy⁻[t] : 0)
    )

    obj_function = JuMP.@objective(
        jump_model,
        Min,
        sum(
            no_load_cost[g] * ug[g, t] +
            Cg[g, t] +
            sum(startup_cost[g][s] * δ_sg[g, s, t] for s in startup_categories) +
            shutdown_cost[g] * wg[g, t] for g in thermal_gen_names, t in time_steps
        ) + (
            use_slack ?
            C_res_penalty *
            sum(slack_reg⁺[t] + slack_reg⁻[t] + slack_spin[t] for t in time_steps) +
            C_ener_penalty * sum(slack_energy⁺[t] + slack_energy⁻[t] for t in time_steps) : 0
        )
    )
end

function _enforce_ha_commitments!(problem::PSI.OperationsProblem{HourAheadUnitCommitmentCC};)
    h = problem.ext["step"]
    DA_obj_dict = problem.ext["UC_obj_dict"]

    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    obj_dict = jump_model.obj_dict

    # Collect HA variables
    ug = obj_dict[:ug]
    reg⁺ = obj_dict[:reg⁺]
    reg⁻ = obj_dict[:reg⁻]
    spin = obj_dict[:spin]
    eb = obj_dict[:eb]

    # In each of these, the dataframe reverse the dimension order.

    # 1) on/off
    ug_DA = PSI.axis_array_to_dataframe(DA_obj_dict[:ug], [:ug])[h:(h + 1), :]
    for g in names(ug_DA), t_da in 1:2
        if ug_DA[t_da, g] == 1
            for t_ha in 1:12
                if JuMP.solver_name(jump_model) == "Xpress"
                    JuMP.@constraint(
                        jump_model,
                        ug[g, t_ha + (t_da - 1) * 12] >= 1
                    )
                else
                    JuMP.fix(ug[g, t_ha + (t_da - 1) * 12], 1.0; force = true)
                end
            end
        end
    end

    # 2) reg up/down and spin for thermal and battery
    reg⁺_DA = PSI.axis_array_to_dataframe(DA_obj_dict[:reg⁺], [:reg⁺])[h:(h + 1), :]
    for g in names(reg⁺_DA), t_da in 1:2
        for t_ha in 1:12
            JuMP.@constraint(
                jump_model,
                reg⁺[g, t_ha + (t_da - 1) * 12] >= reg⁺_DA[t_da, g]
            )
        end
    end
    reg⁻_DA = PSI.axis_array_to_dataframe(DA_obj_dict[:reg⁻], [:reg⁻])[h:(h + 1), :]
    for g in names(reg⁻_DA), t_da in 1:2
        for t_ha in 1:12
            JuMP.@constraint(
                jump_model,
                reg⁻[g, t_ha + (t_da - 1) * 12] >= reg⁻_DA[t_da, g]
            )
        end
    end
    spin_DA = PSI.axis_array_to_dataframe(DA_obj_dict[:spin], [:spin])[h:(h + 1), :]
    for g in names(spin_DA), t_da in 1:2
        for t_ha in 1:12
            JuMP.@constraint(
                jump_model,
                spin[g, t_ha + (t_da - 1) * 12] >= spin_DA[t_da, g]
            )
        end
    end

    # 3) eb at end of hour
    eb_DA = PSI.axis_array_to_dataframe(DA_obj_dict[:eb], [:eb])[h:(h + 1), :]
    for b in names(eb_DA), t_da in 1:2
        JuMP.@constraint(
            jump_model,
            eb[b, t_da * 12] == eb_DA[t_da, b]
        )
    end

end

function get_first_init_conditions(system)
    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))

    # initial conditions
    ug_t0 = Dict(
        g => PSY.get_status(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    time_up_t0 = Dict(
        g => ug_t0[g] * get_time_at_status(get_component(ThermalMultiStart, system, g))
        for g in thermal_gen_names
    )
    time_down_t0 = Dict(
        g =>
            (1 - ug_t0[g]) *
            get_time_at_status(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    # This is just power (Pg), not power above minimum (pg)
    Pg_t0 = Dict(
        g => get_active_power(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )

    return Dict{Symbol, Dict}(
        :ug_t0 => ug_t0,
        :Pg_t0 => Pg_t0,
        :time_up_t0 => time_up_t0,
        :time_down_t0 => time_down_t0
        )
end

function get_hauc_ending_conditions(problem::PSI.OperationsProblem{HourAheadUnitCommitmentCC})
    system = PSI.get_system(problem)
    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    optimization_container = PSI.get_optimization_container(problem)
    obj_dict = PSI.get_jump_model(optimization_container).obj_dict

    ug_df = PSI.axis_array_to_dataframe(obj_dict[:ug], [:ug])[[12], :]
    ug_t0 = Dict(
            g => convert(Bool, ug_df[1, g]) for
            g in thermal_gen_names
        )

    Pg_t0_df = get_thermal_generator_power_dataframe(problem, 12:12, nothing)
    Pg_t0 = Dict(
        g => Pg_t0_df[1, g] for
        g in thermal_gen_names
    )

    ug_df = PSI.axis_array_to_dataframe(obj_dict[:ug], [:ug])[1:12, :]
    time_up_t0 = Dict(
        g => get_new_time_up(problem, ug_df, g)
        for g in thermal_gen_names
    )
    time_down_t0 = Dict(
        g => get_new_time_down(problem, ug_df, g)
        for g in thermal_gen_names
    )

    return Dict{Symbol, Dict}(
        :ug_t0 => ug_t0,
        :Pg_t0 => Pg_t0,
        :time_up_t0 => time_up_t0,
        :time_down_t0 => time_down_t0
        )
end

function get_new_time_down(problem, ug_df, g)
    MINS_IN_HOUR = 60.0
    Δt = Minute(PSI.model_resolution(PSI.get_optimization_container(problem))).value / MINS_IN_HOUR

    if round(ug_df[12, g]) == 0.0
        # If it has down the whole time, add 1 hour to down time
        if all(round.(ug_df[:, g]) .== 0.0)
            return problem.ext["init_conditions"][:time_down_t0][g] + 1.0
        else
            return first(findall(reverse(round.(ug_df[:, g])) .== 1.0)) * Δt
        end
    else
        return 0.0
    end
end

function get_new_time_up(problem, ug_df, g)
    MINS_IN_HOUR = 60.0
    Δt = Minute(PSI.model_resolution(PSI.get_optimization_container(problem))).value / MINS_IN_HOUR

    if round(ug_df[12, g]) == 1.0
        # If it has up the whole time, add 1 hour to up time
        if all(round.(ug_df[:, g]) .== 1.0)
            return problem.ext["init_conditions"][:time_up_t0][g] + 1.0
        else
            return first(findall(reverse(round.(ug_df[:, g])) .== 0.0)) * Δt
        end
    else
        return 0.0
    end
end
