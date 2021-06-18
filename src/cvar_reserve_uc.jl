struct CVaRReserveUnitCommitmentCC <: PSI.PowerSimulationsOperationsProblem end

function PSI.problem_build!(problem::PSI.OperationsProblem{CVaRReserveUnitCommitmentCC};)
    use_storage = problem.ext["use_storage"]
    use_storage_reserves = problem.ext["use_storage_reserves"]
    storage_reserve_names = problem.ext["storage_reserve_names"]
    use_reg = problem.ext["use_reg"]
    use_spin = problem.ext["use_spin"]
    C_RR = problem.ext["C_RR"]
    α = problem.ext["α"]
    C_res_penalty = problem.ext["C_res_penalty"]
    C_ener_penalty = problem.ext["C_ener_penalty"]
    supp_type = problem.ext["supp_type"]

    if use_storage_reserves && !use_storage
        throw(ArgumentError("Can only add storage to reserves if use_storage is true"))
    end
    if use_storage_reserves && !(use_reg || use_spin)
        throw(ArgumentError("Can only use storage reserves if a reserve category is on"))
    end
    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    # Hack to reset the internal JuMP model during development. Do not remove.
    problem.internal.optimization_container.JuMPmodel = JuMP.Model()
    PSI.optimization_container_init!(
        optimization_container,
        PSI.CopperPlatePowerModel,
        system,
    )
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)
    case_initial_time = PSI.get_initial_time(problem)

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
    reg⁺_device_names = get_name.(get_components(ThermalMultiStart, system, x -> !PSY.get_must_run(x)))
    reg⁻_device_names = get_name.(get_components(ThermalMultiStart, system, x -> !PSY.get_must_run(x)))
    spin_device_names = get_name.(get_components(ThermalMultiStart, system, x -> !PSY.get_must_run(x)))

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
    apply_solar!(problem,
        required_reg⁺,
        required_reg⁻,
        required_spin
    )
    pS = jump_model.obj_dict[:pS]
    scenarios = 1:size(pS)[1]

    # -------------------------------------------------------------
    # Variables
    # -------------------------------------------------------------
    ug = JuMP.@variable(
        jump_model,
        ug[g in thermal_gen_names, t in time_steps],
        binary = true
    )
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
    pg = JuMP.@variable(jump_model, pg[g in thermal_gen_names, j in scenarios, t in time_steps] >= 0) # power ABOVE MINIMUM
    pW = JuMP.@variable(jump_model, pW[t in time_steps] >= 0)
    if use_reg
        reg⁺ = JuMP.@variable(
            jump_model,
            reg⁺[
                # g in (use_storage_reserves ? union(reg⁺_device_names, storage_reserve_names) :
                #         reg⁺_device_names),
                g in reg⁺_device_names,
                t in time_steps,
            ] >= 0
        )
        reg⁻ = JuMP.@variable(
            jump_model,
            reg⁻[
                # g in (use_storage_reserves ? union(reg⁻_device_names, storage_reserve_names) :
                #         reg⁻_device_names),
                g in reg⁻_device_names,
                t in time_steps,
            ] >= 0
        )
    end
    if use_spin
        spin = JuMP.@variable(
            jump_model,
            spin[
                g in (use_storage_reserves ? union(spin_device_names, storage_reserve_names) :
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

    if supp_type == "nonspin"
        supp = JuMP.@variable(
            jump_model,
            supp[g in thermal_gen_names, j in scenarios, t in time_steps] >= 0
        )
        # Non-spinning reserve binary variable. Required because of CC exclusion set;
        # otherwise could use 1 - ug
        # Trying for now indexed by j; could change that to just being by g, t
        ng = JuMP.@variable(
            jump_model,
            ng[g in thermal_gen_names, j in scenarios, t in time_steps],
            binary = true
        )
        total_supp = JuMP.@expression(
            jump_model,
            [j in scenarios, t in time_steps],
            sum(supp[g, j, t] for g in thermal_gen_names)
        )
        optimization_container.expressions[:total_supp] = total_supp
    else # generic
        total_supp = JuMP.@variable(
            jump_model,
            total_supp[j in scenarios, t in time_steps] >= 0
        )
    end

    z = JuMP.@variable(jump_model, z[j in scenarios] >= 0) # Eq (25)
    β = JuMP.@variable(jump_model, β)
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

    apply_wind!(problem,
        required_reg⁺,
        required_reg⁻,
        required_spin
    )

    apply_thermal_constraints!(problem,
        spin_device_names
    )
    Cg = optimization_container.expressions[:Cg]

    if use_storage
        apply_storage!(problem, storage_reserve_names)
        pb_in = jump_model.obj_dict[:pb_in]
        pb_out = jump_model.obj_dict[:pb_out]
    end

    if use_reg
        apply_reg_requirements!(problem,
            reg⁺_device_names,
            reg⁻_device_names,
            required_reg⁺,
            required_reg⁻,
            storage_reserve_names
        )
    end

    if use_spin
        apply_spin_requirements!(problem,
            spin_device_names,
            required_spin,
            storage_reserve_names
        )
    end

    # Power balance constraint, hydro included
    power_balance_constraint = JuMP.@constraint(
        jump_model,
        [j in scenarios, t in time_steps],
        sum(pg[g, j, t] + pg_lim[g].min * ug[g, t] for g in thermal_gen_names) +
        pS[j, t] +
        pW[t] + 
        total_hydro[t] +
        (use_slack ? slack_energy⁺[t] : 0) +
        (use_storage ? sum(pb_out[b, t] - pb_in[b, t] for b in storage_names) : 0) ==
        total_load[t] + (use_slack ? slack_energy⁻[t] : 0)
    )

    # Auxiliary variable definition
    auxiliary_constraint = JuMP.@constraint(
        jump_model,
        [j in scenarios],
        z[j] >= sum(total_supp[j, t] for t in time_steps) - β
    )

    obj_function = JuMP.@objective(
        jump_model,
        Min,
        sum(
            no_load_cost[g] * ug[g, t] +
            (1 / length(scenarios)) * sum(Cg[g, j, t] for j in scenarios) +
            sum(startup_cost[g][s] * δ_sg[g, s, t] for s in startup_categories) +
            shutdown_cost[g] * wg[g, t] for g in thermal_gen_names, t in time_steps
        ) +
        C_RR * (β + 1 / (length(scenarios) * (1 - α)) * sum(z[j] for j in scenarios)) +
        (
            use_slack ?
            C_res_penalty *
            sum(slack_reg⁺[t] + slack_reg⁻[t] + slack_spin[t] for t in time_steps) +
            C_ener_penalty * sum(slack_energy⁺[t] + slack_energy⁻[t] for t in time_steps) : 0
        )
    )
end
