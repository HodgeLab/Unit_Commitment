struct StochasticUnitCommitmentCC <: PSI.PowerSimulationsOperationsProblem end

function PSI.problem_build!(problem::PSI.OperationsProblem{StochasticUnitCommitmentCC};)
    use_storage = problem.ext["use_storage"]
    use_storage_reserves = problem.ext["use_storage_reserves"]
    storage_reserve_names = problem.ext["storage_reserve_names"]
    use_reg = problem.ext["use_reg"]
    use_spin = problem.ext["use_spin"]
    use_must_run = problem.ext["use_must_run"]
    C_res_penalty = problem.ext["C_res_penalty"]
    C_ener_penalty = problem.ext["C_ener_penalty"]
    L_REG = problem.ext["L_REG"]
    L_SPIN = problem.ext["L_SPIN"]

    if use_storage_reserves && !use_storage
        throw(ArgumentError("Can only add storage to reserves if use_storage is true"))
    end
    if use_storage_reserves && !(use_reg || use_spin)
        throw(ArgumentError("Can only use storage reserves if a reserve category is on"))
    end
    system = PSI.get_system(problem)
    case_initial_time = PSI.get_initial_time(problem)
    optimization_container = PSI.get_optimization_container(problem)
    # Hack to reset the internal JuMP model during development. Do not remove.
    problem.internal.optimization_container.JuMPmodel = nothing
    PSI.optimization_container_init!(
        optimization_container,
        PSI.CopperPlatePowerModel,
        system,
    )
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    resolution = PSI.model_resolution(optimization_container)
    use_slack = PSI.get_balance_slack_variables(optimization_container.settings)

    # Constants
    MINS_IN_HOUR = 60.0
    Δt = 1

    # -------------------------------------------------------------
    # Collect definitions from PSY model
    # -------------------------------------------------------------
    thermal_gens = get_components(ThermalMultiStart, system)
    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    hydro_gens_names = get_name.(get_components(PSY.HydroGen, system))
    storage_names = PSY.get_name.(get_components(PSY.GenericBattery, system))
    get_rmp_up_limit(g) = PSY.get_ramp_limits(g).up
    get_rmp_dn_limit(g) = PSY.get_ramp_limits(g).down
    ramp_up = Dict(
        g =>
            get_rmp_up_limit(get_component(ThermalMultiStart, system, g)) * MINS_IN_HOUR for g in thermal_gen_names
    )
    ramp_dn = Dict(
        g =>
            get_rmp_dn_limit(get_component(ThermalMultiStart, system, g)) * MINS_IN_HOUR for g in thermal_gen_names
    )
    pg_lim = Dict(
        g => get_active_power_limits(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    pg_power_trajectory = Dict(
        g => get_power_trajectory(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    if use_must_run
        must_run_gen_names =
            get_name.(get_components(ThermalMultiStart, system, x -> PSY.get_must_run(x)))
    end
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
    variable_cost = Dict(
        g => get_variable(get_operation_cost(get_component(ThermalMultiStart, system, g))) for g in thermal_gen_names
    )
    # initial conditions
    ug_t0 = Dict(
        g => PSY.get_status(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    # This is just power (Pg), not power above minimum (pg)
    Pg_t0 = Dict(
        g => get_active_power(get_component(ThermalMultiStart, system, g)) for
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

    reg_reserve_up = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "REG_UP")
    reg_reserve_dn =
        PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, "REG_DN")
    spin_reserve = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "SPIN")
    reg⁺_device_names = get_name.(get_contributing_devices(system, reg_reserve_up))
    reg⁻_device_names = get_name.(get_contributing_devices(system, reg_reserve_dn))
    spin_device_names = get_name.(get_contributing_devices(system, spin_reserve))
    no_res⁺_device_names = setdiff(
        thermal_gen_names,
        union(
            use_reg ? reg⁺_device_names : Vector{String}(),
            use_spin ? spin_device_names : Vector{String}(),
        ),
    )
    no_res⁻_device_names = setdiff(
        thermal_gen_names,
        (use_reg ? reg⁻_device_names : Vector{String}())
    )

    required_spin = get_time_series_values(
        Deterministic,
        spin_reserve,
        "requirement";
        start_time = case_initial_time,
    )
    # -------------------------------------------------------------
    # Time-series data
    # -------------------------------------------------------------
    total_load = get_area_total_time_series(problem, PowerLoad).*1.15
    total_hydro = get_area_total_time_series(problem, HydroGen)
    total_wind = get_area_total_time_series(
        problem,
        RenewableGen;
        filter = x -> get_prime_mover(x) != PrimeMovers.PVe,
    )

    # Begin with solar equations
    apply_solar!(problem)
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
                g in (use_storage_reserves ? union(reg⁺_device_names, storage_reserve_names) :
                      reg⁺_device_names),
                t in time_steps,
            ] >= 0
        )
        reg⁻ = JuMP.@variable(
            jump_model,
            reg⁻[
                g in (use_storage_reserves ? union(reg⁻_device_names, storage_reserve_names) :
                      reg⁻_device_names),
                t in time_steps,
            ] >= 0
        )
    end
    if use_spin
        spin = JuMP.@variable(
            jump_model,
            spin[
                g in spin_device_names,
                t in time_steps,
            ] >= 0
        )
    end
    λ = JuMP.@variable(
        jump_model,
        0 <=
        λ[g in thermal_gen_names, j in scenarios, i in 1:length(variable_cost[g]), t in time_steps] <=
        PSY.get_breakpoint_upperbounds(variable_cost[g])[i]
    )
    if use_slack
        slack_reg⁺ = JuMP.@variable(jump_model, slack_reg⁺[t in time_steps] >= 0)
        slack_reg⁻ = JuMP.@variable(jump_model, slack_reg⁻[t in time_steps] >= 0)
        slack_spin = JuMP.@variable(jump_model, slack_spin[t in time_steps] >= 0)
        slack_energy⁺ = JuMP.@variable(jump_model, slack_energy⁺[t in time_steps] >= 0)
        slack_energy⁻ = JuMP.@variable(jump_model, slack_energy⁻[t in time_steps] >= 0)
    end

    if use_storage
        apply_storage!(problem, storage_reserve_names)
        pb_in = jump_model.obj_dict[:pb_in]
        pb_out = jump_model.obj_dict[:pb_out]
    end

    # -------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------

    # Eq (5) Wind constraint
    wind_constraint =
        JuMP.@constraint(jump_model, [t in time_steps], pW[t] <= total_wind[t])

    # Eq (6) PWL variable cost constraint
    # PWL Cost function auxiliary variables
    lambda_constraint = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, j in scenarios, t in time_steps],
        sum(λ[g, j, i, t] for i in 1:length(variable_cost[g])) ==
        pg_lim[g].min * ug[g, t] + pg[g, j, t]
    )

    Cg = JuMP.@expression(
        jump_model,
        [g in thermal_gen_names, j in scenarios, t in time_steps],
        sum(
            100.0 * λ[g, j, i, t] * PSY.get_slopes(variable_cost[g])[i] for
            i in 1:length(variable_cost[g])
        )
    )
    optimization_container.expressions[:Cg] = Cg

    obj_function = JuMP.@objective(
        jump_model,
        Min,
        sum(
            no_load_cost[g] * ug[g, t] +
            (1 / length(scenarios)) * sum(Cg[g, j, t] for j in scenarios) +
            sum(startup_cost[g][s] * δ_sg[g, s, t] for s in startup_categories) +
            shutdown_cost[g] * wg[g, t] for g in thermal_gen_names, t in time_steps
        ) +
        (
            use_slack ?
            C_res_penalty *
            sum(slack_reg⁺[t] + slack_reg⁻[t] + slack_spin[t] for t in time_steps) +
            C_ener_penalty * sum(slack_energy⁺[t] + slack_energy⁻[t] for t in time_steps) : 0
        )
    )

    # Eq (7) Commitment constraints
    commitment_constraints = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        thermal_gen_names,
        time_steps,
    )
    for g in thermal_gen_names, t in time_steps
        if t == 1  # pg_lib (6)
            commitment_constraints[g, 1] =
                JuMP.@constraint(jump_model, ug[g, 1] - ug_t0[g] == vg[g, 1] - wg[g, 1])
        else
            commitment_constraints[g, t] =
                JuMP.@constraint(jump_model, ug[g, t] - ug[g, t - 1] == vg[g, t] - wg[g, t])
        end
    end

    # Eq (8) Must-run -- added as a constraint for Xpress due to fix not working
    if use_must_run
        if JuMP.solver_name(jump_model) == "Xpress"
            JuMP.@constraint(
                jump_model,
                [g in must_run_gen_names, t in time_steps],
                ug[g, t] >= 1
            )
        else
            for g in must_run_gen_names, t in time_steps
                JuMP.fix(ug[g, t], 1.0; force = true)
            end
        end
    end

    # Eq (9) - (10), up and down time constraints
    # Translated from device_duration_compact_retrospective!
    uptime_constraints = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        thermal_gen_names,
        time_steps,
    )
    downtime_constraints = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        thermal_gen_names,
        time_steps,
    )
    for g in thermal_gen_names, t in time_steps
        time_limits = get_time_limits(get_component(ThermalMultiStart, system, g))

        lhs_on = JuMP.AffExpr(0)
        if t in UnitRange{Int}(
            Int(min(time_limits[:up], length(time_steps))),
            length(time_steps),
        )
            for i in UnitRange{Int}(Int(t - time_limits[:up] + 1), t)
                if i in time_steps
                    JuMP.add_to_expression!(lhs_on, vg[g, i])
                end
            end
        elseif t <= max(0, time_limits[:up] - time_up_t0[g]) && time_up_t0[g] > 0
            JuMP.add_to_expression!(lhs_on, 1)
        else
            continue
        end
        uptime_constraints[g, t] = JuMP.@constraint(jump_model, lhs_on - ug[g, t] <= 0.0)

        lhs_off = JuMP.AffExpr(0)
        if t in UnitRange{Int}(
            Int(min(time_limits[:down], length(time_steps))),
            length(time_steps),
        )
            for i in UnitRange{Int}(Int(t - time_limits[:down] + 1), t)
                if i in time_steps
                    JuMP.add_to_expression!(lhs_off, wg[g, i])
                end
            end
        elseif t <= max(0, time_limits[:down] - time_down_t0[g]) && time_down_t0[g] > 0
            JuMP.add_to_expression!(lhs_off, 1)
        else
            continue
        end
        downtime_constraints[g, t] = JuMP.@constraint(jump_model, lhs_off + ug[g, t] <= 1.0)
    end

    # Eq (11) Start-up lag
    startup_lag_constraints = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        thermal_gen_names,
        startup_categories[1:(end - 1)],
        time_steps,
    )
    for g in thermal_gen_names,
        (si, startup) in enumerate(startup_categories[1:(end - 1)]),
        t in time_steps

        g_startup = PSI._convert_hours_to_timesteps(
            get_start_time_limits(get_component(ThermalMultiStart, system, g)),
            resolution,
        )
        if t >= g_startup[si + 1]
            time_range = UnitRange{Int}(Int(g_startup[si]), Int(g_startup[si + 1] - 1))
            startup_lag_constraints[g, startup, t] = JuMP.@constraint(
                jump_model,
                δ_sg[g, startup, t] <= sum(wg[g, t - i] for i in time_range)
            )
        end

        # Initial start-up type, based on Tight and Compact eq (15) rather than pg_lib (7)
        # Xpress version adds extra constraints because it's not enforcing the fix
        if (g_startup[si + 1] - time_down_t0[g]) < t && t < g_startup[si + 1]
            if JuMP.solver_name(jump_model) == "Xpress"
                JuMP.@constraint(jump_model, δ_sg[g, startup, t] <= 0)
            else
                JuMP.fix(δ_sg[g, startup, t], 0.0; force = true)
            end
        end
    end

    # Eq (12) Start-up exclusive categories
    startup_exclusive_constraints = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        vg[g, t] == sum(δ_sg[g, s, t] for s in startup_categories)
    )

    if use_reg

        apply_reg_requirements!(problem,
            reg⁺_device_names,
            reg⁻_device_names,
            storage_reserve_names
        )

        # Eq (20) Reg up response time
        reg⁺_response_constraints = JuMP.@constraint(
            jump_model,
            [g in reg⁺_device_names, t in time_steps],
            reg⁺[g, t] <= L_REG * ramp_up[g]
        )
        # Eq (21) Reg down response time
        reg⁻_response_constraints = JuMP.@constraint(
            jump_model,
            [g in reg⁻_device_names, t in time_steps],
            reg⁻[g, t] <= L_REG * ramp_dn[g]
        )
    end
    if use_spin
        # Eq (19) Total spin
        spin_constraints = JuMP.@constraint(
            jump_model,
            [t in time_steps],
            sum(spin[g, t] for g in spin_device_names) >=
            required_spin[t] - (use_slack ? slack_spin[t] : 0)
        )
        # Eq (22) Spin response time
        spin_response_constraints = JuMP.@constraint(
            jump_model,
            [g in spin_device_names, t in time_steps],
            spin[g, t] <= L_SPIN * ramp_up[g]
        )
    end

    # Eq (30) Power balance constraint, Eq (45) hydro included
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

    # Eq (31) Max output 1 -- in 3 parts for 3 reserve groupings
    if use_spin
        maxoutput1_constraint_spin = JuMP.@constraint(
            jump_model,
            [g in spin_device_names, j in scenarios, t in time_steps],
            pg[g, j, t] + spin[g, t] <=
            (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
            max(0, (pg_lim[g].max - pg_power_trajectory[g].startup)) * vg[g, t]
        )
    end
    if use_reg
        maxoutput1_constraint_reg = JuMP.@constraint(
            jump_model,
            [g in reg⁺_device_names, j in scenarios, t in time_steps],
            pg[g, j, t] + reg⁺[g, t] <=
            (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
            max(0, (pg_lim[g].max - pg_power_trajectory[g].startup)) * vg[g, t]
        )
    end
    maxoutput1_constraint_none = JuMP.@constraint(
        jump_model,
        [g in no_res⁺_device_names, j in scenarios, t in time_steps],
        pg[g, j, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].startup)) * vg[g, t]
    )

    # Limits on downward reserves 
    reservedn_constraint_reg = JuMP.@constraint(
        jump_model,
        [g in reg⁻_device_names, j in scenarios, t in time_steps],
        pg[g, j, t] - reg⁻[g, t] >= 0
    )

    # Eq (32) Max output 2 -- in 3 parts for 3 reserve groupings
    if use_spin
        maxoutput2_constraint_spin = JuMP.@constraint(
            jump_model,
            [g in spin_device_names, j in scenarios, t in time_steps[1:(end - 1)]],
            pg[g, j, t] + spin[g, t] <=
            (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
            max(0, (pg_lim[g].max - (pg_power_trajectory[g].shutdown <= pg_lim[g].min ? pg_lim[g].max : pg_power_trajectory[g].shutdown))) * wg[g, t + 1]
            # max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, t + 1]
        )
    end
    if use_reg
        maxoutput2_constraint_reg = JuMP.@constraint(
            jump_model,
            [g in reg⁺_device_names, j in scenarios, t in time_steps[1:(end - 1)]],
            pg[g, j, t] + reg⁺[g, t] <=
            (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
            max(0, (pg_lim[g].max - (pg_power_trajectory[g].shutdown <= pg_lim[g].min ? pg_lim[g].max : pg_power_trajectory[g].shutdown))) * wg[g, t + 1]
            # max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, t + 1]
        )
    end
    maxoutput2_constraint_none = JuMP.@constraint(
        jump_model,
        [g in no_res⁺_device_names, j in scenarios, t in time_steps[1:(end - 1)]],
        pg[g, j, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - (pg_power_trajectory[g].shutdown <= pg_lim[g].min ? pg_lim[g].max : pg_power_trajectory[g].shutdown))) * wg[g, t + 1]
        # max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, t + 1]
    )
    # Initial condition ignores t0 reserves, is same for all reserve groups. pg_lib (10)
    # Leaving separate because not indexed over j
    maxoutput2_constraint_all1 = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names],
        ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <=
        (pg_lim[g].max - pg_lim[g].min) * ug_t0[g] -
        max(0, (pg_lim[g].max - (pg_power_trajectory[g].shutdown <= pg_lim[g].min ? pg_lim[g].max : pg_power_trajectory[g].shutdown))) * wg[g, 1]
        # max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, 1]
    )

    # Eq (33) Ramp up -- in 3 parts for 3 reserve groupings
    if use_spin
        rampup_constraint_spin = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
            undef,
            spin_device_names,
            scenarios,
            time_steps,
        )
        for g in spin_device_names, j in scenarios, t in time_steps
            if t == 1  # pg_lib (8)
                rampup_constraint_spin[g, j, 1] = JuMP.@constraint(
                    jump_model,
                    pg[g, j, 1] + spin[g, 1] -
                    ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <= ramp_up[g]
                )
            else
                rampup_constraint_spin[g, j, t] = JuMP.@constraint(
                    jump_model,
                    pg[g, j, t] + spin[g, t] - pg[g, j, t - 1] <= ramp_up[g]
                )
            end
        end
    end
    if use_reg
        rampup_constraint_reg = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
            undef,
            reg⁺_device_names,
            scenarios,
            time_steps,
        )
        for g in reg⁺_device_names, j in scenarios, t in time_steps
            if t == 1  # pg_lib (8)
                rampup_constraint_reg[g, j, 1] = JuMP.@constraint(
                    jump_model,
                    pg[g, j, 1] + reg⁺[g, 1] - ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <=
                    ramp_up[g]
                )
            else
                rampup_constraint_reg[g, j, t] = JuMP.@constraint(
                    jump_model,
                    pg[g, j, t] + reg⁺[g, t] - pg[g, j, t - 1] <= ramp_up[g]
                )
            end
        end
    end
    rampup_constraint_none = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        no_res⁺_device_names,
        scenarios,
        time_steps,
    )
    for g in no_res⁺_device_names, j in scenarios, t in time_steps
        if t == 1  # pg_lib (8)
            rampup_constraint_none[g, j, 1] = JuMP.@constraint(
                jump_model,
                pg[g, j, 1] - ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <= ramp_up[g]
            )
        else
            rampup_constraint_none[g, j, t] =
                JuMP.@constraint(jump_model, pg[g, j, t] - pg[g, j, t - 1] <= ramp_up[g])
        end
    end

    # Eq (34) Ramp down -- in 2 parts for 2 reserve groupings
    if use_reg
        rampdn_constraint_reg = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
            undef,
            reg⁻_device_names,
            scenarios,
            time_steps,
        )
        for g in reg⁻_device_names, j in scenarios, t in time_steps
            if t == 1  # pg_lib (9)
                rampdn_constraint_reg[g, j, 1] = JuMP.@constraint(
                    jump_model,
                    ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) - pg[g, j, 1] - reg⁻[g, t] <=
                    ramp_dn[g]
                )
            else
                rampdn_constraint_reg[g, j, t] = JuMP.@constraint(
                    jump_model,
                    pg[g, j, t - 1] - pg[g, j, t] - reg⁻[g, t] <= ramp_dn[g]
                )
            end
        end
    end
    rampdn_constraint_none = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        no_res⁻_device_names,
        scenarios,
        time_steps,
    )
    for g in no_res⁻_device_names, j in scenarios, t in time_steps
        if t == 1  # pg_lib (9)
            rampdn_constraint_none[g, j, 1] = JuMP.@constraint(
                jump_model,
                ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) - pg[g, j, 1] <= ramp_dn[g]
            )
        else
            rampdn_constraint_none[g, j, t] =
                JuMP.@constraint(jump_model, pg[g, j, t - 1] - pg[g, j, t] <= ramp_dn[g])
        end
    end
    
    # Apply CC constraints
    restrictions = problem.ext["cc_restrictions"]
    cc_constraints = JuMP.@constraint(
        jump_model,
        [k in keys(restrictions), t in time_steps],
        sum(ug[i, t] for i in restrictions[k]) <= 1
    )
end