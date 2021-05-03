struct CVaRUnitCommitmentCC <: PSI.PowerSimulationsOperationsProblem end

function PSI.problem_build!(problem::PSI.OperationsProblem{CVaRUnitCommitmentCC})
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

    # Sets
    scenarios = 1:31

    # Constants
    MINS_IN_HOUR = 60.0
    Δt = 1
    L_REG = 1 / 12 # 5 min
    L_SPIN = 1 / 6 # 10 min
    L_SUPP = 1 / 4 # 15 min, to start
    C_RR = 1000 # Penalty cost of recourse reserve
    α = 0.20 # Risk tolerance level

    # -------------------------------------------------------------
    # Collect definitions from PSY model
    # -------------------------------------------------------------
    thermal_gens = get_components(ThermalMultiStart, system)
    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    hydro_generators = PSY.get_components(PSY.HydroGen, system)
    hydro_gens_names = PSY.get_name.(hydro_generators)
    get_rmp_up_limit(g) = PSY.get_ramp_limits(g).up
    get_rmp_dn_limit(g) = PSY.get_ramp_limits(g).down
    ramp_up = Dict(
        g => get_rmp_up_limit(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    ramp_dn = Dict(
        g => get_rmp_dn_limit(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    pg_lim = Dict(
        g => get_active_power_limits(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    pg_power_trajectory = Dict(
        g => get_power_trajectory(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    must_run_gen_names =
        get_name.(get_components(ThermalMultiStart, system, x -> PSY.get_must_run(x)))
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

    required_reg⁺ =
        get_requirement(reg_reserve_up) * get_time_series_values(
            Deterministic,
            reg_reserve_up,
            "requirement";
            start_time = case_initial_time,
        )
    required_reg⁻ =
        get_requirement(reg_reserve_dn) * get_time_series_values(
            Deterministic,
            reg_reserve_dn,
            "requirement";
            start_time = case_initial_time,
        )
    required_spin =
        get_requirement(spin_reserve) * get_time_series_values(
            Deterministic,
            spin_reserve,
            "requirement";
            start_time = case_initial_time,
        )

    # -------------------------------------------------------------
    # Time-series data
    # -------------------------------------------------------------
    total_load = zeros(length(time_steps))
    total_wind = zeros(length(time_steps))
    total_hydro = zeros(length(time_steps))

    for l in get_components(PowerLoad, system)
        f = get_time_series_values(
            Deterministic,
            l,
            "max_active_power";
            start_time = case_initial_time,
        )
        total_load .+= f * PSY.get_max_active_power(l)
    end

    for l in
        get_components(RenewableGen, system, x -> get_prime_mover(x) != PrimeMovers.PVe)
        f = get_time_series_values(
            Deterministic,
            l,
            "max_active_power";
            start_time = case_initial_time,
        )
        total_wind .+= f * PSY.get_max_active_power(l)
    end

    for l in get_components(HydroGen, system)
        f = get_time_series_values(
            Deterministic,
            l,
            "max_active_power";
            start_time = case_initial_time,
        )
        total_hydro .+= f * PSY.get_max_active_power(l)
    end

    # Populate solar scenarios
    area = PSY.get_component(Area, system, "1")
    area_solar_forecast_scenarios = ones(31, length(time_steps))
    #area_solar_forecast_scenarios = PSY.get_time_series_values(
    #            Scenarios,
    #            area,
    #            "solar_power";
    #            start_time = case_initial_time,
    #) / 100

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
    pg = JuMP.@variable(jump_model, pg[g in thermal_gen_names, t in time_steps] >= 0) # power ABOVE MINIMUM
    pS = JuMP.@variable(jump_model, pS[j in scenarios, t in time_steps] >= 0)
    pW = JuMP.@variable(jump_model, pW[t in time_steps] >= 0)
    reg⁺ = JuMP.@variable(jump_model, reg⁺[g in reg⁺_device_names, t in time_steps] >= 0)
    reg⁻ = JuMP.@variable(jump_model, reg⁻[g in reg⁻_device_names, t in time_steps] >= 0)
    spin = JuMP.@variable(jump_model, spin[g in spin_device_names, t in time_steps] >= 0)
    # Come back to decide contributing device set if needed
    supp⁺ = JuMP.@variable(
        jump_model,
        supp⁺[g in spin_device_names, j in scenarios, t in time_steps] >= 0
    )
    supp⁻ = JuMP.@variable(
        jump_model,
        supp⁻[g in spin_device_names, j in scenarios, t in time_steps] >= 0
    )
    total_supp⁺ =
        JuMP.@variable(jump_model, total_supp⁺[j in scenarios, t in time_steps] >= 0)
    total_supp⁻ =
        JuMP.@variable(jump_model, total_supp⁻[j in scenarios, t in time_steps] >= 0)
    z = JuMP.@variable(jump_model, z[j in scenarios] >= 0) # Eq (25)
    β = JuMP.@variable(jump_model, β)
    λ = JuMP.@variable(
        jump_model,
        0 <=
        λ[g in thermal_gen_names, i in 1:length(variable_cost[g]), t in time_steps] <=
        PSY.get_breakpoint_upperbounds(variable_cost[g])[i]
    )
    # -------------------------------------------------------------
    # Constraints
    # -------------------------------------------------------------

    # Eq (5) Wind constraint
    wind_constraint =
        JuMP.@constraint(jump_model, [t in time_steps], pW[t] <= total_wind[t])

    # Eq (6) PWL variable cost constraint
    # PWL Cost function auxiliary variables
    lambda_bound = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, i in 1:length(variable_cost[g]), t in time_steps],
        0 <= λ[g, i, t] <= PSY.get_breakpoint_upperbounds(variable_cost[g])[i]
    )

    lambda_constraint = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        sum(λ[g, i, t] for i in 1:length(variable_cost[g])) ==
        pg_lim[g].min * ug[g, t] + pg[g, t]
    )

    Cg = JuMP.@expression(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        sum(
            100.0 * λ[g, i, t] * PSY.get_slopes(variable_cost[g])[i] for
            i in 1:length(variable_cost[g])
        )
    )

    obj_function = JuMP.@objective(
        jump_model,
        Min,
        sum(
            no_load_cost[g] * ug[g, t] +
            Cg[g, t] +
            sum(startup_cost[g][s] * δ_sg[g, s, t] for s in startup_categories) +
            shutdown_cost[g] * wg[g, t] for g in thermal_gen_names, t in time_steps
        ) +
        C_RR * (β + 1 / (length(scenarios) * (1 - α)) * sum(z[j] for j in scenarios))
    )

    # Eq (7) Commitment constraints
    commitment_constraints = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        thermal_gen_names,
        time_steps,
    )
    for g in thermal_gen_names, t in time_steps
        if t == 1  # pg_lib (6)
            commitment_constraints[g, 1] = JuMP.@constraint(
                jump_model,
                ug[g, 1] - ug_t0[g] == vg[g, 1] - wg[g, 1]
            )
        else
            commitment_constraints[g, t] = JuMP.@constraint(
                jump_model,
                ug[g, t] - ug[g, t - 1] == vg[g, t] - wg[g, t]
            )
        end
    end

    # Eq (8) Must-run
    for g in must_run_gen_names, t in time_steps
        JuMP.fix(ug[g, t], 1.0)
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
    for g in thermal_gen_names, (si, startup) in enumerate(startup_categories[1:(end - 1)]), t in time_steps
        g_startup = get_start_time_limits(get_component(ThermalMultiStart, system, g))
        if t >= g_startup[si + 1]
            time_range = UnitRange{Int}(
                Int(ceil(g_startup[si])),
                Int(ceil(g_startup[si + 1] - 1))
            )
            startup_lag_constraints[g, startup, t] = JuMP.@constraint(
                jump_model,
                δ_sg[g, startup, t] <= sum(wg[g, t - i] for i in time_range)
            )
        end
        # Initial start-up type, based on Tight and Compact eq (15) rather than pg_lib (7)
        if (g_startup[si + 1] - time_down_t0[g]) < t && t < g_startup[si + 1]
            JuMP.fix(δ_sg[g, startup, t], 0.0)
        end
    end

    # Eq (12) Start-up exclusive categories
    startup_exclusive_constraints = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        vg[g, t] == sum(δ_sg[g, s, t] for s in startup_categories)
    )

    # Eq (17) Total reg up
    reg⁺_constraints = JuMP.@constraint(
        jump_model,
        [t in time_steps],
        sum(reg⁺[g, t] for g in reg⁺_device_names) >= required_reg⁺[t]
    )
    # Eq (18) Total reg down
    reg⁻_constraints = JuMP.@constraint(
        jump_model,
        [t in time_steps],
        sum(reg⁻[g, t] for g in reg⁻_device_names) >= required_reg⁻[t]
    )
    # Eq (19) Total spin
    spin_constraints = JuMP.@constraint(
        jump_model,
        [t in time_steps],
        sum(spin[g, t] for g in spin_device_names) >= required_spin[t]
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
    # Eq (22) Spin response time
    spin_response_constraints = JuMP.@constraint(
        jump_model,
        [g in spin_device_names, t in time_steps],
        spin[g, t] <= L_SPIN * ramp_up[g]
    )

    # Eq (23) Solar scenarios
    solar_constraints = JuMP.@constraint(
        jump_model,
        [j in scenarios, t in time_steps],
        pS[j, t] <= area_solar_forecast_scenarios[j, t]
    )

    # Eq (25) is included in z variable definition
    # Eq (26) Auxiliary variable definition
    auxiliary_constraint = JuMP.@constraint(
        jump_model,
        [j in scenarios],
        z[j] >= sum(total_supp⁺[j, t] + total_supp⁻[j, t] for t in time_steps) - β
    )

    # Eq (27) Total supplemental up defintion
    supp⁺_constraint = JuMP.@constraint(
        jump_model,
        [j in scenarios, t in time_steps],
        sum(supp⁺[g, j, t] for g in spin_device_names) >= total_supp⁺[j, t]
    )
    # Eq (28) Total supplemental down defintion
    supp⁻_constraint = JuMP.@constraint(
        jump_model,
        [j in scenarios, t in time_steps],
        sum(supp⁻[g, j, t] for g in spin_device_names) >= total_supp⁻[j, t]
    )
    # Eq (29) is included in supp_ variable definitions

    # Eq (30) Power balance constraint, Eq (45) hydro included
    power_balance_constraint = JuMP.@constraint(
        jump_model,
        [j in scenarios, t in time_steps],
        sum(pg[g, t] + pg_lim[g].min * ug[g, t] for g in thermal_gen_names) +
        pS[j, t] +
        pW[t] +
        total_supp⁺[j, t] - total_supp⁻[j, t] + total_hydro[t] == total_load[t]
    )

    # Eq (31) Max output 1 -- in 3 parts for 3 reserve groupings
    maxoutput1_constraint_spin = JuMP.@constraint(
        jump_model,
        [g in spin_device_names, j in scenarios, t in time_steps],
        pg[g, t] + spin[g, t] + supp⁺[g, j, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].startup)) * vg[g, t]
    )
    maxoutput1_constraint_reg = JuMP.@constraint(
        jump_model,
        [g in reg⁺_device_names, t in time_steps],
        pg[g, t] + reg⁺[g, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].startup)) * vg[g, t]
    )
    maxoutput1_constraint_none = JuMP.@constraint(
        jump_model,
        [
            g in setdiff(thermal_gen_names, union(reg⁺_device_names, spin_device_names)),
            t in time_steps,
        ],
        pg[g, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].startup)) * vg[g, t]
    )

    # Eq (32) Max output 2 -- in 3 parts for 3 reserve groupings
    maxoutput2_constraint_spin = JuMP.@constraint(
        jump_model,
        [g in spin_device_names, j in scenarios, t in time_steps[1:(end - 1)]],
        pg[g, t] + spin[g, t] + supp⁺[g, j, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, t + 1]
    )
    maxoutput2_constraint_reg = JuMP.@constraint(
        jump_model,
        [g in reg⁺_device_names, t in time_steps[1:(end - 1)]],
        pg[g, t] + reg⁺[g, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, t + 1]
    )
    maxoutput2_constraint_none = JuMP.@constraint(
        jump_model,
        [
            g in setdiff(thermal_gen_names, union(reg⁺_device_names, spin_device_names)),
            t in time_steps[1:(end - 1)],
        ],
        pg[g, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, t + 1]
    )
    # Initial condition ignores t0 reserves, is same for all reserve groups. pg_lib (10)
    # Leaving separate because not indexed over j
    maxoutput2_constraint_all1 = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names],
        ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <=
        (pg_lim[g].max - pg_lim[g].min) * ug_t0[g] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, 1]
    )

    # Eq (33) Ramp up -- in 3 parts for 3 reserve groupings
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
                pg[g, 1] + spin[g, 1] + supp⁺[g, j, 1] -
                ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <= ramp_up[g]
            )
        else
            rampup_constraint_spin[g, j, t] = JuMP.@constraint(
                jump_model,
                pg[g, t] + spin[g, t] + supp⁺[g, j, t] - pg[g, t - 1] <= ramp_up[g]
            )
        end
    end
    rampup_constraint_reg = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        reg⁺_device_names,
        time_steps,
    )
    for g in reg⁺_device_names, t in time_steps
        if t == 1  # pg_lib (8)
            rampup_constraint_reg[g, 1] = JuMP.@constraint(
                jump_model,
                pg[g, 1] + reg⁺[g, 1] - ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <= ramp_up[g]
            )
        else
            rampup_constraint_reg[g, t] = JuMP.@constraint(
                jump_model,
                pg[g, t] + reg⁺[g, t] - pg[g, t - 1] <= ramp_up[g]
            )
        end
    end
    rampup_constraint_none = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        setdiff(thermal_gen_names, union(reg⁺_device_names, spin_device_names)),
        time_steps,
    )
    for g in setdiff(thermal_gen_names, union(reg⁺_device_names, spin_device_names)),
        t in time_steps

        if t == 1  # pg_lib (8)
            rampup_constraint_none[g, 1] = JuMP.@constraint(
                jump_model,
                pg[g, 1] - ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <= ramp_up[g]
            )
        else
            rampup_constraint_none[g, t] =
                JuMP.@constraint(jump_model, pg[g, t] - pg[g, t - 1] <= ramp_up[g])
        end
    end

    # Eq (34) Ramp down -- in 3 parts for 3 reserve groupings
    rampdn_constraint_supp = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        spin_device_names,
        scenarios,
        time_steps,
    )
    for g in spin_device_names, j in scenarios, t in time_steps
        if t == 1  # pg_lib (9)
            rampdn_constraint_supp[g, j, 1] = JuMP.@constraint(
                jump_model,
                ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) - pg[g, 1] - supp⁻[g, j, 1] <=
                ramp_dn[g]
            )
        else
            rampdn_constraint_supp[g, j, t] = JuMP.@constraint(
                jump_model,
                pg[g, t - 1] - pg[g, t] - supp⁻[g, j, t] <= ramp_dn[g]
            )
        end
    end
    rampdn_constraint_reg = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        reg⁻_device_names,
        time_steps,
    )
    for g in reg⁻_device_names, t in time_steps
        if t == 1  # pg_lib (9)
            rampdn_constraint_reg[g, 1] = JuMP.@constraint(
                jump_model,
                ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) - pg[g, 1] - reg⁻[g, t] <= ramp_dn[g]
            )
        else
            rampdn_constraint_reg[g, t] = JuMP.@constraint(
                jump_model,
                pg[g, t - 1] - pg[g, t] - reg⁻[g, t] <= ramp_dn[g]
            )
        end
    end
    rampdn_constraint_none = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        setdiff(thermal_gen_names, union(reg⁻_device_names, spin_device_names)),
        time_steps,
    )
    for g in setdiff(thermal_gen_names, union(reg⁻_device_names, spin_device_names)),
        t in time_steps

        if t == 1  # pg_lib (9)
            rampdn_constraint_none[g, 1] = JuMP.@constraint(
                jump_model,
                ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) - pg[g, 1] <= ramp_dn[g]
            )
        else
            rampdn_constraint_none[g, t] =
                JuMP.@constraint(jump_model, pg[g, t - 1] - pg[g, t] <= ramp_dn[g])
        end
    end

    # Eq (35) Supp reserve response time
    supp⁺_response_constraints = JuMP.@constraint(
        jump_model,
        [g in spin_device_names, j in scenarios, t in time_steps],
        supp⁺[g, j, t] <= L_SUPP * ramp_up[g]
    )

    # Eq (36) Supp reserve response time
    supp⁻_response_constraints = JuMP.@constraint(
        jump_model,
        [g in spin_device_names, j in scenarios, t in time_steps],
        supp⁻[g, j, t] <= L_SUPP * ramp_dn[g]
    )

    # Eq (37) Linked up-reserve response times
    # Only needs to implement the spin and supp⁺ linkage, while reg⁺ group is mutually exclusive
    linked_⁺_response_constraints = JuMP.@constraint(
        jump_model,
        [g in spin_device_names, j in scenarios, t in time_steps],
        L_SUPP / L_SPIN * spin[g, t] + supp⁺[g, j, t] <= L_SUPP * ramp_up[g]
    )

    # Eq (38) Linked down-reserve response times
    # Does not need to be implemented while reg⁻ and supp⁻ groups are mutually exclusive

    # TODO storage equations

    # Apply CC constraints
    restrictions = problem.ext["cc_restrictions"]
    cc_constraints = JuMP.@constraint(
        jump_model,
        [k in keys(restrictions), t in time_steps],
        sum(ug[i, t] for i in restrictions[k]) <= 1
    )
end
