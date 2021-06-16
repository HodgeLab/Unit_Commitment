
function apply_thermal_constraints!(
        problem::PSI.OperationsProblem{T},
        spin_device_names::Vector{String}
    ) where T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC, BasecaseUnitCommitmentCC}
    use_reg = problem.ext["use_reg"]
    use_spin = problem.ext["use_spin"]
    use_must_run = problem.ext["use_must_run"]
    L_REG = problem.ext["L_REG"]
    L_SPIN = problem.ext["L_SPIN"]
    allowable_reserve_prop = problem.ext["allowable_reserve_prop"]

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    resolution = PSI.model_resolution(optimization_container)

    # ------------------------------------------------
    # Input data
    # ------------------------------------------------
    MINS_IN_HOUR = 60.0
    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    pg_lim = Dict(
        g => get_active_power_limits(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
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
    if use_must_run
        must_run_gen_names =
            get_name.(get_components(ThermalMultiStart, system, x -> PSY.get_must_run(x)))
    end
    startup_categories = (:hot, :warm, :cold)
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

    # ------------------------------------------------
    # Collect variables
    # ------------------------------------------------
    pg = jump_model.obj_dict[:pg]
    ug = jump_model.obj_dict[:ug]
    wg = jump_model.obj_dict[:wg]
    vg = jump_model.obj_dict[:vg]
    δ_sg = jump_model.obj_dict[:δ_sg]
    total_reserve⁺ = optimization_container.expressions[:total_reserve⁺]
    total_reserve⁻ = optimization_container.expressions[:total_reserve⁻]

    # Commitment constraints
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

    # Apply CC constraints
    restrictions = problem.ext["cc_restrictions"]
    if :ng in keys(jump_model.obj_dict)
        ng = jump_model.obj_dict[:ng]
        JuMP.@constraint(
            jump_model,
            [k in keys(restrictions), j in scenarios, t in time_steps],
            sum(ug[i, t] + ng[i, j, t] for i in restrictions[k]) <= 1
        )
    else
        cc_constraints = JuMP.@constraint(
            jump_model,
            [k in keys(restrictions), t in time_steps],
            sum(ug[i, t] for i in restrictions[k]) <= 1
        )
    end

    # Must-run -- added as a constraint for Xpress due to fix not working
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

    # Up and down time constraints
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

    # Start-up lag
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

    # Start-up exclusive categories
    startup_exclusive_constraints = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        vg[g, t] == sum(δ_sg[g, s, t] for s in startup_categories)
    )

    # if use_reg
    #     reg⁺ = jump_model.obj_dict[:reg⁺]
    #     reg⁻ = jump_model.obj_dict[:reg⁻]
    #     # Reg up response time
    #     reg⁺_response_constraints = JuMP.@constraint(
    #         jump_model,
    #         [g in reg⁺_device_names, t in time_steps],
    #         reg⁺[g, t] <= L_REG * ramp_up[g]
    #     )
    #     # Reg down response time
    #     reg⁻_response_constraints = JuMP.@constraint(
    #         jump_model,
    #         [g in reg⁻_device_names, t in time_steps],
    #         reg⁻[g, t] <= L_REG * ramp_dn[g]
    #     )
    # end

    if use_spin
        spin = jump_model.obj_dict[:spin]
        # Spin response time
        spin_response_constraints = JuMP.@constraint(
            jump_model,
            [g in spin_device_names, t in time_steps],
            spin[g, t] <= L_SPIN * ramp_up[g]
        )
    end

    # Restrict reserve allocation by % rather than response time
    JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        total_reserve⁺[g, t] <= pg_lim[g].max * allowable_reserve_prop
    )
    JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        total_reserve⁻[g, t] <= pg_lim[g].max * allowable_reserve_prop
    )

    _apply_thermal_scenario_based_constraints!(problem,
        pg_lim,
        ramp_up,
        ramp_dn,
        ug_t0)

    return
end


function _apply_thermal_scenario_based_constraints!(
        problem::PSI.OperationsProblem{T},
        pg_lim::Dict{},
        ramp_up::Dict{},
        ramp_dn::Dict{},
        ug_t0::Dict{}
    ) where T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC}
    use_reg = problem.ext["use_reg"]
    use_spin = problem.ext["use_spin"]

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    scenarios = 1:size(jump_model.obj_dict[:pS])[1]

    # ------------------------------------------------
    # Input data
    # ------------------------------------------------
    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    pg_power_trajectory = Dict(
        g => get_power_trajectory(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    # This is just power (Pg), not power above minimum (pg)
    Pg_t0 = Dict(
        g => get_active_power(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    variable_cost = Dict(
        g => get_variable(get_operation_cost(get_component(ThermalMultiStart, system, g))) for g in thermal_gen_names
    )

    # ------------------------------------------------
    # Collect or define variables
    # ------------------------------------------------
    pg = jump_model.obj_dict[:pg]
    ug = jump_model.obj_dict[:ug]
    wg = jump_model.obj_dict[:wg]
    vg = jump_model.obj_dict[:vg]
    λ = JuMP.@variable(
        jump_model,
        0 <=
        λ[g in thermal_gen_names, j in scenarios, i in 1:length(variable_cost[g]), t in time_steps] <=
        PSY.get_breakpoint_upperbounds(variable_cost[g])[i]
    )
    total_reserve⁺ = optimization_container.expressions[:total_reserve⁺]
    total_reserve⁻ = optimization_container.expressions[:total_reserve⁻]

    # PWL variable cost constraint
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

    # Max output 1
    maxoutput1_constraint = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, j in scenarios, t in time_steps],
        pg[g, j, t] + total_reserve⁺[g, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].startup)) * vg[g, t]
    )

    # Limits on downward reserves 
    reservedn_constraint_reg = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, j in scenarios, t in time_steps],
        pg[g, j, t] - total_reserve⁻[g, t] >= 0
    )

    # Max output 2
    maxoutput2_constraint = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, j in scenarios, t in time_steps[1:(end - 1)]],
        pg[g, j, t] + total_reserve⁺[g, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - (pg_power_trajectory[g].shutdown <= pg_lim[g].min ? pg_lim[g].max : pg_power_trajectory[g].shutdown))) * wg[g, t + 1]
        # max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, t + 1]
    )
    # Initial condition ignores t0 reserves, is same for all reserve groups. pg_lib (10)
    # Leaving separate because not indexed over j
    maxoutput2_constraint_1 = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names],
        ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <=
        (pg_lim[g].max - pg_lim[g].min) * ug_t0[g] -
        max(0, (pg_lim[g].max - (pg_power_trajectory[g].shutdown <= pg_lim[g].min ? pg_lim[g].max : pg_power_trajectory[g].shutdown))) * wg[g, 1]
        # max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, 1]
    )

    # Ramp up
    rampup_constraint = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        thermal_gen_names,
        scenarios,
        time_steps,
    )
    for g in thermal_gen_names, j in scenarios, t in time_steps
        if t == 1  # pg_lib (8)
            rampup_constraint[g, j, 1] = JuMP.@constraint(
                jump_model,
                pg[g, j, 1] + total_reserve⁺[g, 1] - ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <=
                ramp_up[g]
            )
        else
            rampup_constraint[g, j, t] = JuMP.@constraint(
                jump_model,
                pg[g, j, t] + total_reserve⁺[g, t] - pg[g, j, t - 1] <= ramp_up[g]
            )
        end
    end

    # Ramp down
    rampdn_constraint = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        thermal_gen_names,
        scenarios,
        time_steps,
    )
    for g in thermal_gen_names, j in scenarios, t in time_steps
        if t == 1  # pg_lib (9)
            rampdn_constraint[g, j, 1] = JuMP.@constraint(
                jump_model,
                ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) - pg[g, j, 1] - total_reserve⁻[g, 1] <=
                ramp_dn[g]
            )
        else
            rampdn_constraint[g, j, t] = JuMP.@constraint(
                jump_model,
                pg[g, j, t - 1] - pg[g, j, t] - total_reserve⁻[g, t] <= ramp_dn[g]
            )
        end
    end

    if :supp in keys(jump_model.obj_dict)
        supp = jump_model.obj_dict[:supp]
        ng = jump_model.obj_dict[:ng]
        scenarios = 1:size(supp)[2]
        L_SUPP = problem.ext["L_SUPP"]

        # Non-spin constraints
        # Supplemental reserve must exceed minimum power rating
        JuMP.@constraint(
            jump_model,
            [g in thermal_gen_names, j in scenarios, t in time_steps],
            supp[g, j, t] >= pg_lim[g].min * ng[g, j, t]
        )

        # Supplemental only as high as portion of start-up it can achieve in 10 minutes
        JuMP.@constraint(
            jump_model,
            [g in thermal_gen_names, j in scenarios, t in time_steps],
            supp[g, j, t] <= pg_power_trajectory[g].startup * L_SUPP * ng[g, j, t]
        )

        # Supplemental only as high as unit can ramp in 10 minutes
        JuMP.@constraint(
            jump_model,
            [g in thermal_gen_names, j in scenarios, t in time_steps],
            supp[g, j, t] <=  L_SUPP * ramp_up[g]
        )

    end

    return
end

function _apply_thermal_scenario_based_constraints!(
        problem::PSI.OperationsProblem{T},
        pg_lim::Dict{},
        ramp_up::Dict{},
        ramp_dn::Dict{},
        ug_t0::Dict{}
    ) where T <: BasecaseUnitCommitmentCC
    use_reg = problem.ext["use_reg"]
    use_spin = problem.ext["use_spin"]

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)

    # ------------------------------------------------
    # Input data
    # ------------------------------------------------
    thermal_gen_names = get_name.(get_components(ThermalMultiStart, system))
    pg_power_trajectory = Dict(
        g => get_power_trajectory(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    # This is just power (Pg), not power above minimum (pg)
    Pg_t0 = Dict(
        g => get_active_power(get_component(ThermalMultiStart, system, g)) for
        g in thermal_gen_names
    )
    variable_cost = Dict(
        g => get_variable(get_operation_cost(get_component(ThermalMultiStart, system, g))) for g in thermal_gen_names
    )

    # ------------------------------------------------
    # Collect or define variables
    # ------------------------------------------------
    pg = jump_model.obj_dict[:pg]
    ug = jump_model.obj_dict[:ug]
    wg = jump_model.obj_dict[:wg]
    vg = jump_model.obj_dict[:vg]
    λ = JuMP.@variable(
        jump_model,
        0 <=
        λ[g in thermal_gen_names, i in 1:length(variable_cost[g]), t in time_steps] <=
        PSY.get_breakpoint_upperbounds(variable_cost[g])[i]
    )
    total_reserve⁺ = optimization_container.expressions[:total_reserve⁺]
    total_reserve⁻ = optimization_container.expressions[:total_reserve⁻]

    # PWL variable cost constraint
    # PWL Cost function auxiliary variables
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
    optimization_container.expressions[:Cg] = Cg

    # Max output 1
    maxoutput1_constraint = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, t in time_steps],
        pg[g, t] + total_reserve⁺[g, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - pg_power_trajectory[g].startup)) * vg[g, t]
    )

    # Limits on downward reserves 
    if use_reg
        reservedn_constraint_reg = JuMP.@constraint(
            jump_model,
            [g in thermal_gen_names, t in time_steps],
            pg[g, t] - total_reserve⁻[g, t] >= 0
        )
    end

    # Max output 2
    maxoutput2_constraint_reg = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names, t in time_steps[1:(end - 1)]],
        pg[g, t] + total_reserve⁺[g, t] <=
        (pg_lim[g].max - pg_lim[g].min) * ug[g, t] -
        max(0, (pg_lim[g].max - (pg_power_trajectory[g].shutdown <= pg_lim[g].min ? pg_lim[g].max : pg_power_trajectory[g].shutdown))) * wg[g, t + 1]
        # max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, t + 1]
    )
    # Initial condition ignores t0 reserves, is same for all reserve groups. pg_lib (10)
    maxoutput2_constraint_all1 = JuMP.@constraint(
        jump_model,
        [g in thermal_gen_names],
        ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <=
        (pg_lim[g].max - pg_lim[g].min) * ug_t0[g] -
        max(0, (pg_lim[g].max - (pg_power_trajectory[g].shutdown <= pg_lim[g].min ? pg_lim[g].max : pg_power_trajectory[g].shutdown))) * wg[g, 1]
        # max(0, (pg_lim[g].max - pg_power_trajectory[g].shutdown)) * wg[g, 1]
    )

    # Ramp up -- in 3 parts for 3 reserve groupings
    rampup_constraint = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        thermal_gen_names,
        time_steps,
    )
    for g in thermal_gen_names, t in time_steps
        if t == 1  # pg_lib (8)
            rampup_constraint[g, 1] = JuMP.@constraint(
                jump_model,
                pg[g, 1] + total_reserve⁺[g, 1] -
                ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) <= ramp_up[g]
            )
        else
            rampup_constraint[g, t] = JuMP.@constraint(
                jump_model,
                pg[g, t] + total_reserve⁺[g, t] - pg[g, t - 1] <= ramp_up[g]
            )
        end
    end

    # Ramp down
    rampdn_constraint = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        thermal_gen_names,
        time_steps,
    )
    for g in thermal_gen_names, t in time_steps
        if t == 1  # pg_lib (9)
            rampdn_constraint[g, 1] = JuMP.@constraint(
                jump_model,
                ug_t0[g] * (Pg_t0[g] - pg_lim[g].min) - pg[g, 1] - total_reserve⁻[g, 1] <=
                ramp_dn[g]
            )
        else
            rampdn_constraint[g, t] = JuMP.@constraint(
                jump_model,
                pg[g, t - 1] - pg[g, t] - total_reserve⁻[g, t] <= ramp_dn[g]
            )
        end
    end

    # Linked reserve response times
    # Does not need to be implemented while reg⁺ and spin groups are mutually exclusive

    return
end
