
function apply_storage!(problem::PSI.OperationsProblem{T},
    storage_reserve_names::Vector{String}
    ) where T
    use_storage_reserves = problem.ext["use_storage_reserves"]
    use_reg = false # Hack for now to ignore storage for reg
    use_spin = problem.ext["use_spin"]
    L_REG = problem.ext["L_REG"]
    L_SPIN = problem.ext["L_SPIN"]

    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)
    system = PSI.get_system(problem)

    storage_names = PSY.get_name.(get_components(PSY.GenericBattery, system))

    # Battery parameters
    get_scaled_SOC = function(b)
        min = get_state_of_charge_limits(b)[:min] * problem.ext["storage_scale"]
        max = get_state_of_charge_limits(b)[:max] * problem.ext["storage_scale"]
        return (min = min, max = max)
    end
    eb_lim = Dict(
        b => get_scaled_SOC(get_component(GenericBattery, system, b)) for
        b in storage_names
    )
    eb_t0 = Dict(
        b => get_initial_energy(get_component(GenericBattery, system, b)) for
        b in storage_names
    )
    η = Dict(
        b => get_efficiency(get_component(GenericBattery, system, b)) for
        b in storage_names
    )
    pb_in_max = Dict(
        b =>
            get_input_active_power_limits(get_component(GenericBattery, system, b))[:max] *
            problem.ext["storage_scale"]
        for b in storage_names
    )
    pb_out_max = Dict(
        b =>
            get_output_active_power_limits(get_component(GenericBattery, system, b))[:max] *
            problem.ext["storage_scale"]
        for b in storage_names
    )

    # Register storage variables
    # 1==discharge, 0==charge
    ϕb = JuMP.@variable(
        jump_model,
        ϕb[b in storage_names, t in time_steps],
        binary = true
    )
    pb_in = JuMP.@variable(jump_model, pb_in[b in storage_names, t in time_steps] >= 0)
    pb_out =
        JuMP.@variable(jump_model, pb_out[b in storage_names, t in time_steps] >= 0)
    eb = JuMP.@variable(
        jump_model,
        eb_lim[b].min <= eb[b in storage_names, t in time_steps] <= eb_lim[b].max
    )

    # Get pre-registered Variables
    reg⁺ = jump_model.obj_dict[:reg⁺]
    reg⁻ = jump_model.obj_dict[:reg⁻]
    spin = jump_model.obj_dict[:spin]

    # Storage charge/discharge decisions
    storage_charge_constraints = JuMP.@constraint(
        jump_model,
        [b in storage_names, t in time_steps],
        pb_in[b, t] <= pb_in_max[b] * (1 - ϕb[b, t])
    )
    storage_discharge_constraints = JuMP.@constraint(
        jump_model,
        [b in storage_names, t in time_steps],
        pb_out[b, t] <= pb_out_max[b] * ϕb[b, t]
    )
    # Storage energy update
    storage_energy_balance = JuMP.Containers.DenseAxisArray{JuMP.ConstraintRef}(
        undef,
        storage_names,
        time_steps,
    )
    for b in storage_names, t in time_steps
        if t == 1
            storage_energy_balance[b, 1] = JuMP.@constraint(
                jump_model,
                eb[b, 1] ==
                eb_t0[b] + η[b].in * pb_in[b, 1] - (1 / η[b].out) * pb_out[b, 1]
            )
        else
            storage_energy_balance[b, t] = JuMP.@constraint(
                jump_model,
                eb[b, t] ==
                eb[b, t - 1] + η[b].in * pb_in[b, t] - (1 / η[b].out) * pb_out[b, t]
            )
        end
    end

    if use_storage_reserves
        # Storage energy satisfies reserve deployment period
        storage_⁺_response_constraints = JuMP.@constraint(
            jump_model,
            [b in storage_reserve_names, t in time_steps],
            η[b].out * (eb[b, t] - eb_lim[b].min) >=
            (use_reg ? L_REG * reg⁺[b, t] : 0) +
            (use_spin ? L_SPIN * spin[b, t] : 0)
        )
        if use_reg
            storage_⁻_response_constraints = JuMP.@constraint(
                jump_model,
                [b in storage_reserve_names, t in time_steps],
                (1 / η[b].in) * (eb_lim[b].max - eb[b, t]) >= L_REG * reg⁻[b, t]
            )
        end

        # Power limits on reserves storage can provide
        storage_⁺_reserve_constraints = JuMP.@constraint(
            jump_model,
            [b in storage_reserve_names, t in time_steps],
            (use_reg ? reg⁺[b, t] : 0) +
            (use_spin ? spin[b, t] : 0) <=
            pb_out_max[b] - pb_out[b, t] + pb_in[b, t]
        )
        if use_reg
            storage_⁻_reserve_constraints = JuMP.@constraint(
                jump_model,
                [b in storage_reserve_names, t in time_steps],
                reg⁻[b, t] <= pb_in_max[b] - pb_in[b, t] + pb_out[b, t]
            )
        end
    end

    return
end