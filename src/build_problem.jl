function apply_cc_constraints!(problem)
    optimization_container = PSI.get_optimization_container(problem)
    restrictions = problem.ext["cc_restrictions"]
    commitment_variables = PSI.get_variable(optimization_container, :On__ThermalMultiStart)
    time_steps = PSI.model_time_steps(optimization_container)
    constraint = PSI.add_cons_container!(
        optimization_container,
        :CC_constraint,
        collect(keys(restrictions)),
        time_steps,
    )
    jump_model = PSI.get_jump_model(optimization_container)
    for t in time_steps, (k, v) in restrictions
        constraint[k, t] =
            JuMP.@constraint(jump_model, sum(commitment_variables[i, t] for i in v) <= 1)
    end
    return
end

function apply_must_run_constraints!(problem)
    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    must_run_gens =
        [g for g in PSY.get_components(ThermalMultiStart, system, x -> PSY.get_must_run(x))]
    commitment_variables = PSI.get_variable(optimization_container, :On__ThermalMultiStart)
    for t in time_steps, g in get_name.(must_run_gens)
        JuMP.fix(commitment_variables[g, t], 1.0)
    end
end

function PSI.problem_build!(problem::PSI.OperationsProblem{MultiStartUnitCommitmentCC})
    PSI.build_impl!(
        PSI.get_optimization_container(problem),
        PSI.get_template(problem),
        PSI.get_system(problem),
    )
    apply_cc_constraints!(problem)
    apply_must_run_constraints!(problem)
end
