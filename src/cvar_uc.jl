struct CVaRUnitCommitmentCC <: PSI.PowerSimulationsOperationsProblem end

function PSI.problem_build!(problem::PSI.OperationsProblem{CVaRUnitCommitmentCC})
    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    time_steps = PSI.model_time_steps(optimization_container)
    jump_model = PSI.get_jump_model(optimization_container)

    # Here is where you build the JuMP model

end
