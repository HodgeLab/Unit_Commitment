include("src/Unit_commitment.jl")
using PowerSimulations
using PowerSystems
using Dates
# using PowerGraphics

## Local
# using Xpress
# solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.1) # MIPRELSTOP was  0.0001
## Eagle
using Gurobi
solver = optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.1)

output_path = "./results/CVaR"
## Jose
# system_file_path = "/Users/jdlara/cache/blue_texas/"
## Kate
system_file_path = "data/"

system_da = System(joinpath(system_file_path, "DA_sys.json"); time_series_read_only = true)
# system_ha = System("data/HA_sys.json"; time_series_read_only = true)
# system_ed = System("data/RT_sys.json"; time_series_read_only = true)

# Jose's tune-ups for the HA UC
for system in [system_da] # [system_da, system_ha, system_ed]
    for g in get_components(
        ThermalMultiStart,
        system,
        x -> get_prime_mover(x) == ThermalFuels.NATURAL_GAS,
    )
        lims = get_ramp_limits(g)
        set_time_limits!(g, (up = 1.1 * lims.up, down = 1.1 * lims.down))
        set_ramp_limits!(g, (up = 1.25 * lims.up, down = 1.25 * lims.down))
    end

    for g in get_components(
        ThermalMultiStart,
        system,
        x -> get_prime_mover(x) == ThermalFuels.COAL,
    )
        lims = get_ramp_limits(g)
        set_ramp_limits!(g, (up = 1.25 * lims.up, down = 1.25 * lims.down))
    end
    for g in get_components(
        RenewableDispatch,
        system,
        x -> get_prime_mover(x) != PrimeMovers.PVe,
    )
        set_available!(g, true)
    end

    for g in get_components(HydroGen, system)
        set_available!(g, true)
    end

    s = get_component(VariableReserve{ReserveUp}, system, "REG_UP")
    req = get_requirement(s)
    set_requirement!(s, req * 1.5)

    s = get_component(VariableReserve{ReserveDown}, system, "REG_DN")
    req = get_requirement(s)
    set_requirement!(s, req * 1.5)
end

# Set all CC's to start off
for g in get_components(
    ThermalMultiStart,
    system_da,
    x -> get_prime_mover(x) in [PrimeMovers.CT, PrimeMovers.CC],
)
    set_status!(g, false)
    set_active_power!(g, 0.0)
end

template_dauc = OperationsProblemTemplate(CopperPlatePowerModel)
# template_hauc = OperationsProblemTemplate(CopperPlatePowerModel)
# template_ed = OperationsProblemTemplate(CopperPlatePowerModel)
for template in [template_dauc] # [template_dauc, template_hauc, template_ed]
    set_device_model!(template, RenewableDispatch, RenewableFullDispatch)
    set_device_model!(template, PowerLoad, StaticPowerLoad)
    # Use FixedOutput instead of HydroDispatchRunOfRiver to get consistent results because model might decide to curtail wind vs. hydro (same cost)
    set_device_model!(template, HydroDispatch, FixedOutput)
    set_service_model!(template, ServiceModel(VariableReserve{ReserveUp}, RangeReserve))
    set_service_model!(template, ServiceModel(VariableReserve{ReserveDown}, RangeReserve))
    set_device_model!(template, GenericBattery, BookKeepingwReservation)
end

set_device_model!(template_dauc, ThermalMultiStart, ThermalMultiStartUnitCommitment)
# ignore HA for now
# set_device_model!(template_ed, ThermalMultiStart, ThermalRampLimited)

UC = OperationsProblem(
    CVaRUnitCommitmentCC,
    template_dauc,
    system_da,
    optimizer = solver,
    initial_time = DateTime("2018-04-01T00:00:00"),
    optimizer_log_print = true,
    balance_slack_variables = true,
)
UC.ext["cc_restrictions"] =
    JSON.parsefile(joinpath(system_file_path, "cc_restrictions.json"))

# Build and solve the standalone problem
build!(UC; output_dir = output_path, serialize = false) # Can add balance_slack_variables (load shedding and curtailment), use serialize=true to get OptimizationModel.json to debug
solve!(UC)
problem_results = ProblemResults(UC)
write_to_CSV(problem_results, output_path)
