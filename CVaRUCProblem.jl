# To run: julia --project CVaRUCProblem.jl true true true true

include("src/Unit_commitment.jl")
using PowerSimulations
using PowerSystems
using Dates
using CSV
using DataFrames
using PowerGraphics
plotlyjs()

## Local
using Xpress
solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.1) # MIPRELSTOP was  0.0001
## Eagle
# using Gurobi
# solver = optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.1)

initial_time = "2018-04-20T00:00:00"
use_storage = isempty(ARGS) ? false : (ARGS[1] == "true")
use_storage_reserves = isempty(ARGS) ? false : (ARGS[2] == "true")
use_reg = isempty(ARGS) ? false : (ARGS[3] == "true")
use_spin = isempty(ARGS) ? false : (ARGS[4] == "true")
optional_title = (use_storage ? " stor" : "") *
                (use_storage_reserves ? " storres" : "") *
                (use_reg ? " reg" : "") *
                (use_spin ? " spin" : "")

output_path = "./results/CVaR/" * split(initial_time, "T")[1] * optional_title * "/"
if !isdir(output_path)
    mkpath(output_path)
end

## Jose
# system_file_path = "/Users/jdlara/cache/blue_texas/"
## Kate
system_file_path = "data/"

system_da = System(joinpath(system_file_path, "DA_sys_84_scenarios.json"); time_series_read_only = true)
# system_ha = System("data/HA_sys.json"; time_series_read_only = true)
# system_ed = System("data/RT_sys.json"; time_series_read_only = true)

# Jose's tune-ups for the HA UC
for system in [system_da] # [system_da, system_ha, system_ed]
    appply_manual_data_updates!(system)
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
    initial_time = DateTime(initial_time),
    optimizer_log_print = true,
    balance_slack_variables = true,
)
UC.ext["cc_restrictions"] =
    JSON.parsefile(joinpath(system_file_path, "cc_restrictions.json"))
UC.ext["use_storage"] = use_storage
UC.ext["use_storage_reserves"] = use_storage_reserves
UC.ext["use_reg"] = use_reg
UC.ext["use_spin"] = use_spin

# Build and solve the standalone problem
build!(UC; output_dir = output_path, serialize = false) # use serialize=true to get OptimizationModel.json to debug
solve!(UC)

plot_fuel(UC; 
    case_initial_time = DateTime(initial_time), 
    storage = use_storage,
    scenario = 1,
    save=output_path
)

plot_fuel(UC; 
    case_initial_time = DateTime(initial_time), 
    storage = use_storage,
    scenario = 80,
    save=output_path
)

write_to_CSV(UC, output_path)
