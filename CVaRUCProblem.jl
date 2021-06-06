# To run: julia --project CVaRUCProblem.jl true true true true true true 1000 0.25 0.80 Power

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

# May 1st: under low tail in afternoon (Day 121)
# June 20th: Low tail in morning (Day 171)
initial_time = "2018-06-20T00:00:00"
use_storage = isempty(ARGS) ? true : parse(Bool, ARGS[1])
use_storage_reserves = isempty(ARGS) ? true : parse(Bool, ARGS[2])
use_reg = isempty(ARGS) ? true : parse(Bool, ARGS[3])
use_spin = isempty(ARGS) ? true : parse(Bool, ARGS[4])
use_must_run = isempty(ARGS) ? true : parse(Bool, ARGS[5])
use_nuclear = isempty(ARGS) ? true : parse(Bool, ARGS[6])
C_RR = isempty(ARGS) ? 5000 : parse(Float64, ARGS[7]) # Penalty cost of recourse reserve
L_SUPP = isempty(ARGS) ? 1 / 4 : parse(Float64, ARGS[8]) # 15 min response time, to start
α = isempty(ARGS) ? 0.8 : parse(Float64, ARGS[9]) # Risk tolerance level
formulation = isempty(ARGS) ? "Power" : ARGS[10]
scenarios = 31

if formulation == "Power"
    custom_problem = CVaRPowerUnitCommitmentCC
elseif formulation == "Reserve"
    custom_problem = CVaRReserveUnitCommitmentCC
else
    throw(ArgumentError("Unrecognized formulation type. Given " * formulation))
end

optional_title =
    # (use_storage ? " stor" : "") *
    # (use_storage_reserves ? " storres" : "") *
    # (use_reg ? " reg" : "") *
    # (use_spin ? " spin" : "") *
    # (!use_must_run ? " no must run" : "") *
    " C_RR " * string(C_RR) * " alpha " * string(α)

output_path = "./results/" * string(scenarios) * " scenarios/CVaR/" * formulation * "/" * split(initial_time, "T")[1] * optional_title * "/"
if !isdir(output_path)
    mkpath(output_path)
end

## Jose
# system_file_path = "/Users/jdlara/cache/blue_texas/"
## Kate
system_file_path = "data/"

system_da = System(
    joinpath(system_file_path, "DA_sys_" * string(scenarios) * "_scenarios.json");
    time_series_read_only = true,
)
# system_ha = System("data/HA_sys.json"; time_series_read_only = true)
# system_ed = System("data/RT_sys.json"; time_series_read_only = true)

# Jose's tune-ups for the HA UC
for system in [system_da] # [system_da, system_ha, system_ed]
    apply_manual_data_updates!(system, use_nuclear)
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
    custom_problem,
    template_dauc,
    system_da,
    optimizer = solver,
    initial_time = DateTime(initial_time),
    optimizer_log_print = true,
    balance_slack_variables = false,
)
UC.ext["cc_restrictions"] =
    JSON.parsefile(joinpath(system_file_path, "cc_restrictions.json"))
UC.ext["use_storage"] = use_storage
UC.ext["use_storage_reserves"] = use_storage_reserves
UC.ext["use_reg"] = use_reg
UC.ext["use_spin"] = use_spin
UC.ext["use_must_run"] = use_must_run
UC.ext["C_RR"] = C_RR
UC.ext["L_SUPP"] = L_SUPP
UC.ext["α"] = α
UC.ext["C_res_penalty"] = 5000
UC.ext["C_ener_penalty"] = 100000
UC.ext["L_REG"] = 1 / 12 # 5 min
UC.ext["L_SPIN"] = 1 / 6 # 10 min

# Build and solve the standalone problem
build!(UC; output_dir = output_path, serialize = false) # use serialize=true to get OptimizationModel.json to debug
(status, solvetime) = @timed solve!(UC)

if status.value == 0
    write_to_CSV(UC, output_path; time=solvetime)

    for scenario in 1:scenarios
        plot_fuel(
            UC;
            scenario = scenario,
            save_dir = output_path,
        )
    end
end
