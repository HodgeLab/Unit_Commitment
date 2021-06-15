# To run: julia --project DeterministicProblem.jl true true true true true true

include("src/Unit_commitment.jl")
using PowerSimulations
using PowerSystems
using Dates
using CSV
using HDF5
using DataFrames
using PowerGraphics
plotlyjs()

## Local
using Xpress
solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.1) # MIPRELSTOP was  0.0001
## Eagle
# using Gurobi
# solver = optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.1)


use_storage = isempty(ARGS) ? true : parse(Bool, ARGS[1])
use_storage_reserves = isempty(ARGS) ? true : parse(Bool, ARGS[2])
use_solar_reserves = isempty(ARGS) ? true : parse(Bool, ARGS[3])
use_spin = isempty(ARGS) ? true : parse(Bool, ARGS[4])
use_must_run = isempty(ARGS) ? true : parse(Bool, ARGS[5])
use_nuclear = isempty(ARGS) ? true : parse(Bool, ARGS[6])
scenarios = 31

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

optional_title =
    (use_storage ? " stor" : "") *
    (use_storage_reserves ? " storres" : "") *
    (use_solar_reserves ? " solres" : "")

days_per_month = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
for month in 1:12
    for d in 1:days_per_month[month]
        initial_time = "2018-" * (month < 10 ? "0" * string(month) : string(month)) *
            "-" * string(d) * "T00:00:00"

        output_path = "./results/" * string(scenarios) * " scenarios/Deterministic/" * split(initial_time, "T")[1] * optional_title * "/"

        if !isdir(output_path)
            mkpath(output_path)
        end

        # Jose's tune-ups for the HA UC
        for system in [system_da] # [system_da, system_ha, system_ed]
            apply_manual_data_updates!(system, use_nuclear)
        end

        UC = OperationsProblem(
            BasecaseUnitCommitmentCC,
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
        UC.ext["use_solar_reserves"] = use_solar_reserves
        UC.ext["use_reg"] = true
        UC.ext["use_spin"] = use_spin
        UC.ext["use_must_run"] = use_must_run
        UC.ext["C_res_penalty"] = 5000
        UC.ext["C_ener_penalty"] = 100000
        UC.ext["L_REG"] = 1 / 12 # 5 min
        UC.ext["L_SPIN"] = 1 / 6 # 10 min
        UC.ext["load_scale"] = 1
        UC.ext["solar_scale"] = 1
        UC.ext["storage_scale"] = 15
        UC.ext["solar_reg_prop"] = 1
        UC.ext["wind_reg_prop"] = 1
        UC.ext["solar_spin_prop"] = 1
        UC.ext["wind_spin_prop"] = 1
        UC.ext["renewable_reg_prop"] = 1
        UC.ext["renewable_spin_prop"] = 1

        # Build and solve the standalone problem
        build!(UC; output_dir = output_path, serialize = false) # Can add balance_slack_variables (load shedding and curtailment), use serialize=true to get OptimizationModel.json to debug
        (status, solvetime) = @timed solve!(UC)

        if status.value == 0
            write_to_CSV(
                UC,
                system_file_path,
                output_path;
                time=solvetime
            )

            plot_fuel(
                UC;
                save_dir = output_path,
                scenario = nothing
            )

            plot_reserve(
                UC,
                "REG_UP";
                save_dir = output_path,
                scenario = nothing
            )

            plot_reserve(
                UC,
                "REG_DN";
                save_dir = output_path,
                scenario = nothing
            )

            plot_reserve(
                UC,
                "SPIN";
                save_dir = output_path,
                scenario = nothing
            )

        end
    end
end