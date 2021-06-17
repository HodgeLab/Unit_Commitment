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
use_solar_reg = isempty(ARGS) ? true : parse(Bool, ARGS[3])
use_solar_spin = isempty(ARGS) ? true : parse(Bool, ARGS[4])
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

template_dauc = OperationsProblemTemplate(CopperPlatePowerModel)

set_device_model!(template_dauc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_dauc, PowerLoad, StaticPowerLoad)
# Use FixedOutput instead of HydroDispatchRunOfRiver to get consistent results because model might decide to curtail wind vs. hydro (same cost)
set_device_model!(template_dauc, HydroDispatch, FixedOutput)
set_service_model!(template_dauc, ServiceModel(VariableReserve{ReserveUp}, RangeReserve))
set_service_model!(template_dauc, ServiceModel(VariableReserve{ReserveDown}, RangeReserve))
set_device_model!(template_dauc, GenericBattery, BookKeepingwReservation)

set_device_model!(template_dauc, ThermalMultiStart, ThermalMultiStartUnitCommitment)
# ignore HA for now
# set_device_model!(template_ed, ThermalMultiStart, ThermalRampLimited)

optional_title =
    (use_storage ? " stor" : "") *
    (use_storage_reserves ? " storres" : "") *
    (use_solar_reg ? " solreg" : "") *
    (use_solar_spin ? " solspin" : "")

scenario_plot_dict = Dict{String,Vector{Int64}}(
    "2018-03-15T00:00:00" => [30, 29],
    "2018-03-27T00:00:00" => [31, 13],
    "2018-04-15T00:00:00" => [27, 30],
    "2018-05-17T00:00:00" => [31, 19],
    "2018-07-22T00:00:00" => [1, 5],
    "2018-07-24T00:00:00" => [30, 13],
    "2018-08-15T00:00:00" => [3, 29],
    "2018-09-21T00:00:00" => [28, 23],
    "2018-09-24T00:00:00" => [28, 14],
    "2018-10-08T00:00:00" => [31, 23],
    "2018-11-09T00:00:00" => [3, 18],
    "2018-12-07T00:00:00" => [3, 24],
    "2018-12-26T00:00:00" => [25, 13],
)

# days_per_month = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
# for month in 1:12
#     for d in 1:days_per_month[month]
#         initial_time = "2018-" * (month < 10 ? "0" * string(month) : string(month)) *
#             "-" * string(d) * "T00:00:00"

for initial_time in keys(scenario_plot_dict)

        output_path = "./results/" * string(scenarios) * " scenarios/Deterministic/" * split(initial_time, "T")[1] * optional_title * "/"

        # if !isdir(output_path)
        #     mkpath(output_path)
        # end

        apply_manual_data_updates!(system_da, use_nuclear, joinpath(system_file_path, "initial_on_" * split(initial_time, "T")[1] * ".csv"))

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
        UC.ext["storage_reserve_names"] = ["EXPOSE_STORAGE"]
        UC.ext["use_wind_reserves"] = false
        UC.ext["use_solar_reg"] = use_solar_reg
        UC.ext["use_solar_spin"] = use_solar_spin
        UC.ext["use_reg"] = true
        UC.ext["use_spin"] = true
        UC.ext["use_must_run"] = use_must_run
        UC.ext["C_res_penalty"] = 5000 * get_base_power(system_da)
        UC.ext["C_ener_penalty"] = 9000 * get_base_power(system_da)
        UC.ext["L_REG"] = 1 / 12 # 5 min
        UC.ext["L_SPIN"] = 1 / 6 # 10 min
        UC.ext["load_scale"] = 1
        UC.ext["solar_scale"] = 1
        UC.ext["storage_scale"] = 1
        UC.ext["solar_reg_prop"] = 1
        UC.ext["wind_reg_prop"] = 1
        UC.ext["solar_spin_prop"] = 1
        UC.ext["wind_spin_prop"] = 1
        UC.ext["renewable_reg_prop"] = 1
        UC.ext["renewable_spin_prop"] = 1
        UC.ext["supp_type"] = "generic"
        UC.ext["allowable_reserve_prop"] = 0.2 # Can use up to 20% total for all reserves

        # Build and solve the standalone problem
        build!(UC; output_dir = output_path, serialize = false) # Can add balance_slack_variables (load shedding and curtailment), use serialize=true to get OptimizationModel.json to debug
        (status, solvetime) = @timed solve!(UC)

        if status.value == 0
            hour = 3
            print(initial_time * " solved")
            save_as_initial_condition(UC,
                joinpath(system_file_path, "initial_on_" * split(initial_time, "T")[1] * ".csv"),
                hour
            )

            # write_to_CSV(
            #     UC,
            #     system_file_path,
            #     output_path;
            #     time=solvetime
            # )

            # plot_fuel(
            #     UC;
            #     save_dir = output_path,
            #     scenario = nothing
            # )

            # plot_reserve(
            #     UC,
            #     "REG_UP";
            #     save_dir = output_path,
            #     scenario = nothing
            # )

            # plot_reserve(
            #     UC,
            #     "REG_DN";
            #     save_dir = output_path,
            #     scenario = nothing
            # )

            # plot_reserve(
            #     UC,
            #     "SPIN";
            #     save_dir = output_path,
            #     scenario = nothing
            # )

        end
#     end
end