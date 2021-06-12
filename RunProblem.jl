# To run: julia --project RunProblem.jl D 2018-05-17T00:00:00 true true true true true true 1000 0.80
# D for deterministic, S for stochastic, C for CVaR

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

# April 15: Totally clear (Day 105)
# May 14th: under low tail midday, in low tail in afternoon (Day 134)
# June 13th: Low tail in morning (Day 164)
# May 17th: Day 137 (low in afternoon, but still mostly in range)
formulation = isempty(ARGS) ? "D" : ARGS[1]
initial_time = isempty(ARGS) ? "2018-05-17T00:00:00" : ARGS[2]
use_storage = isempty(ARGS) ? true : parse(Bool, ARGS[3])
use_storage_reserves = isempty(ARGS) ? true : parse(Bool, ARGS[4])
use_solar_reserves = isempty(ARGS) ? true : parse(Bool, ARGS[5])
use_spin = isempty(ARGS) ? true : parse(Bool, ARGS[6])
use_must_run = isempty(ARGS) ? true : parse(Bool, ARGS[7])
use_nuclear = isempty(ARGS) ? true : parse(Bool, ARGS[8])
C_RR = isempty(ARGS) ? 5000 : parse(Float64, ARGS[9]) # Penalty cost of recourse reserve
α = isempty(ARGS) ? 0.8 : parse(Float64, ARGS[10]) # Risk tolerance level
scenarios = 31

if formulation == "D"
    formulation_dir = "Deterministic"
    custom_problem = BasecaseUnitCommitmentCC
elseif formulation == "C"
    formulation_dir = "CVAR"
    custom_problem = CVaRReserveUnitCommitmentCC
elseif formulation == "S"
    formulation_dir = "Stochastic"
    custom_problem = StochasticUnitCommitmentCC
else throw(ArgumentError("Formulation key unrecognized"))
end

optional_title =
    (use_storage ? " stor" : "") *
    (use_storage_reserves ? " storres" : "") *
    (use_solar_reserves ? " solres" : "") *
    (formulation == "C" ? " C_RR " * string(C_RR) * " alpha " * string(α) : "")

output_path = "./results/" * string(scenarios) * " scenarios/" * formulation_dir * 
    "/" * split(initial_time, "T")[1] * optional_title * "/"
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

apply_manual_data_updates!(system_da, use_nuclear, system_file_path)

template_dauc = OperationsProblemTemplate(CopperPlatePowerModel)
set_device_model!(template_dauc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_dauc, PowerLoad, StaticPowerLoad)
# Use FixedOutput instead of HydroDispatchRunOfRiver to get consistent results because model might decide to curtail wind vs. hydro (same cost)
set_device_model!(template_dauc, HydroDispatch, FixedOutput)
set_service_model!(template_dauc, ServiceModel(VariableReserve{ReserveUp}, RangeReserve))
set_service_model!(template_dauc, ServiceModel(VariableReserve{ReserveDown}, RangeReserve))
set_device_model!(template_dauc, GenericBattery, BookKeepingwReservation)

set_device_model!(template_dauc, ThermalMultiStart, ThermalMultiStartUnitCommitment)

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
UC.ext["storage_reserve_names"] = ["EXPOSE_STORAGE"]
UC.ext["use_solar_reserves"] = use_solar_reserves
UC.ext["use_reg"] = true
UC.ext["use_spin"] = use_spin
UC.ext["use_must_run"] = use_must_run
UC.ext["C_RR"] = C_RR
UC.ext["α"] = α
UC.ext["C_res_penalty"] = 5000
UC.ext["C_ener_penalty"] = 100000
UC.ext["L_REG"] = 1 / 12 # 5 min
UC.ext["L_SPIN"] = 1 / 6 # 10 min

# Build and solve the standalone problem
build!(UC; output_dir = output_path, serialize = false) # use serialize=true to get OptimizationModel.json to debug
(status, solvetime) = @timed solve!(UC)

if status.value == 0
    write_to_CSV(
        UC,
        system_file_path,
        output_path;
        time=solvetime
    )

    for scenario in 1:(formulation == "D" ? 1 : scenarios)
        plot_fuel(
            UC;
            scenario = (formulation == "D" ? nothing : scenario),
            save_dir = output_path,
        )

        plot_reserve(
                UC,
                "REG_UP";
                use_solar_reserves = use_solar_reserves,
                save_dir = output_path,
                scenario = (formulation == "D" ? nothing : scenario)
            )

            plot_reserve(
                UC,
                "REG_DN";
                use_solar_reserves = use_solar_reserves,
                save_dir = output_path,
                scenario = (formulation == "D" ? nothing : scenario)
            )
    end

    plot_reserve(
        UC,
        "SPIN";
        use_solar_reserves = use_solar_reserves,
        save_dir = output_path,
        scenario = nothing
    )
end
