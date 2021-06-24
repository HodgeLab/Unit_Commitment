include("src/Unit_commitment.jl")
plotlyjs()

## Local
# using Xpress
# solver = optimizer_with_attributes(Xpress.Optimizer, "MIPRELSTOP" => 0.01) # MIPRELSTOP was  0.0001
## Eagle
using Gurobi
solver = optimizer_with_attributes(Gurobi.Optimizer, "MIPGap" => 0.01)

############################## First Stage Problem Definition ##############################
formulation = isempty(ARGS) ? "D" : ARGS[1]
initial_time = isempty(ARGS) ? "2018-03-15T00:00:00" : ARGS[2]
use_storage = isempty(ARGS) ? true : parse(Bool, ARGS[3])
use_storage_reserves = isempty(ARGS) ? true : parse(Bool, ARGS[4])
use_solar_reg = isempty(ARGS) ? true : parse(Bool, ARGS[5])
use_solar_spin = isempty(ARGS) ? true : parse(Bool, ARGS[6])
use_must_run = isempty(ARGS) ? true : parse(Bool, ARGS[7])
use_nuclear = isempty(ARGS) ? true : parse(Bool, ARGS[8])
C_RR = isempty(ARGS) ? 5000 : parse(Float64, ARGS[9]) # Penalty cost of recourse reserve
α = isempty(ARGS) ? 0.8 : parse(Float64, ARGS[10]) # Risk tolerance level
supp_type = isempty(ARGS) ? "generic" : ARGS[11]
scenarios = 31

scenario_plot_dict = Dict{String, Vector{Int64}}(
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

if initial_time in keys(scenario_plot_dict)
    plot_scenarios = scenario_plot_dict[initial_time]
else
    plot_scenarios = 1:scenarios
end

if formulation == "D"
    formulation_dir = "Deterministic"
    custom_problem = BasecaseUnitCommitmentCC
elseif formulation == "C"
    formulation_dir = "CVAR"
    custom_problem = CVaRReserveUnitCommitmentCC
    if !(supp_type in ["generic", "nonspin"])
        throw(ArgumentError("Supp reserves must be generic or nonspin"))
    end
elseif formulation == "S"
    formulation_dir = "Stochastic"
    custom_problem = StochasticUnitCommitmentCC
else
    throw(ArgumentError("Formulation key unrecognized"))
end

optional_title =
    (use_storage ? " stor" : "") *
    (use_storage_reserves ? " storres" : "") *
    (use_solar_reg ? " solreg" : "") *
    (use_solar_spin ? " solspin" : "") *
    (formulation == "C" ? " C_RR " * string(C_RR) * " alpha " * string(α) : "") *
    (formulation == "C" ? " " * supp_type : "")

output_path =
    "./results/" *
    string(scenarios) *
    " scenarios/" *
    formulation_dir *
    "/" *
    split(initial_time, "T")[1] *
    optional_title *
    "/"
if !isdir(output_path)
    mkpath(output_path)
end

## Jose
# system_file_path = "/Users/jdlara/Dropbox/texas_data"
# simulation_folder = mktempdir()
## Kate
system_file_path = "data/"
simulation_folder = output_path

system_da = System(
    joinpath(system_file_path, "DA_sys_" * string(scenarios) * "_scenarios.json");
    time_series_read_only = true,
)

initial_cond_file = joinpath("data/", "initial_on_" * split(initial_time, "T")[1] * ".csv")
if !isfile(initial_cond_file)
    initial_cond_file = joinpath("data/", "initial_on.csv")
end

apply_manual_data_updates!(system_da, use_nuclear, initial_cond_file)

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
UC.ext["use_wind_reserves"] = false
UC.ext["use_solar_reg"] = use_solar_reg
UC.ext["use_solar_spin"] = use_solar_spin
UC.ext["use_reg"] = true
UC.ext["use_spin"] = true
UC.ext["use_must_run"] = use_must_run
UC.ext["C_RR"] = C_RR * get_base_power(system_da)
UC.ext["α"] = α
UC.ext["C_res_penalty"] = 5000 * get_base_power(system_da)
UC.ext["C_ener_penalty"] = 9000 * get_base_power(system_da)
UC.ext["L_REG"] = 1 / 12 # 5 min
UC.ext["L_SPIN"] = 1 / 6 # 10 min
UC.ext["L_SUPP"] = 1 / 6 # 10 min
UC.ext["load_scale"] = 1
UC.ext["solar_scale"] = 1
UC.ext["storage_scale"] = 1
UC.ext["solar_reg_prop"] = 1
UC.ext["solar_spin_prop"] = 1
UC.ext["wind_reg_prop"] = 1
UC.ext["wind_spin_prop"] = 1
UC.ext["renewable_reg_prop"] = 1
UC.ext["renewable_spin_prop"] = 1
UC.ext["supp_type"] = supp_type
UC.ext["allowable_reserve_prop"] = 0.2 # Can use up to 20% total for all reserves

#################################### Stage 2 problem Definition, HAUC ########################
system_ha = System(
    joinpath(system_file_path, "HA_sys_UC_experiment.json");
    time_series_read_only = true,
)

apply_manual_data_updates!(system_ha, use_nuclear, initial_cond_file)

template_hauc = OperationsProblemTemplate(CopperPlatePowerModel)
set_device_model!(template_hauc, RenewableDispatch, RenewableFullDispatch)
set_device_model!(template_hauc, PowerLoad, StaticPowerLoad)
# Use FixedOutput instead of HydroDispatchRunOfRiver to get consistent results because model might decide to curtail wind vs. hydro (same cost)
set_device_model!(template_hauc, HydroDispatch, FixedOutput)
set_service_model!(template_hauc, ServiceModel(VariableReserve{ReserveUp}, RangeReserve))
set_service_model!(template_hauc, ServiceModel(VariableReserve{ReserveDown}, RangeReserve))
set_device_model!(template_hauc, GenericBattery, BookKeepingwReservation)
set_device_model!(template_hauc, ThermalMultiStart, ThermalMultiStartUnitCommitment)

#################################### Solve Stage 1 Problem ################################

# Build and solve the Stage 1
build!(UC; output_dir = output_path, serialize = false) # use serialize=true to get OptimizationModel.json to debug
(status, solvetime) = @timed solve!(UC)

if status.value != 0
    throw(ErrorException("DAUC failed"))
end

#################################### Solve Stage 2 Problem ################################

for h in 1:24

    hauc_initial_time = DateTime(initial_time) + Hour(h - 1)

    HAUC = OperationsProblem(
        HourAheadUnitCommitmentCC,
        template_hauc,
        system_ha,
        optimizer = solver,
        initial_time = hauc_initial_time,
        optimizer_log_print = false,
        balance_slack_variables = true,
        system_to_file = false,
    )
    for (k, v) in UC.ext
        HAUC.ext[k] = v
    end
    HAUC.ext["UC_obj_dict"] = PSI.get_jump_model(PSI.get_optimization_container(UC)).obj_dict
    HAUC.ext["step"] = h
    
    build!(HAUC; output_dir = output_path, serialize = false) # use serialize=true to get OptimizationModel.json to debug
    (status, solvetime) = @timed solve!(HAUC)

    if status.value != 0
        throw(ErrorException("HAUC failed at step " * string(h)))
    end

    write_to_CSV(HAUC, output_path; append = (h == 1 ? false : true))

end

#################################### Export  ################################

# Record stage 1 outputs
write_to_CSV(UC, output_path)
write_summary_stats(UC, output_path, solvetime)

for scenario in (formulation == "D" ? [nothing] : plot_scenarios)
    plot_fuel(UC; scenario = scenario, save_dir = output_path)

    plot_reserve(UC, "SPIN"; save_dir = output_path, scenario = scenario)

    plot_reserve(UC, "REG_UP"; save_dir = output_path, scenario = scenario)

    plot_reserve(UC, "REG_DN"; save_dir = output_path, scenario = scenario)
end

# Stage 2 outputs
    
# Plot and Save
