function PG.plot_fuel(
    problem::PSI.OperationsProblem{T};
    kwargs...,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
        StochasticUnitCommitmentCC,
    },
}
    title = get(kwargs, :title, "DAUC Fuel")
    save_dir = get(kwargs, :save_dir, nothing)
    scenario = kwargs[:scenario]
    time_steps = get(kwargs, :time_steps, nothing)
    charge_line = get(kwargs, :charge_line, false)

    p = PG._empty_plot()
    backend = Plots.backend()

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    case_initial_time = PSI.get_initial_time(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    if isnothing(time_steps)
        time_steps = PSI.model_time_steps(optimization_container)
    end
    use_slack =
        PSI.get_balance_slack_variables(problem.internal.optimization_container.settings)
    storage = problem.ext["use_storage"]

    total_load =
        get_area_total_time_series(problem, PowerLoad)[time_steps] .*
        problem.ext["load_scale"]
    total_hydro = get_area_total_time_series(problem, HydroGen)[time_steps]
    total_wind = get_area_total_time_series(
        problem,
        RenewableGen;
        filter = x -> get_prime_mover(x) != PrimeMovers.PVe,
    )[time_steps]

    solar_forecast = _get_solar_forecast(problem, time_steps; kwargs...)

    gen = get_generation_data(
        problem,
        total_wind,
        total_hydro,
        scenario,
        solar_forecast,
        time_steps;
        kwargs...,
    )
    cat = make_fuel_dictionary(system)

    fuel = my_categorize_data(gen.data, cat, use_slack, storage)

    # Hack to make nuclear on the bottom
    cat_names = intersect(PG.CATEGORY_DEFAULT, keys(fuel))
    cat_names = [cat_names[2], cat_names[1], cat_names[3:end]...]
    fuel_agg = PG.combine_categories(fuel; names = cat_names)

    y_label = get(kwargs, :y_label, "Generation (GW)")

    seriescolor = get(kwargs, :seriescolor, PG.match_fuel_colors(fuel_agg, backend))
    DataFrames.rename!(fuel_agg, Dict("Natural gas" => "NG-CT (solo)"))
    DataFrames.rename!(fuel_agg, Dict("NG-CT" => "NG-CT (in CC train)"))
    p = plot_dataframe(
        fuel_agg,
        gen.time;
        seriescolor = seriescolor,
        y_label = y_label,
        title = nothing,
        stack = true,
        set_display = false,
        stair = true,
        kwargs...,
    )

    kwargs = Dict{Symbol, Any}((k, v) for (k, v) in kwargs if k ∉ [:nofill, :seriescolor])
    kwargs[:linestyle] = get(kwargs, :linestyle, :dash)
    kwargs[:linewidth] = get(kwargs, :linewidth, 3)

    # Add load line
    load_agg =
        DataFrames.DataFrame(Dict(:Load => total_load)) .* get_base_power(system) ./ 1000
    p = plot_dataframe(
        p,
        load_agg,
        gen.time;
        seriescolor = ["black"],
        y_label = y_label,
        title = nothing,
        stack = true,
        nofill = true,
        set_display = false,
        stair = true,
        kwargs...,
    )

    # Add charging load line
    if charge_line && storage
        load_and_charging = load_agg .+ sum.(eachrow(gen.data[:pb_in]))
        DataFrames.rename!(load_and_charging, Symbol.(["Load + charging"]))
        col = PG.match_fuel_colors(fuel_agg[!, ["Storage"]], backend)
        p = plot_dataframe(
            p,
            load_and_charging,
            gen.time;
            seriescolor = [col],
            y_label = y_label,
            title = nothing,
            stack = true,
            nofill = true,
            set_display = false,
            stair = true,
            kwargs...,
        )
    end

    # Overwrite x axis label
    layout_kwargs =
        Dict{Symbol, Any}(:xaxis => Plots.PlotlyJS.attr(; title = "Time of Day"))
    Plots.PlotlyJS.relayout!(p, Plots.PlotlyJS.Layout(; layout_kwargs...))

    if !isnothing(save_dir)
        title = replace(title, " " => "_")
        for format in ("png", "pdf")
            fname = _get_save_path(problem, title, format, save_dir; kwargs...)
            # Overwrite existing plots
            if isfile(fname)
                rm(fname)
            end
            PG.save_plot(p, fname, backend; kwargs...)
        end
    end
    return p
end

function PG.get_generation_data(
    problem::PSI.OperationsProblem{T},
    total_wind,
    total_hydro,
    scenario,
    solar_forecast,
    time_steps;
    kwargs...,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
        StochasticUnitCommitmentCC,
    },
}
    storage = problem.ext["use_storage"]
    use_slack =
        PSI.get_balance_slack_variables(problem.internal.optimization_container.settings)

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)

    var_names = Vector{Symbol}()
    if storage
        push!(var_names, :pb_in, :pb_out)
    end
    if use_slack
        push!(var_names, :slack_energy⁺)
    end
    variables = Dict{Symbol, DataFrames.DataFrame}()
    for v in var_names
        variables[v] =
            PSI.axis_array_to_dataframe(jump_model.obj_dict[v], [v])[time_steps, :]
    end

    variables[:pg] = get_thermal_generator_power_dataframe(problem, time_steps, scenario)

    # Select single solar scenario
    variables[:pS] = _get_solar_realization(problem, time_steps; kwargs...)
    variables[:pW] =
        PSI.axis_array_to_dataframe(jump_model.obj_dict[:pW], [:pW])[time_steps, :]
    variables[:pH] = DataFrames.DataFrame(Dict(:pH => total_hydro))

    variables[:curt] = DataFrames.DataFrame(
        :curt =>
            total_wind - variables[:pW][!, 1] + solar_forecast - variables[:pS][!, 1],
    )

    # Scale from 100 MW to GW
    for v in keys(variables)
        variables[v] .*= get_base_power(system) ./ 1000
    end

    timestamps = get_timestamps(problem)[time_steps]
    return PG.PGData(variables, timestamps)
end

function _scenario_in_3D_array_to_dataframe(
    input_array::JuMP.Containers.DenseAxisArray{},
    scenario,
    time_steps,
)
    result = Array{Float64, 2}(undef, length(time_steps), length(input_array.axes[1]))
    names = Array{Symbol, 1}(undef, length(input_array.axes[1]))
    for t in time_steps, (ix, name) in enumerate(input_array.axes[1])
        result[t, ix] = PSI._jump_value(input_array[name, scenario, t])
        names[ix] = Symbol(name)
    end

    return DataFrames.DataFrame(result, names)
end

function my_categorize_data(
    data::Dict{Symbol, DataFrames.DataFrame},
    cat::Dict,
    use_slack,
    storage,
)
    category_dataframes = Dict{String, DataFrames.DataFrame}()
    var_types = Dict([("ThermalMultiStart", :pg)])
    if storage
        var_types["GenericBattery"] = :pb_out
    end
    for (category, list) in cat
        category_df = DataFrames.DataFrame()
        for atuple in list
            if haskey(var_types, atuple[1])
                category_data = data[var_types[atuple[1]]]
                colname =
                    typeof(names(category_data)[1]) == String ? "$(atuple[2])" :
                    Symbol(atuple[2])
                DataFrames.insertcols!(
                    category_df,
                    (colname => category_data[:, colname]),
                    makeunique = true,
                )
            end
        end
        category_dataframes[string(category)] = category_df
    end
    category_dataframes["Wind"] = data[:pW]
    category_dataframes["Hydropower"] = data[:pH]
    category_dataframes["PV"] = data[:pS]
    category_dataframes["Curtailment"] = data[:curt]
    if use_slack
        category_dataframes["Unserved Energy"] = data[:slack_energy⁺]
    end

    return category_dataframes
end

function _get_solar_forecast(
    problem::PSI.OperationsProblem{T},
    time_steps;
    kwargs...,
) where {T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC}}
    scenario = kwargs[:scenario]
    solar_scale = problem.ext["solar_scale"]

    system = PSI.get_system(problem)
    case_initial_time = PSI.get_initial_time(problem)

    area = PSY.get_component(Area, system, "FarWest")
    scenario_forecast = permutedims(
        PSY.get_time_series_values(
            Scenarios,
            area,
            "solar_power";
            start_time = case_initial_time,
        ) .* solar_scale ./ 100,
    )[
        scenario,
        time_steps,
    ]
    return scenario_forecast
end

function _get_solar_forecast(
    problem::PSI.OperationsProblem{T},
    time_steps;
    kwargs...,
) where {T <: BasecaseUnitCommitmentCC}
    solar_scale = problem.ext["solar_scale"]

    forecast =
        get_area_total_time_series(
            problem,
            RenewableGen;
            filter = x -> get_prime_mover(x) == PrimeMovers.PVe && get_available(x),
        )[time_steps] .* solar_scale

    return forecast
end

function _get_solar_realization(
    problem::PSI.OperationsProblem{T},
    time_steps;
    kwargs...,
) where {T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC}}
    scenario = kwargs[:scenario]

    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    solar =
        PSI.axis_array_to_dataframe(jump_model.obj_dict[:pS], [:pS])[time_steps, [scenario]]
    return solar
end

function _get_solar_realization(
    problem::PSI.OperationsProblem{T},
    time_steps;
    kwargs...,
) where {T <: BasecaseUnitCommitmentCC}
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    solar = PSI.axis_array_to_dataframe(jump_model.obj_dict[:pS], [:pS])[time_steps, :]
    return solar
end

function _get_save_path(
    problem::PSI.OperationsProblem{T},
    title,
    format,
    save_dir;
    kwargs...,
) where {T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC}}
    scenario = kwargs[:scenario]
    fname = joinpath(save_dir, "$title Scenario $scenario.$format")
    return fname
end

function _get_save_path(
    problem::PSI.OperationsProblem{T},
    title,
    format,
    save_dir;
    kwargs...,
) where {T <: BasecaseUnitCommitmentCC}
    fname = joinpath(save_dir, "$title.$format")
    return fname
end

# --------------------------------------------------
# Modified plot_fuel for stage 2 using PowerSimulations
# --------------------------------------------------

function my_plot_fuel(res::PSI.SimulationProblemResults, system::PSY.System, storage::Bool; kwargs...)
    save_dir = get(kwargs, :save_dir, nothing)
    time_steps = get(kwargs, :time_steps, nothing)
    use_slack = get(kwargs, :use_slack, true)
    title = get(kwargs, :title, "HAUC Fuel")
    charge_line = get(kwargs, :charge_line, false)

    p = PG._empty_plot()
    backend = Plots.backend()

    timestamps = get_realized_timestamps(res)
    gen = get_generation_data(res)
    if isnothing(time_steps)
        time_steps = 1:length(timestamps)
    end

    # Scale to GW
    for k in keys(gen.data)
        gen.data[k][:, setdiff(names(gen.data[k]), ["DateTime"])] .*=
            get_base_power(system) ./ 1000
    end

    categories = make_fuel_dictionary(system)

    fuel = categorize_data(gen.data, categories; slacks = use_slack, curtailment = true)

    # Hack to make nuclear on the bottom
    cat_names = intersect(PG.CATEGORY_DEFAULT, keys(fuel))
    cat_names = [cat_names[2], cat_names[1], cat_names[3:end]...]
    fuel_agg = PG.combine_categories(fuel; names = cat_names)

    y_label = "Generation (GW)"

    seriescolor = PG.match_fuel_colors(fuel_agg, backend)
    DataFrames.rename!(fuel_agg, Dict("Natural gas" => "NG-CT (solo)"))
    DataFrames.rename!(fuel_agg, Dict("NG-CT" => "NG-CT (in CC train)"))
    p = plot_dataframe(
        fuel_agg,
        gen.time[time_steps];
        seriescolor = seriescolor,
        y_label = y_label,
        title = nothing,
        stack = true,
        set_display = false,
    )

    kwargs = Dict{Symbol, Any}((k, v) for (k, v) in kwargs if k ∉ [:nofill, :seriescolor])
    kwargs[:linestyle] = :dash
    kwargs[:linewidth] = 3

    # Can also try system instead of res, I made up this timestamps call
    load = get_load_data(res; initial_time = first(timestamps))
    load_agg = PG.combine_categories(load.data)
    load_agg .*= get_base_power(system) ./ 1000
    DataFrames.rename!(load_agg, [:Load])

    p = plot_dataframe(
        p,
        load_agg,
        gen.time[time_steps];
        seriescolor = ["black"],
        y_label = y_label,
        title = nothing,
        stack = true,
        nofill = true,
        set_display = false,
        kwargs...,
    )

    # Add charging load line
    if charge_line && storage
        pin = read_realized_variables(res, names = [:Pin__GenericBattery])[:Pin__GenericBattery]
        pin = sum.(eachrow(pin[:, setdiff(names(pin), ["DateTime"])])) .* get_base_power(system) ./ 1000
        load_and_charging = load_agg .+ pin
        DataFrames.rename!(load_and_charging, Symbol.(["Load + charging"]))
        col = PG.match_fuel_colors(fuel_agg[!, ["Storage"]], backend)
        p = plot_dataframe(
            p,
            load_and_charging,
            gen.time;
            seriescolor = [col],
            y_label = y_label,
            title = nothing,
            stack = true,
            nofill = true,
            set_display = false,
            kwargs...,
        )
    end

    # Overwrite x axis label
    layout_kwargs =
        Dict{Symbol, Any}(:xaxis => Plots.PlotlyJS.attr(; title = "Time of Day"))
    Plots.PlotlyJS.relayout!(p, Plots.PlotlyJS.Layout(; layout_kwargs...))

    if !isnothing(save_dir)
        title = replace(title, " " => "_")
        for format in ("png", "pdf")
            fname = joinpath(save_dir, "$title.$format")
            # Overwrite existing plots
            if isfile(fname)
                rm(fname)
            end
            PG.save_plot(p, fname, backend; kwargs...)
        end
    end
    return p
end
