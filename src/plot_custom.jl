function PG.plot_fuel(problem::PSI.OperationsProblem{CVaRUnitCommitmentCC}; kwargs...)
    title = get(kwargs, :title, "Fuel")
    save_dir = get(kwargs, :save_dir, nothing)
    storage = get(kwargs, :storage, true)
    scenario = get(kwargs, :scenario, 1)
    case_initial_time = get(kwargs, :case_initial_time, nothing)
    time_steps = get(kwargs, :time_steps, nothing)

    p = PG._empty_plot()
    backend = Plots.backend()

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    if isnothing(time_steps)
        time_steps = PSI.model_time_steps(optimization_container)
    end

    total_load = get_area_total_time_series(problem, PowerLoad)[time_steps]
    total_hydro = get_area_total_time_series(problem, HydroGen)[time_steps]
    total_wind = get_area_total_time_series(
        problem,
        RenewableGen;
        filter = x -> get_prime_mover(x) != PrimeMovers.PVe,
    )[time_steps]

    area = PSY.get_component(Area, system, "1")
    scenario_forecast = permutedims(
        PSY.get_time_series_values(
            Scenarios,
            area,
            "solar_power";
            start_time = case_initial_time,
        ) ./ 100,
    )[
        scenario,
        time_steps,
    ]

    gen = get_generation_data(
        problem,
        total_wind,
        total_hydro,
        scenario,
        scenario_forecast,
        time_steps;
        kwargs...,
    )
    cat = make_fuel_dictionary(system)

    fuel = my_categorize_data(gen.data, cat; kwargs...)

    # Hack to make nuclear on the bottom and curtailment on top
    cat_names = intersect(PG.CATEGORY_DEFAULT, keys(fuel))
    cat_names = [cat_names[2], cat_names[1], cat_names[3:end-2]..., cat_names[end], cat_names[end-1]]
    fuel_agg = PG.combine_categories(fuel; names = cat_names)

    y_label = get(kwargs, :y_label, "Generation (GW)")

    seriescolor = get(kwargs, :seriescolor, PG.match_fuel_colors(fuel_agg, backend))
    DataFrames.rename!(fuel_agg, Dict("Imports/Exports" => "Supp⁺"))
    p = plot_dataframe(
        fuel_agg,
        gen.time;
        seriescolor = seriescolor,
        y_label = y_label,
        title = title,
        stack = true,
        set_display = false,
        kwargs...,
    )

    kwargs = Dict{Symbol, Any}((k, v) for (k, v) in kwargs if k ∉ [:nofill, :seriescolor])
    kwargs[:linestyle] = get(kwargs, :linestyle, :dash)
    kwargs[:linewidth] = get(kwargs, :linewidth, 3)

    # Add load line
    load_agg = DataFrames.DataFrame(Dict(:Load => total_load)) .*
        get_base_power(system) ./ 1000
    p = plot_dataframe(
        p,
        load_agg,
        gen.time;
        seriescolor = ["black"],
        y_label = y_label,
        title = title,
        stack = true,
        nofill = true,
        set_display = false,
        kwargs...,
    )

    # Add charging load line
    if storage
        load_and_charging = load_agg .+ sum.(eachrow(gen.data[:pb_in]))
        DataFrames.rename!(load_and_charging, Symbol.(["Load + charging"]))
        col = PG.match_fuel_colors(fuel_agg[!, ["Storage"]], backend)
        p = plot_dataframe(
            p,
            load_and_charging,
            gen.time;
            seriescolor = [col],
            y_label = y_label,
            title = title,
            stack = true,
            nofill = true,
            set_display = false,
            kwargs...,
        )
    end

    if !isnothing(save_dir)
        title = replace(title, " " => "_")
        format = get(kwargs, :format, "png")
        fname = joinpath(save_dir, "$title Scenario $scenario.$format")
        # Overwrite existing plots
        if isfile(fname)
            rm(fname)
        end
        PG.save_plot(p, fname, backend)
    end
    return p
end

function PG.get_generation_data(
    problem::PSI.OperationsProblem{CVaRUnitCommitmentCC},
    total_wind,
    total_hydro,
    scenario,
    scenario_forecast,
    time_steps;
    kwargs...,
)
    curtailment = get(kwargs, :curtailment, true)
    storage = get(kwargs, :storage, true)

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)

    # Power variable names
    var_names = Vector{Symbol}()
    push!(var_names, :pg)
    if storage
        push!(var_names, :pb_in, :pb_out)
    end

    variables = Dict{Symbol, DataFrames.DataFrame}()
    for v in var_names
        variables[v] = PSI.axis_array_to_dataframe(jump_model.obj_dict[v], [v])[time_steps, :]
        if v == :pg
            Pg = PSI.axis_array_to_dataframe(jump_model.obj_dict[:ug], [:ug])[time_steps, :]
            for n in names(Pg)
                Pg[!, n] .*=
                    get_active_power_limits(get_component(ThermalMultiStart, system, n)).min
            end
            variables[v] .+= Pg
        end
    end
    # Select single solar scenario
    if curtailment
        variables[:pS] =
        PSI.axis_array_to_dataframe(jump_model.obj_dict[:pS], [:pS])[time_steps, [scenario]]
        variables[:pW] =
            PSI.axis_array_to_dataframe(jump_model.obj_dict[:pW], [:pW])[time_steps, :]
    else
        variables[:pW] = DataFrames.DataFrame(Dict(:pW => total_wind))
        variables[:pS] = DataFrames.DataFrame(Dict(:pW => scenario_forecast))
    end
    variables[:pH] = DataFrames.DataFrame(Dict(:pH => total_hydro))

    # Supp is 3D transformed to 2D; select single scenario out
    x = PSI.axis_array_to_dataframe(jump_model.obj_dict[:supp⁺], [:supp⁺])
    variables[:supp⁺] = (x[x[!, :S1] .== scenario, names(x) .!= "S1"])[time_steps, :]

    variables[:curt] = DataFrames.DataFrame(
        :curt =>
            total_wind - variables[:pW][!, 1] + scenario_forecast -
            variables[:pS][!, 1],
    )

    # Scale from 100 MW to GW
    for v in keys(variables)
        variables[v] .*= get_base_power(system) ./ 1000
    end

    timestamps = get_timestamps(problem)[time_steps]
    return PG.PGData(variables, timestamps)
end

function my_categorize_data(
    data::Dict{Symbol, DataFrames.DataFrame},
    cat::Dict;
    curtailment = true,
    slacks = true,
    kwargs...,
)
    storage = get(kwargs, :storage, true)
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
    # Hack to match color, will be renamed
    category_dataframes["Imports/Exports"] = data[:supp⁺]
    category_dataframes["Curtailment"] = data[:curt]

    return category_dataframes
end
