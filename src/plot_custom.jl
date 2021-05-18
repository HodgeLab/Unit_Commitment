function PG.plot_fuel(problem::PSI.OperationsProblem{CVaRUnitCommitmentCC}; kwargs...)
    p = PG._empty_plot()
    set_display = get(kwargs, :set_display, true)
    title = get(kwargs, :title, "Fuel")

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)

    total_load = get_area_total_time_series(problem, PowerLoad)
    total_hydro = get_area_total_time_series(problem, HydroGen)
    total_wind = get_area_total_time_series(
        problem,
        RenewableGen;
        filter = x -> get_prime_mover(x) != PrimeMovers.PVe,
    )

    gen = get_generation_data(problem, total_wind, total_hydro; kwargs...)
    cat = make_fuel_dictionary(system)

    fuel = my_categorize_data(gen.data, cat; kwargs...)

    # Hack to make nuclear on the bottom
    cat_names = intersect(PG.CATEGORY_DEFAULT, keys(fuel))
    cat_names = [cat_names[2], cat_names[1], cat_names[3:end]...]
    fuel_agg = PG.combine_categories(fuel; names = cat_names)

    y_label = get(
        kwargs,
        :y_label,
        "Generation (GW)",
    )

    seriescolor = get(kwargs, :seriescolor, PG.match_fuel_colors(fuel_agg, backend))
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
    load_agg = (PSI.axis_array_to_dataframe(jump_model.obj_dict[:pW], [:pW]) .= total_load) .* 
        get_base_power(system) ./ 1000
    DataFrames.rename!(load_agg, Symbol.(["Load"]))
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

    if set_display
        backend == Plots.PlotlyJSBackend() && Plots.PlotlyJS.plot(p)
        display(p)
    end
    return p
end

function PG.get_generation_data(
    problem::PSI.OperationsProblem{CVaRUnitCommitmentCC},
    total_wind,
    total_hydro;
    kwargs...,
)
    curtailment = get(kwargs, :curtailment, true)
    storage = get(kwargs, :storage, true)

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)

    # Power variable names
    var_names = Vector{Symbol}()
    push!(var_names, :pg, :pS, :pW)
    if storage
        push!(var_names, :pb_in, :pb_out)
    end

    variables = Dict{Symbol, DataFrames.DataFrame}()
    for v in var_names
        variables[v] = PSI.axis_array_to_dataframe(jump_model.obj_dict[v], [v])
        if v == :pg
            Pg = PSI.axis_array_to_dataframe(jump_model.obj_dict[:ug], [:ug])
            for n in names(Pg)
                Pg[!, n] .*=
                    get_active_power_limits(get_component(ThermalMultiStart, system, n)).min
            end
            variables[v] .+= Pg
        end
    end
    # Hack to get shape right
    variables[:pH] =
        PSI.axis_array_to_dataframe(jump_model.obj_dict[:pW], [:pW]) .= total_hydro

    if curtailment
        variables[:pW_curt] =
            PSI.axis_array_to_dataframe(jump_model.obj_dict[:pW], [:pW]) .-= total_wind
    end

    # Scale from 100 MW to GW
    for v in keys(variables)
        variables[v] .*= get_base_power(system) ./ 1000
    end

    timestamps = get_timestamps(problem)
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
    if curtailment
        category_dataframes["Curtailment"] = data[:pW_curt]
    end

    return category_dataframes
end
