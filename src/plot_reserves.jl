function plot_reserve(
    problem::PSI.OperationsProblem{T}, 
    reserve_name::String;
    use_solar_reserves = true,
    kwargs...) where T <: Union{CVaRReserveUnitCommitmentCC, BasecaseUnitCommitmentCC, StochasticUnitCommitmentCC}
    title = get(kwargs, :title, reserve_name)
    save_dir = get(kwargs, :save_dir, nothing)
    time_steps = get(kwargs, :time_steps, nothing)

    p = PG._empty_plot()
    backend = Plots.backend()

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    case_initial_time = PSI.get_initial_time(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    if isnothing(time_steps)
        time_steps = PSI.model_time_steps(optimization_container)
    end
    use_slack = PSI.get_balance_slack_variables(problem.internal.optimization_container.settings)
    storage = problem.ext["use_storage"]

    sym_dict = Dict{String, Symbol}()
    if reserve_name == "REG_UP"
        reserve = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, reserve_name)
        sym_dict["reserve"] = :reg⁺
        if use_slack sym_dict["slack"] = :slack_reg⁺ end
        if use_solar_reserves sym_dict["solar"] = :reg⁺_S end
        if :supp⁺ in keys(jump_model.obj_dict) sym_dict["supp"] = :supp⁺ end
    elseif reserve_name == "REG_DN"
        reserve = PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, reserve_name)
        sym_dict["reserve"] = :reg⁻
        if use_slack sym_dict["slack"] = :slack_reg⁻ end
        if use_solar_reserves sym_dict["solar"] = :reg⁻_S end
        if :supp⁻ in keys(jump_model.obj_dict) sym_dict["supp"] = :supp⁻ end
    elseif reserve_name == "SPIN"
        reserve = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, reserve_name)
        sym_dict["reserve"] = :spin
        if use_slack sym_dict["slack"] = :spin end
    else
        throw(ArgumentError("Allowable reserve names are REG_UP, REG_DN, or SPIN"))
    end

    device_names = get_name.(get_contributing_devices(system, reserve))
    required_reserve = get_time_series_values(
        Deterministic,
        reserve,
        "requirement";
        start_time = case_initial_time,
    )[time_steps]

    gen = get_reserve_data(
        problem,
        sym_dict,
        time_steps;
        kwargs...
    )
    cat = make_fuel_dictionary(system)

    reserves = my_categorize_reserves(gen.data,
        cat,
        sym_dict,
        storage
    )

    # Hack to make nuclear on the bottom and curtailment on top
    cat_names = intersect(PG.CATEGORY_DEFAULT, keys(reserves))
    cat_names = [cat_names[2], cat_names[1], cat_names[3:end-2]..., cat_names[end], cat_names[end-1]]
    reserves_agg = PG.combine_categories(reserves; names = cat_names)

    y_label = get(kwargs, :y_label, "Reserves (MW)")

    seriescolor = get(kwargs, :seriescolor, PG.match_fuel_colors(reserves_agg, backend))
    if "supp" in keys(sym_dict)
        DataFrames.rename!(reserves_agg, Dict("Imports/Exports" => (reserve_name == "REG_UP" ? "Supp⁺" : "Supp⁻")))
    end
    p = plot_dataframe(
        reserves_agg,
        gen.time;
        seriescolor = seriescolor,
        y_label = y_label,
        title = nothing,
        stack = true,
        set_display = false,
        kwargs...,
    )

    kwargs = Dict{Symbol, Any}((k, v) for (k, v) in kwargs if k ∉ [:nofill, :seriescolor])
    kwargs[:linestyle] = get(kwargs, :linestyle, :dash)
    kwargs[:linewidth] = get(kwargs, :linewidth, 3)

    # Add reserve line
    res_req = DataFrames.DataFrame(Dict(:Requirement => required_reserve)) .*
        get_base_power(system)

    p = plot_dataframe(
        p,
        res_req,
        gen.time;
        seriescolor = ["black"],
        y_label = y_label,
        title = nothing,
        stack = true,
        nofill = true,
        set_display = false,
        kwargs...,
    )

    if !isnothing(save_dir)
        title = replace(title, " " => "_")
        for format in ("png", "pdf")
            fname = _get_reserve_save_path(problem, title, format, save_dir; kwargs...)
            # Overwrite existing plots
            if isfile(fname)
                rm(fname)
            end
            PG.save_plot(p, fname, backend)
        end
    end
    return p
end

function get_reserve_data(
    problem::PSI.OperationsProblem{T},
    sym_dict::Dict,
    time_steps::UnitRange{Int64};
    kwargs...,
) where T <: Union{CVaRReserveUnitCommitmentCC, BasecaseUnitCommitmentCC, StochasticUnitCommitmentCC}
    scenario = kwargs[:scenario]

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)

    variables = Dict{Symbol, DataFrames.DataFrame}()
    variables[sym_dict["reserve"]] = PSI.axis_array_to_dataframe(jump_model.obj_dict[sym_dict["reserve"]], [sym_dict["reserve"]])[time_steps, :]
    if "slack" in keys(sym_dict)
        variables[sym_dict["slack"]] = PSI.axis_array_to_dataframe(jump_model.obj_dict[sym_dict["slack"]], [sym_dict["slack"]])[time_steps, :]
    end

    if "solar" in keys(sym_dict)
        if length(size(jump_model.obj_dict[sym_dict["solar"]])) == 2 # 2D
            variables[sym_dict["solar"]] = PSI.axis_array_to_dataframe(jump_model.obj_dict[sym_dict["solar"]], [sym_dict["solar"]])[time_steps, [scenario]]
        else # 1D
            variables[sym_dict["solar"]] = PSI.axis_array_to_dataframe(jump_model.obj_dict[sym_dict["solar"]], [sym_dict["solar"]])[time_steps, :]
        end
    end

    if "supp" in keys(sym_dict)
        # Supp is 3D transformed to 2D; select single scenario out
        variables[sym_dict["supp"]] = _scenario_in_3D_array_to_dataframe(
            jump_model.obj_dict[sym_dict["supp"]],
            scenario,
            time_steps
        )
    end

    # Scale from 100 MW to MW
    for v in keys(variables)
        variables[v] .*= get_base_power(system)
    end

    timestamps = get_timestamps(problem)[time_steps]
    return PG.PGData(variables, timestamps)
end

function my_categorize_reserves(
    data::Dict{Symbol, DataFrames.DataFrame},
    cat::Dict,
    sym_dict::Dict, 
    storage::Bool
)
    category_dataframes = Dict{String, DataFrames.DataFrame}()
    var_types = Dict([("ThermalMultiStart", sym_dict["reserve"])])
    if storage
        var_types["GenericBattery"] = sym_dict["reserve"]
    end
    for (category, list) in cat
        category_df = DataFrames.DataFrame()
        for atuple in list
            if haskey(var_types, atuple[1])
                category_data = data[var_types[atuple[1]]]
                colname =
                    typeof(names(category_data)[1]) == String ? "$(atuple[2])" :
                    Symbol(atuple[2])
                if colname in names(category_data)
                    DataFrames.insertcols!(
                        category_df,
                        (colname => category_data[:, colname]),
                        makeunique = true,
                    )
                end
            end
        end
        category_dataframes[string(category)] = category_df
    end
    if "solar" in keys(sym_dict)
        category_dataframes["PV"] = data[sym_dict["solar"]]
    end
    if "supp" in keys(sym_dict)
        # Hack to match color, will be renamed
        category_dataframes["Imports/Exports"] = data[sym_dict["supp"]]
    end
    if "slack" in keys(sym_dict)
        category_dataframes["Unserved Energy"] = data[sym_dict["slack"]]
    end

    return category_dataframes
end

function _get_reserve_save_path(
    problem::PSI.OperationsProblem{T},
    title,
    format,
    save_dir;
    kwargs...
    ) where T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC, BasecaseUnitCommitmentCC}
    fname = joinpath(save_dir, "$title.$format")
    return fname
end
