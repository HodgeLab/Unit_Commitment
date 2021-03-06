function plot_reserve(
    problem::PSI.OperationsProblem{T},
    reserve_name::String;
    kwargs...,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
        StochasticUnitCommitmentCC,
    },
}
    title = get(kwargs, :title, reserve_name)
    save_dir = get(kwargs, :save_dir, nothing)

    time_steps = get(kwargs, :time_steps, nothing)
    use_solar_reg = problem.ext["use_solar_reg"]
    use_solar_spin = problem.ext["use_solar_spin"]
    use_wind_reserves = problem.ext["use_wind_reserves"]

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

    sym_dict = Dict{String, Symbol}()
    if reserve_name == "REG_UP"
        reserve =
            PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, reserve_name)
        sym_dict["reserve"] = :reg⁺
        if use_slack
            sym_dict["slack"] = :slack_reg⁺
        end
        if use_solar_reg
            sym_dict["solar"] = :reg⁺_S
        end
        if use_wind_reserves
            sym_dict["wind"] = :reg⁺_W
        end
    elseif reserve_name == "REG_DN"
        reserve =
            PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, reserve_name)
        sym_dict["reserve"] = :reg⁻
        if use_slack
            sym_dict["slack"] = :slack_reg⁻
        end
        if use_solar_reg
            sym_dict["solar"] = :reg⁻_S
        end
        if use_wind_reserves
            sym_dict["wind"] = :reg⁻_W
        end
    elseif reserve_name == "SPIN"
        reserve =
            PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, reserve_name)
        sym_dict["reserve"] = :spin
        if use_solar_spin
            sym_dict["solar"] = :spin_S
        end
        if use_wind_reserves
            sym_dict["wind"] = :spin_W
        end
        if use_slack
            sym_dict["slack"] = :slack_spin
        end
        if :total_supp in keys(jump_model.obj_dict) ||
           :total_supp in keys(optimization_container.expressions)
            sym_dict["supp"] = :total_supp
        end
    else
        throw(ArgumentError("Allowable reserve names are REG_UP, REG_DN, or SPIN"))
    end

    required_reserve = get_time_series_values(
        Deterministic,
        reserve,
        "requirement";
        start_time = case_initial_time,
    )[time_steps]

    gen = get_reserve_data(problem, sym_dict, time_steps; kwargs...)
    cat = make_fuel_dictionary(system)

    reserves = my_categorize_reserves(gen.data, cat, sym_dict, storage)

    res_req =
        DataFrames.DataFrame(Dict(:Requirement => required_reserve)) .*
        get_base_power(system) ./ 1000

    p = _plot_reserve_internal(reserves, res_req, gen.time, sym_dict, use_slack; kwargs...)

    if !isnothing(save_dir)
        title = replace(title, " " => "_")
        for format in ("png", "pdf")
            fname = _get_reserve_save_path(problem, title, format, save_dir; kwargs...)
            # Overwrite existing plots
            if isfile(fname)
                rm(fname)
            end
            PG.save_plot(p, fname, Plots.backend(); kwargs...)
        end
    end

    return p
end

function _plot_reserve_internal(
    reserves,
    res_req,
    timestamps,
    sym_dict,
    use_slack;
    kwargs...,
)
    p = PG._empty_plot()
    backend = Plots.backend()

    # Hack to make nuclear on the bottom and curtailment on top
    cat_names = intersect(PG.CATEGORY_DEFAULT, keys(reserves))
    cat_names = [cat_names[2], cat_names[1], cat_names[3:end]...]
    reserves_agg = PG.combine_categories(reserves; names = cat_names)

    y_label = "Reserves (GW)"

    seriescolor = PG.match_fuel_colors(reserves_agg, backend)
    if "supp" in keys(sym_dict)
        DataFrames.rename!(reserves_agg, Dict("Imports/Exports" => "Supplemental"))
    end
    if use_slack
        DataFrames.rename!(reserves_agg, Dict("Unserved Energy" => "Unserved Reserves"))
    end
    DataFrames.rename!(reserves_agg, Dict("Natural gas" => "NG-CT (solo)"))
    DataFrames.rename!(reserves_agg, Dict("NG-CT" => "NG-CT (in CC train)"))
    p = plot_dataframe(
        reserves_agg,
        timestamps;
        seriescolor = seriescolor,
        y_label = y_label,
        title = nothing,
        stack = true,
        set_display = false,
        stair = true,
        kwargs...,
    )

    if !isnothing(res_req)
        kwargs =
            Dict{Symbol, Any}((k, v) for (k, v) in kwargs if k ∉ [:nofill, :seriescolor])
        kwargs[:linestyle] = get(kwargs, :linestyle, :dash)
        kwargs[:linewidth] = get(kwargs, :linewidth, 3)

        # Add reserve line
        p = plot_dataframe(
            p,
            res_req,
            timestamps;
            seriescolor = ["black"],
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

    return p
end

function get_reserve_data(
    problem::PSI.OperationsProblem{T},
    sym_dict::Dict,
    time_steps::UnitRange{Int64};
    kwargs...,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
        StochasticUnitCommitmentCC,
    },
}
    scenario = kwargs[:scenario]

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)

    variables = Dict{Symbol, DataFrames.DataFrame}()
    variables[sym_dict["reserve"]] = PSI.axis_array_to_dataframe(
        jump_model.obj_dict[sym_dict["reserve"]],
        [sym_dict["reserve"]],
    )[
        time_steps,
        :,
    ]
    if "slack" in keys(sym_dict)
        variables[sym_dict["slack"]] = PSI.axis_array_to_dataframe(
            jump_model.obj_dict[sym_dict["slack"]],
            [sym_dict["slack"]],
        )[
            time_steps,
            :,
        ]
    end

    if "solar" in keys(sym_dict)
        if length(size(jump_model.obj_dict[sym_dict["solar"]])) == 2 # 2D
            variables[sym_dict["solar"]] = PSI.axis_array_to_dataframe(
                jump_model.obj_dict[sym_dict["solar"]],
                [sym_dict["solar"]],
            )[
                time_steps,
                [scenario],
            ]
        else # 1D
            variables[sym_dict["solar"]] = PSI.axis_array_to_dataframe(
                jump_model.obj_dict[sym_dict["solar"]],
                [sym_dict["solar"]],
            )[
                time_steps,
                :,
            ]
        end
    end

    if "wind" in keys(sym_dict)
        variables[sym_dict["wind"]] = PSI.axis_array_to_dataframe(
            jump_model.obj_dict[sym_dict["wind"]],
            [sym_dict["wind"]],
        )[
            time_steps,
            :,
        ]
    end

    if "supp" in keys(sym_dict)
        # Supp is 3D transformed to 2D; select single scenario out
        if problem.ext["supp_type"] == "generic" # generic version
            variables[sym_dict["supp"]] =
                PSI.axis_array_to_dataframe(jump_model.obj_dict[:total_supp], [:total_supp])[
                    time_steps,
                    [scenario],
                ]
        else # Nonspin version, total_supp in obj_dict
            variables[sym_dict["supp"]] = _scenario_in_3D_array_to_dataframe(
                jump_model.obj_dict[:supp],
                scenario,
                time_steps,
            )
        end
    end

    # Scale from 100 MW to GW
    for v in keys(variables)
        variables[v] .*= get_base_power(system) ./ 1000
    end

    timestamps = get_timestamps(problem)[time_steps]
    return PG.PGData(variables, timestamps)
end

function my_categorize_reserves(
    data::Dict{Symbol, DataFrames.DataFrame},
    cat::Dict,
    sym_dict::Dict,
    storage::Bool,
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
    if "wind" in keys(sym_dict)
        category_dataframes["Wind"] = data[sym_dict["wind"]]
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
    kwargs...,
) where {T <: Union{CVaRReserveUnitCommitmentCC, StochasticUnitCommitmentCC}}
    scenario = kwargs[:scenario]
    fname = joinpath(save_dir, "$title Scenario $scenario.$format")
    return fname
end

function _get_reserve_save_path(
    problem::PSI.OperationsProblem{T},
    title,
    format,
    save_dir;
    kwargs...,
) where {T <: BasecaseUnitCommitmentCC}
    fname = joinpath(save_dir, "$title.$format")
    return fname
end

################################# Stage 2 ###############################################

function plot_stage2_reserves(
    res::PSI.SimulationProblemResults,
    system::PSY.System,
    reserve_name::String;
    kwargs...,
)
    title = get(kwargs, :title, reserve_name)
    save_dir = get(kwargs, :save_dir, nothing)
    use_slack = get(kwargs, :use_slack, true)

    sym_dict = Dict{String, Symbol}()
    if reserve_name == "REG_UP"
        reserve_component =
            PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, reserve_name)
        sym_dict["reserve"] = :REG_UP__VariableReserve_ReserveUp
        if use_slack
            sym_dict["slack"] = :γ⁺__REG_UP
        end
    elseif reserve_name == "REG_DN"
        reserve_component =
            PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, reserve_name)
        sym_dict["reserve"] = :REG_DN__VariableReserve_ReserveDown
        if use_slack
            sym_dict["slack"] = :γ⁺__REG_DN
        end
    elseif reserve_name == "SPIN"
        reserve_component =
            PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, reserve_name)
        sym_dict["reserve"] = :SPIN__VariableReserve_ReserveUp
        if use_slack
            sym_dict["slack"] = :γ⁺__SPIN
        end
    else
        throw(ArgumentError("Allowable reserve names are REG_UP, REG_DN, or SPIN"))
    end

    gen = get_reserve_data(res, sym_dict, system)
    categories = make_fuel_dictionary(system)

    reserves = my_categorize_stage2_reserves(gen.data, categories, sym_dict)

    timestamps = get_realized_timestamps(res)
    hour_timestamps = collect(first(timestamps):Hour(1):last(timestamps))

    required_reserve = zeros(length(timestamps))
    for i in 1:length(hour_timestamps)
        required_reserve[((i - 1) * 12 + 1):(i * 12)] =
            get_time_series_values(
                Deterministic,
                reserve_component,
                "requirement";
                start_time = hour_timestamps[i],
                len = 12,
            ) .* get_base_power(system) ./ 1000
    end
    res_req = DataFrames.DataFrame(Dict(:Requirement => required_reserve))

    p = _plot_reserve_internal(reserves, res_req, gen.time, sym_dict, use_slack; kwargs...)

    if !isnothing(save_dir)
        title = replace(title, " " => "_")
        for format in ("png", "pdf")
            fname = joinpath(save_dir, "$title.$format")
            # Overwrite existing plots
            if isfile(fname)
                rm(fname)
            end
            PG.save_plot(p, fname, Plots.backend(); kwargs...)
        end
    end

    return p
end

function get_reserve_data(
    res::PSI.SimulationProblemResults,
    sym_dict::Dict,
    system::PSY.System,
)
    variables = Dict{Symbol, DataFrames.DataFrame}()
    variables[sym_dict["reserve"]] =
        read_realized_variables(res, names = [sym_dict["reserve"]])[sym_dict["reserve"]]
    if "slack" in keys(sym_dict)
        variables[sym_dict["slack"]] =
            read_realized_variables(res, names = [sym_dict["slack"]])[sym_dict["slack"]]
    end

    # Scale from 100 MW to GW
    for v in keys(variables)
        variables[v][:, setdiff(names(variables[v]), ["DateTime"])] .*=
            get_base_power(system) ./ 1000
    end

    timestamps = get_realized_timestamps(res)
    return PG.PGData(variables, timestamps)
end

function my_categorize_stage2_reserves(
    data::Dict{Symbol, DataFrames.DataFrame},
    categories::Dict,
    sym_dict::Dict,
)
    category_dataframes = Dict{String, DataFrames.DataFrame}()
    var_types = Dict([("ThermalMultiStart", sym_dict["reserve"])])
    var_types["GenericBattery"] = sym_dict["reserve"]
    var_types["RenewableDispatch"] = sym_dict["reserve"]
    for (category, list) in categories
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
    if "slack" in keys(sym_dict)
        category_dataframes["Unserved Energy"] = data[sym_dict["slack"]]
    end

    return category_dataframes
end
