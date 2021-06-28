function plot_charging(
    problem::PSI.OperationsProblem{T};
    kwargs...,
) where {
    T <: Union{
        CVaRReserveUnitCommitmentCC,
        BasecaseUnitCommitmentCC,
        StochasticUnitCommitmentCC,
    },
}
    title = get(kwargs, :title, "charging")
    save_dir = get(kwargs, :save_dir, nothing)
    time_steps = get(kwargs, :time_steps, nothing)

    system = PSI.get_system(problem)
    optimization_container = PSI.get_optimization_container(problem)
    jump_model = PSI.get_jump_model(optimization_container)
    if isnothing(time_steps)
        time_steps = PSI.model_time_steps(optimization_container)
    end

    pin =
        sum.(eachrow(PSI.axis_array_to_dataframe(jump_model.obj_dict[:pb_in], [:pb_in])[time_steps, :] .*
        get_base_power(system) ./ 1000))
    pout =
        sum.(eachrow(PSI.axis_array_to_dataframe(jump_model.obj_dict[:pb_out], [:pb_out])[time_steps, :] .*
        get_base_power(system) ./ 1000))
    timestamps = get_timestamps(problem)[time_steps]
    net_pout = pout .- pin

    trace1 = Plots.PlotlyJS.scatter(; x = timestamps, y = net_pout,
        mode = "lines", line_shape = "hv"
    )

    plot_kwargs =
    Dict{Symbol, Any}(((k, v) for (k, v) in kwargs if k in PG.SUPPORTED_EXTRA_PLOT_KWARGS))

    layout = Plots.PlotlyJS.Layout(;
    xaxis=Plots.PlotlyJS.attr(title="Time of Day"),
    yaxis=Plots.PlotlyJS.attr(title="Discharge Power (GW)"),)
    p = Plots.PlotlyJS.plot([trace1], layout; plot_kwargs...)

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

################################# Stage 2 ###############################################

function plot_charging(
    res::PSI.SimulationProblemResults,
    system::PSY.System;
    kwargs...,
)
    title = get(kwargs, :title, "charging")
    save_dir = get(kwargs, :save_dir, nothing)

    pout = read_realized_variables(res, names = [:Pout__GenericBattery])[:Pout__GenericBattery]
    timestamps = pout[!, :DateTime]
    pout = sum.(eachrow(pout[:, setdiff(names(pout), ["DateTime"])])) .* get_base_power(system) ./ 1000

    pin = read_realized_variables(res, names = [:Pin__GenericBattery])[:Pin__GenericBattery]
    pin = sum.(eachrow(pin[:, setdiff(names(pin), ["DateTime"])])) .* get_base_power(system) ./ 1000
    net_pout = pout .- pin

    trace1 = Plots.PlotlyJS.scatter(; x = timestamps, y = net_pout,
        mode = "lines", line_shape = "hv"
    )

    plot_kwargs =
    Dict{Symbol, Any}(((k, v) for (k, v) in kwargs if k in PG.SUPPORTED_EXTRA_PLOT_KWARGS))

    layout = Plots.PlotlyJS.Layout(;
    xaxis=Plots.PlotlyJS.attr(title="Time of Day"),
    yaxis=Plots.PlotlyJS.attr(title="Discharge Power (GW)"),)
    p = Plots.PlotlyJS.plot([trace1], layout; plot_kwargs...)

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
