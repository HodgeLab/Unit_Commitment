function apply_manual_data_updates!(system, use_nuclear, system_file_path)
    for g in get_components(
        RenewableDispatch,
        system,
        x -> get_prime_mover(x) != PrimeMovers.PVe,
    )
        set_available!(g, true)
    end

    for g in get_components(HydroGen, system)
        set_available!(g, true)
    end

    # Set all CC's to start off
    for g in get_components(
        ThermalMultiStart,
        system,
        x -> get_prime_mover(x) in [PrimeMovers.CT, PrimeMovers.CC],
    )
        set_status!(g, false)
        set_active_power!(g, 0.0)
    end

    if !use_nuclear
        for g in get_components(
            ThermalMultiStart,
            system,
            x -> get_fuel(x) == ThermalFuels.NUCLEAR
        )
            set_status!(g, false)
            set_active_power!(g, 0.0)
            set_must_run!(g, false)
        end
    end

    # Overwrite with selected inital on conditions, from running 5/17/2018
    initial_on = CSV.read(joinpath(system_file_path, "initial_on.csv"), DataFrame)
    for g in get_components(
        ThermalMultiStart,
        system
    )
        new_on = Bool(initial_on[1, get_name(g)])
        set_status!(g, new_on)
        if !new_on
            set_active_power!(g, 0.0)
        else
            set_active_power!(g, (get_active_power_limits(g).max-10^(-10))/get_base_power(g)*get_base_power(system))
        end
    end
end
