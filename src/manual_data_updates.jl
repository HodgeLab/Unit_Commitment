function appply_manual_data_updates!(system, use_nuclear)
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

    # g = get_component(ThermalMultiStart, system, "SANDY_CREEK_J04")
    # set_status!(g, false)
    # set_active_power!(g, 0.0)
end
