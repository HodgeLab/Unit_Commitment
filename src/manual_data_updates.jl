function apply_manual_data_updates!(system, use_nuclear, initial_cond_file)
    for g in get_components(
        RenewableDispatch,
        system,
        x -> get_prime_mover(x) != PrimeMovers.PVe,
    )
        set_available!(g, true)
    end

    if !use_nuclear
        for g in get_components(
            ThermalMultiStart,
            system,
            x -> get_fuel(x) == ThermalFuels.NUCLEAR,
        )
            set_status!(g, false)
            set_active_power!(g, 0.0)
            set_must_run!(g, false)
        end
    end

    # Overwrite with initial conditions, tailored by day
    initial_on = CSV.read(initial_cond_file, DataFrames.DataFrame)
    for g in PowerSystems.get_components(PowerSystems.ThermalMultiStart, system)
        new_on = Bool(initial_on[1, get_name(g)])
        set_status!(g, new_on)
        if !new_on
            set_active_power!(g, 0.0)
        else
            set_active_power!(
                g,
                (get_active_power_limits(g).max - 10^(-10)) / get_base_power(g) *
                get_base_power(system),
            )
        end
    end
end
