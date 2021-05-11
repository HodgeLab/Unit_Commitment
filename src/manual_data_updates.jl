function appply_manual_data_updates!(system)
    for g in get_components(
        ThermalMultiStart,
        system,
        x -> get_prime_mover(x) == ThermalFuels.NATURAL_GAS,
    )
        lims = get_ramp_limits(g)
        set_time_limits!(g, (up = 1.1 * lims.up, down = 1.1 * lims.down))
        set_ramp_limits!(g, (up = 1.25 * lims.up, down = 1.25 * lims.down))
    end

    for g in get_components(
        ThermalMultiStart,
        system,
        x -> get_prime_mover(x) == ThermalFuels.COAL,
    )
        lims = get_ramp_limits(g)
        set_ramp_limits!(g, (up = 1.25 * lims.up, down = 1.25 * lims.down))
    end
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

    # Removing additional reg requirements because unneeded for this
    # s = get_component(VariableReserve{ReserveUp}, system, "REG_UP")
    # req = get_requirement(s)
    # set_requirement!(s, req * 1.5)

    # s = get_component(VariableReserve{ReserveDown}, system, "REG_DN")
    # req = get_requirement(s)
    # set_requirement!(s, req * 1.5)
end
