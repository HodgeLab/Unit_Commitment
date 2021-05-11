function appply_manual_data_updates!(system)
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
end