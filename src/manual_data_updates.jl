function set_storage_reserve_SOC_to_max!(system, storage_reserve_names)
    for stor in storage_reserve_names
        g = get_component(GenericBattery, system, stor)
        set_initial_energy!(g, get_state_of_charge_limits(g)[:max])
    end
end

function set_initial_SOC!(system, initial_soc)
    storage_names = PSY.get_name.(get_components(PSY.GenericBattery, system))
    for stor in storage_names
        g = get_component(GenericBattery, system, stor)
        set_initial_energy!(g, get_state_of_charge_limits(g)[:max] * initial_soc)
    end
end

# Updates for stage 1
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
        new_on = Bool(round(initial_on[1, get_name(g)]))
        set_status!(g, new_on)
        if !new_on
            set_active_power!(g, 0.0)
        else
            set_active_power!(g, get_active_power_limits(g).max)
        end
    end
end

# Updates for stage 2
function add_inverter_based_reserves!(
    system,
    use_solar_reg,
    use_solar_spin,
    use_storage_reserves,
    storage_reserve_names,
)
    # Add solar and the selected battery to the contributing device set
    reg_reserve_up = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "REG_UP")
    reg_reserve_dn =
        PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, "REG_DN")
    spin_reserve = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "SPIN")
    for g in get_components(
        RenewableDispatch,
        system,
        x -> get_prime_mover(x) == PrimeMovers.PVe,
    )
        if use_solar_reg
            add_service!(g, reg_reserve_up, system)
            add_service!(g, reg_reserve_dn, system)
        end
        if use_solar_spin
            add_service!(g, spin_reserve, system)
        end
    end
    if use_storage_reserves
        for stor in storage_reserve_names
            add_service!(get_component(GenericBattery, system, stor), spin_reserve, system)
        end
    end
end

function add_to_reserve_contributing_devices!(
    system
)
    reg_reserve_up = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "REG_UP")
    reg_reserve_dn =
        PSY.get_component(PSY.VariableReserve{PSY.ReserveDown}, system, "REG_DN")
    spin_reserve = PSY.get_component(PSY.VariableReserve{PSY.ReserveUp}, system, "SPIN")

    # Relax regulation response time to 100 hours to make non-binding
    set_time_frame!(reg_reserve_up, 6000.0)
    set_time_frame!(reg_reserve_dn, 6000.0)

    default_reg⁺_device_names = get_name.(get_contributing_devices(system, reg_reserve_up))
    default_reg⁻_device_names = get_name.(get_contributing_devices(system, reg_reserve_dn))
    default_spin_device_names = get_name.(get_contributing_devices(system, spin_reserve))
    desired_reg⁺_device_names =
        get_name.(get_components(ThermalMultiStart, system, x -> !PSY.get_must_run(x)))
    desired_reg⁻_device_names =
        get_name.(get_components(ThermalMultiStart, system, x -> !PSY.get_must_run(x)))
    desired_spin_device_names =
        get_name.(get_components(ThermalMultiStart, system, x -> !PSY.get_must_run(x)))

    for g in get_components(
        ThermalMultiStart,
        system,
        x -> get_name(x) in desired_reg⁺_device_names && !(get_name(x) in default_reg⁺_device_names),
    )
        add_service!(g, reg_reserve_up, system)
    end
    for g in get_components(
        ThermalMultiStart,
        system,
        x -> get_name(x) in desired_reg⁻_device_names && !(get_name(x) in default_reg⁻_device_names),
    )
        add_service!(g, reg_reserve_dn, system)
    end
    for g in get_components(
        ThermalMultiStart,
        system,
        x -> get_name(x) in desired_spin_device_names && !(get_name(x) in default_spin_device_names),
    )
        add_service!(g, spin_reserve, system)
    end
end
