"""
    find_critical_press(p, B_eq;
                        δ_lo=0.0, δ_hi=1.0,
                        pers_target=1.0,
                        tol=1e-3,
                        tspan=(0.,500.),
                        t_perturb=250.0,
                        cb=cb)

Given parameter tuple `p` and equilibrium `B_eq`, find by bisection
the largest relative press δ in [δ_lo,δ_hi] such that, when we
bump every consumer threshold by (1+δ)*ξ, the community still
retains fraction ≥ `pers_target` of its species.

Returns the critical δ (within `tol`) at which persistence just dips
below `pers_target`.
"""
function find_critical_press(p, B_eq;
                             δ_lo=0.0, δ_hi=1.0,
                             pers_target=1.0,
                             tol=1e-3,
                             tspan=(0.,500.), t_perturb=250.0,
                             cb)
    R, C = p[1], p[2]

    # helper: compute persistence fraction under relative press δ
    simulate_pers = δ -> begin
        # build relative press
        xi_cons = p[4]
        xi2     = xi_cons .* (1 .+ δ)
        p2      = (R, C, p[3], xi2, p[5], p[6], p[7], p[8])

        _,_,_, pers, _, _ = simulate_press_perturbation(
            B_eq, p2, tspan, t_perturb, 0.0;
            full_or_simple=true,
            cb=cb,
            species_specific_perturbation=false
        )
        return pers
    end

    # ensure that at δ_lo we are safe, at δ_hi we have dropped below
    pers_lo = simulate_pers(δ_lo)
    pers_hi = simulate_pers(δ_hi)
    if pers_lo < pers_target
        error("Even at δ_lo=$δ_lo persistence = $pers_lo < $pers_target")
    end

    # expand δ_hi until persistence < target
    while pers_hi ≥ pers_target
        δ_hi *= 2
        pers_hi = simulate_pers(δ_hi)
        if δ_hi > 1e3
            error("Could not find upper bound for critical δ")
        end
    end

    # now bisection in [δ_lo, δ_hi]
    while δ_hi - δ_lo > tol
        δ_mid = (δ_lo + δ_hi)/2
        if simulate_pers(δ_mid) ≥ pers_target
            δ_lo = δ_mid
        else
            δ_hi = δ_mid
        end
    end

    return δ_lo   # largest δ that still meets persistence ≥ pers_target
end
cb_no_trigger36 = build_callbacks(36, EXTINCTION_THRESHOLD)
find_critical_press(p, B_eq; cb = cb_no_trigger36)