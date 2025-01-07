function generate_competition_matrix(S::Int, mu::Float64, symmetric::Bool; check_condition=true, condition_limit_number=100.0)
    # This function returns (V, μ_matrix) given S, mu, and whether we want a symmetric scenario or not.
    #
    # If symmetric=true:
    #   (I+μ) has diagonal = 1.0, off-diagonal = μ (perfectly uniform scenario).
    #
    # If symmetric=false:
    #   (I+μ) has diagonal = 1.0, off-diagonal random in [0,1], then scaled so mean=μ.
    #   If condition number > condition_limit_number, reduce asymmetry by blending off-diagonal values towards μ.
    #
    # After constructing (I+μ), we invert it to get V, and compute μ_matrix = (I+μ)-I.
    # Check condition number if requested.

    I_plus_mu = Matrix{Float64}(undef, S, S)

    if symmetric
        # Symmetric scenario: diagonal=1.0, off-diagonal=μ
        for i in 1:S, j in 1:S
            I_plus_mu[i,j] = (i == j) ? 1.0 : mu
        end
    else
        # Asymmetric scenario:
        # diagonal=1.0
        # off-diagonal: random values scaled to have mean=μ
        for i in 1:S, j in 1:S
            if i == j
                I_plus_mu[i,j] = 1.0
            else
                I_plus_mu[i,j] = rand()
            end
        end

        # Scale off-diagonal to have mean=μ
        off_diag_vals = [I_plus_mu[i,j] for i in 1:S, j in 1:S if i!=j]
        current_mean = mean(off_diag_vals)
        if current_mean != 0.0
            scaling_factor = mu / current_mean
            for i in 1:S, j in 1:S
                if i != j
                    I_plus_mu[i,j] *= scaling_factor
                end
            end
        else
            # If current_mean=0 (extremely unlikely), set all off-diag to μ
            for i in 1:S, j in 1:S
                if i != j
                    I_plus_mu[i,j] = mu
                end
            end
        end
    end

    if check_condition
        # Check condition number and if too large, reduce asymmetry if not symmetric
        max_attempts = 5
        attempts = 0
        while attempts < max_attempts
            cnum = cond(I_plus_mu)
            println("Condition number of (I+μ): ", cnum)
            if cnum <= condition_limit_number
                break
            end

            if symmetric
                # If symmetric and condition is large, there's not much we can do without changing parameters.
                println("Warning: Condition number still high in symmetric scenario. Consider changing μ or S.")
                break
            else
                # Reduce asymmetry by bringing off-diagonal values closer to μ
                # e.g., blend each off-diag value with μ: new_val = (val + μ)/2
                println("Condition > $condition_limit_number, reducing asymmetry by blending towards μ...")
                for i in 1:S, j in 1:S
                    if i != j
                        I_plus_mu[i,j] = (I_plus_mu[i,j] + mu)/2
                    end
                end
            end

            attempts += 1
        end
    end

    println("Final condition number of (I+μ): ", cond(I_plus_mu))
    if cond(I_plus_mu) > condition_limit_number && !symmetric
        println("Warning!!!!!!!!!: Condition number still high in asymmetric scenario. Consider changing μ or S.")
    end

    # Invert I+μ
    # Compute V = (I+μ)^{-1}
    V = inv(I_plus_mu)

    # Compute μ_matrix = (I+μ)-I
    μ_matrix = Matrix{Float64}(undef, S, S)
    for i in 1:S, j in 1:S
        if i == j
            μ_matrix[i,j] = I_plus_mu[i,j] - 1.0
        else
            μ_matrix[i,j] = I_plus_mu[i,j]
        end
    end

    return V, μ_matrix
end