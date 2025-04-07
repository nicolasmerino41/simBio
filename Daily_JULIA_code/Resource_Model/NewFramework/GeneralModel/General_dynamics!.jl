function general_dynamics!(du, u, p, t)
    # Unpack parameters:
    S, R, K_i, r_i, mu, nu, A, C, m_alpha, K_alpha = p
    
    # Partition the state vector.
    H = @view u[1:S]
    P = @view u[S+1:S+R]
    X = vcat(H, P)  # Combined state vector

    # Herbivore (and omnivore) dynamics:
    # Compute the total trophic effect on herbivores using the merged matrix.
    # effective_interaction = (A[1:S, :] * X) # We've embbedded nu into A already
    # # Compute competition among herbivores.
    # competition_effect = C * H # We've embedded mu into C already
    # Total effective density experienced by each herbivore.
    # effective_density = H .+ effective_interaction .+ competition_effect
    effective_density = (1 - mu) .* H .+ (C * H) .+ (A[1:S, :] * X)

    # Compute herbivore growth rate.
    duH = H .* r_i .* (1 .- effective_density ./ K_i)

    # Predator dynamics:
    effective_interaction_pred = A[S+1:S+R, :] * X
    duP = P .* m_alpha .* ((effective_interaction_pred .- P) ./ K_alpha .- 1)

    # Write the derivatives into du.
    du[1:S] .= duH
    du[S+1:S+R] .= duP
end
