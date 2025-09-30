# SPDX-License-Identifier: GPL-3.0-or-later

using Symbolics

function check_formulas(nb_phases::Integer, nb_players::Integer)
    # Define symbols
    @variables B_0 beta_0 phi delta c i j k q r d n

    # Define constants
    d = nb_phases
    n = nb_players

    # Define functions
    B(phi) = B_0 * (1 + cos(phi)) / 2
    beta(phi) = beta_0 * (1 + cos(phi)) / 2

    # Define fixation probabilities
    function rho_CC(phi)
        return 1 / (1 + sum([exp(k * (k + 1 - n) * delta / (n - 1) * (B(phi) - B_0))
                             for k in 1:(n - 1)]))
    end
    function rho_NC(phi)
        return 1 / (1 + sum([exp(delta / (n - 1) *
                                 (j^2 * (beta(phi) - 1 / 2 * B_0) +
                                  j * (1 / 2 * B_0 - (n - 1) * beta(phi) + (n - 1) * c)))
                             for j in 1:(n - 1)]))
    end
    omega = 1 / exp(delta * ((n - 1) * c - (n - 2) / 2 * B_0))
    rho_CN(phi) = rho_NC(phi) / omega

    # Define matrices
    B_1 = [i != j ?
           rho_CC(floor(d / 2) - abs(floor(d / 2) - abs(i - j))) :
           1 - sum([rho_CC(floor(d / 2) - abs(floor(d / 2) - abs(delta_phi_prime)))
                    for delta_phi_prime in 1:(d - 1)])
           for i in 1:d, j in 1:d]
    B_2 = [rho_CN(floor(d / 2) - abs(floor(d / 2) - abs(i - j))) for i in 1:d, j in 1:d]
    B_3 = [rho_NC(floor(d / 2) - abs(floor(d / 2) - abs(i - j))) for i in 1:d, j in 1:d]
    B_4 = [1 / n for i in 1:d, j in 1:d]

    # Define full Markov transition matrix
    M = mortar([[B_1, B_3] [B_2, B_4]])

    # Define stationary distribution
    # Tripp (2024) shows s_1, s_2 independent of q, so f
    q = 1
    s_1 = sum([M[r, q] for r in (d + 1):(2d)]) /
          (d * sum([M[q, r] + M[r, q] for r in (d + 1):(2d)]))
    s_2 = sum([M[q, r] for r in (d + 1):(2d)]) /
          (d * sum([M[q, r] + M[r, q] for r in (d + 1):(2d)]))
    return s = hcat([s_1 * ones(d), s_2 * ones(d)])

    ## Check that the explicitly calculated stationary distribution matches the eigenvector
    #eigen = eigvals(M)[1]
    #if s != eigen
    #    throw(ErrorException(\"Explicitly calculated stationary distribution differs from eigenvector\"))

    ## Check that fraction of communicative simplifies correctly
    #if not isequal(simplify(d*s_1, threaded=true), 1/(1+1/omega))
    #    throw(ErrorException(\"Cooperative fraction does not simplify to known correct answer 1/(1+1/omega)))
end
