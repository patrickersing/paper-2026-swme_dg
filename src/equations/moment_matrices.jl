
# Compute the necessary tensors for the shallow water moment equations
# in the case where they have not been linearized.

# Usage hints
#
# using StaticArrays
# n = 4
# A = compute_A_tensor(n)
# B = compute_B_tensor(n)
# C = compute_C_matrix(n)

using StaticArrays: MArray, SArray

# Compute the Gaunt coefficients for the Legendre triple product via
# a closed form strategy for the scaled Legendre polynomials shifted
# to the interval [0, 1]. Said shifted and scaled polynomials satisfy
# the recursion relation
#    phi_0 = 1, phi_1 = 1 - 2 * x
#    (j + 1) * phi_{j+1} = (2 * j + 1) * (1 - 2 * x) * phi_j - j * phi_{j-1}
# for j = 1, 2, ..., N
# This, effectively, computes the Wigner-3j coefficient in a numerically efficient way.
#
# The explicit formula comes from the (surprisingly easy to understand) classic reference
# - John Couch Adams (1878)
#   On the expression of the product of any two Legendre’s coefficients by means
#   of a series of Legendre’s coefficients
#   [DOI:10.1098/rspl.1878.0016](https://doi.org/10.1098/rspl.1878.0016)
#
# In this case we are computing
#   A_ijk = (2i + 1) integral_0^1 phi_i phi_j phi_k  d zeta
#   i, j, k = 1, ..., n
#
# OBS! Adjusted the indexing to match how Julian does it. The standard form would
# include the zero index
function compute_A_tensor(n::Int, RealT = Float64)
    # save preliminary values in pre-allocated memory
    A = zero(MArray{Tuple{n, n, n}, RealT})

    # build the three term inner products
    for i in 1:n, j in 1:n, k in 1:n
        # Use integer division so that the variable type of `s`, `i`, `j`, and `k` all match
        s = div(i + j + k, 2)
        if rem(i + j + k, 2) == 1 || abs(i - j) > k || k > i + j
            A[i, j, k] = 0
        else
            # assumes unnormalized basis functions
            scaling = (2 * i + 1) * (-1)^(i + j + k) / (i + j + k + 1)
            A[i, j, k] = scaling * Adams_A(s - i) * Adams_A(s - j) * Adams_A(s - k) /
                         Adams_A(s)
        end
    end # i, j, k

    return SArray{Tuple{n, n, n}, RealT}(A)
end

# Auxiliary function for triple product matrix construction
# that computes the quantity A(n) = (2n-1)!! / n! from the
# Adams paper with n >= 0. Here we call the function Adams_A
# to avoid name conflicts.
function Adams_A(n::Int)
    # Initial value fulfills the n == 0 branch so no need to code it
    A = 1

    if n < 0
        A = 0
    end

    if n >= 1
        for j in 1:n
            A *= (2 * j - 1)
        end # j
        A /= factorial(n)
    end

    return A
end

# Create the B tensor for the moment equations in a (somewhat) clever way.
# We use recursion relations for the Legendre polynomials and their derivative
# to build the necessary components of the integral
#
#   B_ijk = (2i + 1) integral_0^1 phi'_i (integral_0^zeta phi_j ds) phi_k  d zeta
#   i, j, k = 1, ..., N
#
# where phi'_i is a polynomial of degree at most `N-1`
# and \phi_j and \phi_k are polynomials of degree at most `N`.
# Therefore, `integral_0^zeta phi_j ds` is a polynomial of at most `N+1`.
# So, the highest degree of the integrand is a polynomial
# of degree `(N-1) + (N+1) + N = 3N`. We approximate the integrals with
# Legendre-Gauss quadrature which is exact for polynomials up to degree `2M+1`.
# Thus, if we select `M = ceil( (3N - 1) / 2 ) + 1` nodes, then the quadrature
# will be exact for the construction of the B tensor.
function compute_B_tensor(n::Int, RealT = Float64)
    # save preliminary values in pre-allocated memory
    B = zero(MArray{Tuple{n, n, n}, RealT})

    # Given the number of moments `n` determine how many Legendre-Gauss
    # nodes are required and compute them
    M = Int(ceil(1.5f0 * n - 0.5f0) + 1)
    lg_nodes, lg_weights = gauss_nodes_weights(M, RealT)

    # Precompute the values of \phi and \phi' evaluated at the mapped LG nodes
    phi_values = Array{RealT, 2}(undef, n, M)
    phi_prime_values = Array{RealT, 2}(undef, n, M)
    for j in 1:n, m in 1:M
        x = 0.5f0 * (lg_nodes[m] + 1)
        phi_values[j, m], phi_prime_values[j, m] = shifted_legendre_polynomial_and_derivative(j,
                                                                                              x)
    end # j, m

    # Precompute all the internal integrals
    phi_j_integrals = Array{RealT, 2}(undef, n, M)
    compute_inner_integrals!(phi_j_integrals, lg_nodes, lg_weights)

    # Compute everything using the precomputed work arrays
    for i in 1:n, j in 1:n, k in 1:n
        # Legendre-Gauss quadrature loop for the outer integral
        res = 0
        for m in 1:M
            # Outer LG quadrature node and weight
            w_m = lg_weights[m]

            res += phi_prime_values[i, m] * phi_j_integrals[j, m] * phi_values[k, m] * w_m
        end # m
        # Scale by (2i+1) / 2 (extra 1/2 from mapping integral from [0, 1] -> [-1, 1])
        B[i, j, k] = 0.5f0 * (2 * i + 1) * res
    end # i, j, k

    return SArray{Tuple{n, n, n}, RealT}(B)
end

# Helper function to precompute and store every instance of the interior integral
# when computing the B tensor
#
#   integral_0^zeta phi_j(s) ds
#       = integral_0^(0.5*(xi_m + 1)) phi_j(s) d s
#       = 0.25 * (xi_m + 1) integral_-1^1 phi_j(0.25 (xi_m + 1) (eta_p + 1)) d eta
#
function compute_inner_integrals!(phi_j_integrals, lg_nodes, lg_weights)
    # save preliminary values in pre-allocated memory
    fill!(phi_j_integrals, 0)

    n, M = size(phi_j_integrals)

    for j in 1:n, m in 1:M
        # Outer LG quadrature node and weight
        xi_m = lg_nodes[m]

        # For a given `m` and `j` first compute the inner integral
        # Technically, this integral could use fewer LG nodes, but ignore for now
        for p in 1:M
            # Inner LG quadrature node and weight
            eta_p = lg_nodes[p]
            w_p = lg_weights[p]

            # Compute \phi_j at the (mapped) point 0.25 (xi_m + 1) (eta_p + 1)
            x = 0.25f0 * (xi_m + 1) * (eta_p + 1)
            phi_j_value, _ = shifted_legendre_polynomial_and_derivative(j, x)

            # Add result into the sum
            phi_j_integrals[j, m] += phi_j_value * w_p
        end # p
        # Jacobian from integration moving form [0, (xi_m+1)/2] -> [-1, 1]
        phi_j_integrals[j, m] *= 0.25f0 * (xi_m + 1)
    end # j, m

    return nothing
end

# Create the C matrix for the moment equations if one wants to include friction.
# We use recursion relations for the Legendre polynomial derivatives
# to build the necessary components of the integral
#
#   C_ij = (2i + 1) integral_0^1 phi'_i phi'_j  d zeta
#   i, j = 1, ..., N
#
# where phi'_i and \phi'_j are polynomials of degree at most `N-1`.
# So, the highest degree of the integrand is a polynomial
# of degree `(N-1) + (N-1) = 2N - 2`. We approximate the integrals with
# Legendre-Gauss quadrature which is exact for polynomials up to degree `2M+1`.
# Thus, if we select `M = ceil( (2N - 3) / 2 ) + 1` nodes, then the quadrature
# will be exact for the construction of the C matrix.
#
# It is worth noting that the shifted and unnormalized Legendre polynomials
# are no longer orthogonal on [0, 1], so this matrix is not diagonal.
# If one uses a normalized basis then we get orthogonality back as is a classic
# result for the Legendre polynomials also on [-1, 1], e.g.,
# Abramowitz & Stegun, Handbook of Mathematical Functions, Chapter 22
function compute_C_matrix(n::Int, RealT = Float64)
    # save preliminary values in pre-allocated memory
    C = zero(MArray{Tuple{n, n}, RealT})

    # Given the number of moments `n` determine how many Legendre-Gauss
    # nodes are required and compute them
    M = Int(ceil(n - 1.5f0) + 1)
    lg_nodes, lg_weights = gauss_nodes_weights(M, RealT)

    # Compute everything on the fly (might be slow)
    for i in 1:n, j in 1:n
        # Legendre-Gauss quadrature loop for the outer integral
        res = 0
        for m in 1:M
            # Outer LG quadrature node and weight
            xi_m = lg_nodes[m]
            w_m = lg_weights[m]

            # Compute phi'_i and \phi'_j at the mapped node 0.5 (xi_m + 1)
            x = 0.5f0 * (xi_m + 1)
            _, phi_i_prime = shifted_legendre_polynomial_and_derivative(i, x)
            _, phi_j_prime = shifted_legendre_polynomial_and_derivative(j, x)

            # Add result into the sum
            res += phi_i_prime * phi_j_prime * w_m
        end # m
        # Scale by (2i+1) / 2 (extra 1/2 from mapping integral from [0, 1] -> [-1, 1])
        # TODO: Potentially move the (2i+1) factor out of this routine as it scales
        # all the friction terms
        C[i, j] = 0.5f0 * (2 * i + 1) * res
    end # i, j

    return SArray{Tuple{n, n}, RealT}(C)
end

# The tensor B and matrix C are built from the shifted and unnormalized Legendre polynomials.
# They can be created from the Rodrigues' formula
#
#   phi_j = 1/j! d^j / dx^j (x - x^2)^j
#
# for j = 1, ..., N and x in [0, 1]. A more convenient way to build the polynomials
# is from the three term recurrence relation
#
#   phi_0 = 1,
#   phi_1 = 1 - 2x,
#   phi_j = ((2j-1) (1 - 2x) phi_{j-1} - (j-1) phi_{j-2}) / (j), j = 2, ..., N
#
# The derivatives of the shifted and unnormalized Legendre polynomials can also
# be built from a three term recurrence relation
#
#   phi'_0 = 0,
#   phi'_1 = -2,
#   phi'_j = phi'_{j-2} - 2 (2j-1) phi_{j-1}, j = 2, ..., N
#
# This routine computes and returns phi_N and phi'_N evaluated at a point `x`.
function shifted_legendre_polynomial_and_derivative(N::Int, x::Real)
    if N == 0
        poly = 1
        deriv = 0
    elseif N == 1
        poly = 1 - 2 * x
        deriv = -2
    else
        poly_Nm2 = 1
        poly_Nm1 = 1 - 2 * x
        deriv_Nm2 = 0
        deriv_Nm1 = -2

        poly = 0
        deriv = 0
        for j in 2:N
            poly = ((2 * j - 1) * (1 - 2 * x) * poly_Nm1 - (j - 1) * poly_Nm2) / j
            deriv = deriv_Nm2 - 2 * (2 * j - 1) * poly_Nm1
            poly_Nm2 = poly_Nm1
            poly_Nm1 = poly
            deriv_Nm2 = deriv_Nm1
            deriv_Nm1 = deriv
        end
    end

    return poly, deriv
end

#############################################################################
# Routines taken from Trixi.jl to create the Legendre-Gauss quadrature.
# Can be removed when this functionality is moved into TrixiShallowWater.jl
#############################################################################

"""
    gauss_nodes_weights(n_nodes::Integer, RealT = Float64)

Computes nodes ``x_j`` and weights ``w_j`` for the Gauss-Legendre quadrature.
This implements algorithm 23 "LegendreGaussNodesAndWeights" from the book

- David A. Kopriva, (2009).
  Implementing spectral methods for partial differential equations:
  Algorithms for scientists and engineers.
  [DOI:10.1007/978-90-481-2261-5](https://doi.org/10.1007/978-90-481-2261-5)
"""
function gauss_nodes_weights(n_nodes::Integer, RealT = Float64)
    n_iterations = 20
    tolerance = 2 * eps(RealT) # Relative tolerance for Newton iteration

    # Initialize output
    nodes = ones(RealT, n_nodes)
    weights = zeros(RealT, n_nodes)

    # Get polynomial degree for convenience
    N = n_nodes - 1
    if N == 0
        nodes .= 0
        weights .= 2
        return nodes, weights
    elseif N == 1
        nodes[1] = -sqrt(one(RealT) / 3)
        nodes[end] = -nodes[1]
        weights .= 1
        return nodes, weights
    else # N > 1
        # Use symmetry property of the roots of the Legendre polynomials
        for i in 0:(div(N + 1, 2) - 1)
            # Starting guess for Newton method
            nodes[i + 1] = -cospi(one(RealT) / (2 * N + 2) * (2 * i + 1))

            # Newton iteration to find root of Legendre polynomial (= integration node)
            for k in 0:n_iterations
                poly, deriv = legendre_polynomial_and_derivative(N + 1, nodes[i + 1])
                dx = -poly / deriv
                nodes[i + 1] += dx
                if abs(dx) < tolerance * abs(nodes[i + 1])
                    break
                end

                if k == n_iterations
                    @warn "`gauss_nodes_weights` Newton iteration did not converge"
                end
            end

            # Calculate weight
            poly, deriv = legendre_polynomial_and_derivative(N + 1, nodes[i + 1])
            weights[i + 1] = (2 * N + 3) / ((1 - nodes[i + 1]^2) * deriv^2)

            # Set nodes and weights according to symmetry properties
            nodes[N + 1 - i] = -nodes[i + 1]
            weights[N + 1 - i] = weights[i + 1]
        end

        # If odd number of nodes, set center node to origin (= 0.0) and calculate weight
        if n_nodes % 2 == 1
            poly, deriv = legendre_polynomial_and_derivative(N + 1, zero(RealT))
            nodes[div(N, 2) + 1] = 0
            weights[div(N, 2) + 1] = (2 * N + 3) / deriv^2
        end

        return nodes, weights
    end
end

"""
    legendre_polynomial_and_derivative(N::Int, x::Real)

Computes the Legendre polynomial of degree `N` and its derivative at `x`.
This implements algorithm 22 "LegendrePolynomialAndDerivative" from the book

- David A. Kopriva, (2009).
  Implementing spectral methods for partial differential equations:
  Algorithms for scientists and engineers.
  [DOI:10.1007/978-90-481-2261-5](https://doi.org/10.1007/978-90-481-2261-5)
"""
function legendre_polynomial_and_derivative(N::Int, x::Real)
    RealT = typeof(x)
    if N == 0
        poly = one(RealT)
        deriv = zero(RealT)
    elseif N == 1
        poly = x
        deriv = one(RealT)
    else
        poly_Nm2 = one(RealT)
        poly_Nm1 = copy(x)
        deriv_Nm2 = zero(RealT)
        deriv_Nm1 = one(RealT)

        poly = zero(RealT)
        deriv = zero(RealT)
        for i in 2:N
            poly = ((2 * i - 1) * x * poly_Nm1 - (i - 1) * poly_Nm2) / i
            deriv = deriv_Nm2 + (2 * i - 1) * poly_Nm1
            poly_Nm2 = poly_Nm1
            poly_Nm1 = poly
            deriv_Nm2 = deriv_Nm1
            deriv_Nm1 = deriv
        end
    end

    # Normalize
    poly = poly * sqrt(N + 0.5)
    deriv = deriv * sqrt(N + 0.5)

    return poly, deriv
end

#################################################################################
# Testing the moment matrices
#################################################################################

# TODO: incorporate this into the test section

# To validate the moment matrices verify the identity:
# B_ijk / (2i + 1) + B_kji / (2k + 1) + A_kji / (2k + 1) = 0, for all i, j, k = 1, ..., N
# which follows from integration by parts and the boundary conditions of the Legendre polynomials 
# on the interval [0, 1]. 
# TODO: reference paper

function test_moment_matrices(A, B, n)
    for i in 1:n, j in 1:n, k in 1:n
        res = B[i, j, k] / (2 * i + 1) + B[k, j, i] / (2 * k + 1) + A[k, j, i] / (2 * k + 1)
        #@assert res < 1e-12
        println("i = $i, j = $j, k = $k, res = $res")
    end
    return nothing
end

# Directly compute the A tensor from B using the relation
#  A_ijk / (2i + 1) = - B_ijk / (2i + 1) - B_kji / (2k + 1)
function compute_A_tensor_from_B(B::SArray{Tuple{N, N, N}, RealT}) where {N, RealT}
    A = zero(MArray{Tuple{N, N, N}, RealT})

    for i in 1:N, j in 1:N, k in 1:N
        A[i, j, k] = -(2i + 1) * (B[i, j, k] / (2i + 1) + B[k, j, i] / (2k + 1))
    end # i, j, k

    return SArray{Tuple{N, N, N}, RealT}(A)
end
