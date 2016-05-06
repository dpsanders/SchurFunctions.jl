

# All the interesting symmetric function bases are
# naturally indexed by partitions, but we'd like to
# store them as arrays of coefficients.  We could use
# indexing by arbitrary sets, but instead we'll index
# by integers, using an explicit map between partitions
# and integers.

## Number of integer partitions
doc"""$p(n,k)$ is number of integer partitions of $n$ with largest part at most $k$
"""
function p(n, k)
     (k == 1 || n == 0) && return 1
     if n >= k
         p(n, k - 1) + p(n - k, k)
     else
         p(n, n)
     end
end


# Should memoize the result of calculating the p's
# and populate pnums more carefully to avoid redundant evaluations of p
const max_n = 25
const pnums = [p(i, j) for i = 1:max_n, j = 1:max_n]

p(n) = n < max_n ? pnums[n, n] : p(n, n)

doc"""
    int2partition(k, n)

Return the $k$th partition in the list of partitions of $n$.

Partitions are enumerated in lexicographic order, e.g.
for $n=5$ there are $p(5)=7$ partitions:

- 11111
- 2111
- 221
- 331
- 32
- 41
- 5
"""

function int2partition(k::Integer, n::Integer)
    # n < 0 && throw(ArgumentError("n must be ≥ 0"))

    n == 0 && return Array(Int, 0)

   if k > pnums[n, n]  # If k too big, take it modulo pnums[n,n]
       k = mod1(k + pnums[n, n] - 1, pnums[n, n])
   end

   # Determine the largest part of P
   i = 1
   while pnums[n, i] < k
       i += 1
   end

   i == 1 && return ones(Int, n)

   unshift!(int2partition(k - pnums[n, i-1], n - i), i )   # concatenates i onto the front
end

doc"""`partition2int` is the inverse function of `int2partition`:
given a partition $P$, return $k$ such that `int2partition(k, sum(P)) == P`
"""
function partition2int(P::Vector{Int})
   tn = sum(P)
   k = 1

   for i in P
       i < 2 && break

       k += pnums[tn, i-1]
       tn -= i
   end

   k
end


"Augmented monomial to power-sum"
const augmented_to_power_sum = Matrix{Int}[]

function setup_augmented_to_power_sum()
    for i = 1:max_n
        num_partitions = pnums[i, i]
        push!(augmented_to_power_sum, zeros(Int, num_partitions, num_partitions))
        augmented_to_power_sum[i][end, end] = 1
    end
end

setup_augmented_to_power_sum()


doc"""
Compute the augmented monomial symmetric function indexed by
the partition P as an array of coefficients in the p basis (power-sum basis)

We use the fact that
`am[l1, l2, ..., lk] = p[lk]*am[l1, l2, ..., lk]
    - sum(i=1..k-1) am[l1, l2, ..., l(i-1), li+lk, l(i+1), ..., lk]`

The results are stored in the table augmented_to_power_sum, to reduce redundancy.

A vector of zeros is made, and then each term is added in as it is computed.
"""

function augmon2psum(P::Vector{Int})
   psize = sum(P)
   accum = zeros(Int, pnums[psize, psize])  # accumulator

   k = partition2int(P)
   augmented_to_power_sum[psize][k, k] == 1 && return augmented_to_power_sum[psize][:, k]

   if length(P) == 1
       accum[end] = 1
       return accum
   end


   firstterm = augmon2psum(P[1:end-1])
   # firsterm is of form p_(something) * am(something)

   for i = 1:length(firstterm)  # iterating over the thing that was a partition of the smaller object
       new_partition = unshift!(int2partition(i, psize - P[end]), P[end])
       sort!(new_partition, rev=true)

       accum[ partition2int(new_partition) ] += firstterm[i]
   end


   tosum = zeros(Int, pnums[psize, psize])
   # use to collect similar terms like m[3,1,1,1] vs m[1,3,1,1]

   for i = 1:length(P) - 1
       P[i] += P[end]
       tosum[partition2int(sort(P[1:end-1],rev=true))] -= 1
       P[i] -= P[end]
   end

   for i=1:length(tosum)
       if tosum[i] != 0
           accum += tosum[i] * augmon2psum(int2partition(i, psize))
       end
   end

   augmented_to_power_sum[psize][:,k] = accum

   accum
end



doc"""Macdonald $z$-function weight"""
function macz_weight(n)
    macz = zeros(Int128, pnums[n, n])

    for i = 1:pnums[n,n]
        P = int2partition(i, n)

        powers = zeros(Int, n)
        for j in P
            powers[j] += 1  # finds representation of partition using powers as  1^4 2^6 etc.
        end

        product = Int128(1)
        for k in powers
            product *= factorial(k)
        end

        for k in P
            product *= k
        end

        #macz[i] = prod(map(factorial, powers)) * prod(P)
        #@show macz[i] - product
        macz[i] = product
    end

    macz
end

doc"""Calculate weighted dot product of `v1` and `v2` with given `weight` vector."""
function weighted_dot_product{T}(v1::Vector{T}, v2::Vector{T}, weight::Vector)

    total = zero(T)

    for i in 1:length(v1)
        total += v1[i] * v2[i] * weight[i]
    end

    total
end


doc"""Calculate the coefficients of the Schur function of size n in terms of the augmented monomial basis"""
function Schur_function(n)

    macz = macz_weight(n)

    @time augmon2psum(ones(Int, n)) # Side-effect is populating the entire augmented_to_power_sum matrix

    Schur = convert(Array{Rational{Int128}}, augmented_to_power_sum[n])

    # Orthogonalize with respect to all earlier vectors using Gram-Schmidt:

    for i = 1:pnums[n, n]
        for j = 1:i-1

            Schur[:, i] -= weighted_dot_product(Schur[:, i], Schur[:, j], macz) * Schur[:, j]

        end

        # normalize:
        normsq = weighted_dot_product(Schur[:, i], Schur[:, i], macz)
        norm = round(Int, sqrt(normsq))  # we know that normsq must actually be an integer

        Schur[:, i] /= norm

    end

    Schur
end

doc"""Calculate the character table for the group S_n"""
function character_table(n)

    χ = map(Int, inv(Schur_function(n)))

end
