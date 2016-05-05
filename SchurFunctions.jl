

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
    n < 0 && throw(ArgumentError("n must be ≥ 0"))
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

   [i; int2partition(k - pnums[n, i-1], n - i)]  # concatenates i onto the front
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
A2PStored = Matrix{Int}[]
for i = 1:max_n
    num_partitions = pnums[i, i]
    push!(A2PStored, zeros(Int, num_partitions, num_partitions))
    A2PStored[i][end, end] = 1
end

doc"""
Compute the augmented monomial symmetric function indexed by
the partition P as an array of coefficients in the p basis (power-sum basis)

We use the fact that
`am[l1, l2, ..., lk] = p[lk]*am[l1, l2, ..., lk]
    - sum(i=1..k-1)am[l1, l2, ..., l(i-1), li+lk, l(i+1), ..., lk]`

The results are stored in the table A2PStored, to reduce redundancy.

A vector of zeros is made, and then each term is added in as it is computed.
"""

function augmon2psum(P::Vector{Int})
   psize = sum(P)
   accum = zeros(Int, pnums[psize, psize])  # accumulator

   k = partition2int(P)
   A2PStored[psize][k, k] == 1 && return A2PStored[psize][:, k]

   if length(P) == 1
       accum[end] = 1
       return accum
   end


   firstterm = augmon2psum(P[1:end-1])
   # firsterm is of form p_(something) * m̃_(something)

   for i = 1:length(firstterm)  # iterating over the thing that was a partition of the smaller object
       new_partition = [P[end]; int2partition(i, psize  -P[end])]
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

   A2PStored[psize][:,k] = accum

   accum
end


function weighted_dot_product(v1, v2, weight)

    total = zero(promote(v1[1], weight[1])[1])

    for i in 1:length(v1)
        total += v1[i]*v2[i]*weight[i]
    end

    #@show total
    total
end

function character_table(n)
    macz = Array(Int, pnums[n, n])  # Macdonald's z function

    for i = 1:pnums[n,n]
        P = int2partition(i, n)
        count = zeros(Int, n)

        for j in P
            count[j] += 1 # turns into representation of partition as 1^4.2^6 etc.
        end
        macz[i] = prod(map(factorial, count)) * prod(P)
    end
    # maczM = diagm(macz)



    # A side-effect of computing
    #     augmon2psum(ones(Int,n))
    # is that it populates the entire A2PStored matrix

    @time augmon2psum(ones(Int, n))

    #@show A2PStored[n]


    Schur = convert(Array{Rational{Int128}}, A2PStored[n])
    macz_weight = macz # convert(Array{Rational{Int128}}, macz)

    #@show Schur
    #@show macz_weight


    for i = 1:pnums[n, n]
        # Orthogonalize with respect to all earlier vectors using Gram-Schmidt:
        for j = 1:i-1

            #@show i, j, weighted_dot_product(Schur[:, i], Schur[:, j], macz_weight)

            Schur[:, i] -= weighted_dot_product(Schur[:, i], Schur[:, j], macz_weight) * Schur[:,j]
        	# Schur[:, i] -= ( Schur[:, i]' * maczM * Schur[:,j] )[1] * Schur[:,j]
        end

        # Normalize
        # norm = round(Int, sqrt(Schur[:, i]' * maczM * Schur[:, i]))

        #@show Schur[:, i]
        dot = weighted_dot_product(Schur[:, i], Schur[:, i], macz_weight)
        #@show dot
        norm = trunc(Int, sqrt(dot))  # we know that it's actually an integer

        Schur[:,i] /= norm
    end

    #@show Schur

    # normalize rows by first element to get character table:

    χ = copy(Schur)
    for i in 1:size(χ, 1)
        first = χ[i, 1]
        χ[i,:] /= first
    end
    #diagm(vec(Schur[:,1]) .^ (-1)) * Schur

    χ
end
