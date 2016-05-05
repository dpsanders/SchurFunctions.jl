using FactCheck
include("SchurFunctions.jl")

facts("Partition map thing") do
    n = 5
    partition_list = [int2partition(k, n) for k in 1:p(n)]
    @fact partition_list --> Vector{Int}[
         [1,1,1,1,1],
         [2,1,1,1],
         [2,2,1],
         [3,1,1],
         [3,2],
         [4,1],
         [5],
         ]

    @fact map(partition2int, partition_list) --> collect(1:7)

end
