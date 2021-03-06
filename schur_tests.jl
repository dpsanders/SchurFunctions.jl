
using FactCheck
include("SchurFunctions.jl")

facts("Enumerating partitions") do
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

    @fact p(30) --> 5604
    @fact p(30,5) - p(30,4) --> 377  # exactly 5 parts
    @fact p(100) --> 190569292
end

facts("Character tables") do
    @fact character_table(3) --> [  1  -1   1
                                    2   0  -1
                                    1   1   1  ]

    @fact character_table(6) --> [1  -1   1  -1   1  -1   1  -1   1   1  -1
  5  -3   1   1   2   0  -1  -1  -1   0   1
  9  -3   1  -3   0   0   0   1   1  -1   0
  5  -1   1   3  -1  -1   2   1  -1   0   0
 10  -2  -2   2   1   1   1   0   0   0  -1
 16   0   0   0  -2   0  -2   0   0   1   0
  5   1   1  -3  -1   1   2  -1  -1   0   0
 10   2  -2  -2   1  -1   1   0   0   0   1
  9   3   1   3   0   0   0  -1   1  -1   0
  5   3   1  -1   2   0  -1   1  -1   0  -1
  1   1   1   1   1   1   1   1   1   1   1]


end
