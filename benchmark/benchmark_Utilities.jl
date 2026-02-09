using VegSci
using BenchmarkTools
using DataFrames

bench_df = DataFrame(Dict("Rows" => [],
                          "Columns" => [],
                          "Size" => [],
                          "Time" => [],
                          "Memory" => [],
                          "Allocations" => []))

for i in collect(1:7)

    rown = 10
    coln = 10^i

    btime_i = @benchmark VegSci.generate_test_array(rown = rown, coln = coln, 
                                                    meancoloccs = 10, 
                                                    rowprefix = "Species", colprefix = "Releve", 
                                                    sparse_array = true)

    mean_bench = mean(btime_i)

    mean_time = mean_bench.time
    mean_memory = mean_bench.memory
    allocations = mean_bench.allocs

    bench_df_i = DataFrame(Dict("Rows" => rown,
                                "Columns" => coln,
                                "Size" => rown * coln,
                                "Time" => mean_time,
                                "Memory" => mean_memory,
                                "Allocations" => allocations))

    bench_df = vcat(bench_df, bench_df_i)

end