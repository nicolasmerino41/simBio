npp_DA_relative_to_1000 = deepcopy(npp_DA)
max_npp = maximum(npp_DA_relative_to_1000[.!isnan.(npp_DA_relative_to_1000)])

for i in idx
    npp_DA_relative_to_1000[i] = npp_DA[i] / max_npp * 1000
end

# serialize("Objects\\npp_DA_relative_to_1000.jls", npp_DA_relative_to_1000)