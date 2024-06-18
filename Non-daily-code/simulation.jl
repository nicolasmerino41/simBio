carrying_capacity = rand(1000, 1000)
rate = rand(1000, 1000, 80)
abundance = rand(1000, 1000)
growth(abundance, rate, carrying_capacity) = abundance * rate * (1 - (abundance / carrying_capacity))

using ThreadsX # use ThreadsX.map! instead of map! to run threaded

function run(abundance, rate, carrying_capacity)
    n_timesteps = size(rate, 3)
    for t in 1:n_timesteps
        rate_t = view(rate, :, :, t)
        map!(abundance, abundance, rate_t, carrying_capacity) do a, r, c
            a + growth(a, r, c)
        end
    end
end

@time run(abundance, rate, carrying_capacity) # first run is slower because of compilation
run(abundance, rate, carrying_capacity)
@profview run(abundance, rate, carrying_capacity)



