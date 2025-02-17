using CUDA
CUDA.functional()  # Returns true if CUDA is working

A = CUDA.rand(1000, 1000)  # Create a large matrix on GPU
B = CUDA.rand(1000, 1000)

C = A * B  # GPU-accelerated matrix multiplication

cpu_array = rand(1000, 1000)  # Create array on CPU
gpu_array = CuArray(cpu_array)  # Move to GPU

result = Array(gpu_array)  # Move back to CPU

function gpu_add_kernel(a, b, c)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(a)
        c[i] = a[i] + b[i]
    end
    return
end

# Allocate GPU memory
N = 100_000
a = CUDA.rand(N)
b = CUDA.rand(N)
c = CUDA.zeros(N)

# Launch the kernel with threads and blocks
@cuda threads=256 blocks=ceil(Int, N/256) gpu_add_kernel(a, b, c)

CUDA.@sync C = A * B

function shared_memory_example(A, B, C)
    shared_A = CuDynamicSharedArray(Float32, (blockDim().x,))
    shared_B = CuDynamicSharedArray(Float32, (blockDim().x,))
    
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(A)
        shared_A[threadIdx().x] = A[i]
        shared_B[threadIdx().x] = B[i]
        sync_threads()  # Ensure all threads sync before using shared memory
        C[i] = shared_A[threadIdx().x] + shared_B[threadIdx().x]
    end
    return
end
