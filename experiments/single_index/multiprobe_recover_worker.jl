include("multiprobe_worker.jl")
for line in readlines(ARGS[2])
    isempty(strip(line)) && continue
    mp_worker(ARGS[1],parse(Int,line),ARGS[3])
end
