# Gaussian precision refits on known graphs, with a max-det reference.
#
#   julia --project=benchmark/refit benchmark/refit/run.jl OUTDIR [--quick] [--worker k/K]
#
# Start Julia with one thread. Completed cases and fitted graphs are reused.

using CovarianceSelection, ProximalBase, ProximalSolvers, CliqueTrees, Krylov
using CliqueTrees.Multifrontal: ChordalCholesky, complete!, ndz, nlz
using LinearAlgebra, Printf, Random, SparseArrays, Statistics

BLAS.set_num_threads(1)

const METHODS = ("ChordalCompletion", "DualCoordinate", "PrimalCoordinate", "RefitADMM", "cliquetrees")
const GRAPHS = ("band", "grid", "er", "glasso-sparse", "glasso-dense")
const COLUMNS = ("case_id", "p", "n", "graph", "edges", "max_degree", "chordal", "seed",
  "quic_lambda", "quic_iterations", "quic_converged", "quic_status", "method", "settings",
  "admm_a", "admm_maxiter", "calibration_case", "calibration_target", "calibration_kkt",
  "calibration_matched", "calibration_status", "repetition", "time_s", "bytes", "iterations",
  "converged", "kkt", "objective", "off_graph_max", "embedding_fill", "status")

clean(x) = replace(string(x), '\t' => ' ', '\n' => ' ', '\r' => ' ')
function failure(err)
  err isa InterruptException && throw(err)
  message = first(split(sprint(showerror, err), '\n'))
  message = replace(message, r"(?:[A-Za-z]:)?/[^\s,;)]+" => "<path>")
  "$(nameof(typeof(err))): $(clean(message))"
end

function atomic_write(f, file, worker)
  mkpath(dirname(file))
  tmp = file * ".$(worker).tmp"
  open(f, tmp, "w")
  mv(tmp, file; force=true)
end

case_id(p, n, graph) = "p$(p)-n$(n)-$(graph)"
case_seed(p, n, graph) = 20261028 + 100000 * p + 10 * n + findfirst(==(graph), GRAPHS)

function sample_covariance(p, n, graph)
  rng = Xoshiro(case_seed(p, n, graph))
  Ω = zeros(p, p)
  for j in 1:p, i in j+1:p
    if rand(rng) < 4 / (p - 1)
      Ω[i, j] = Ω[j, i] = (rand(rng, Bool) ? 1 : -1) * (0.2 + 0.2 * rand(rng))
    end
  end
  for j in 1:p
    Ω[j, j] = sum(abs, view(Ω, :, j)) + 0.2
  end
  Σ = inv(cholesky(Symmetric(Ω)))
  scales = sqrt.(diag(Σ))
  Σ ./= scales * scales'
  X = randn(rng, n, p) * cholesky(Symmetric(Σ)).U
  Matrix(Symmetric(X' * X / n))
end

function fixed_graph(p, graph)
  E = falses(p, p)
  if graph == "band"
    for j in 1:p, i in j+1:min(p, j+3)
      E[i, j] = E[j, i] = true
    end
  elseif graph == "grid"
    rows = isqrt(p)
    while p % rows != 0
      rows -= 1
    end
    cols = p ÷ rows
    for c in 1:cols, r in 1:rows
      j = r + (c - 1) * rows
      if r < rows
        E[j+1, j] = E[j, j+1] = true
      end
      if c < cols
        E[j+rows, j] = E[j, j+rows] = true
      end
    end
  elseif graph == "er"
    rng = Xoshiro(case_seed(p, p, graph) + 1)
    for j in 1:p, i in j+1:p
      E[i, j] = E[j, i] = rand(rng) < 4 / (p - 1)
    end
  else
    error("unknown graph")
  end
  E
end

function graph_data(S, n, graph, outdir, worker)
  p = size(S, 1)
  startswith(graph, "glasso") || return (fixed_graph(p, graph), (NaN, NaN, missing, "not_applicable"))
  file = joinpath(outdir, "graphs", case_id(p, n, graph) * ".tsv")
  if isfile(file)
    lines = readlines(file)
    meta = split(lines[1], '\t')
    E = falses(p, p)
    for line in lines[3:end]
      i, j = parse.(Int, split(line, '\t'))
      E[i, j] = E[j, i] = true
    end
    return E, (parse(Float64, meta[2]), parse(Int, meta[3]), meta[4] == "true", String(meta[5]))
  end
  λ = (graph == "glasso-sparse" ? 1.0 : 0.5) * sqrt(log(p) / n)
  X = Matrix{Float64}(I, p, p)
  offdiag = ones(p, p) - I
  result = solve!(X, LogDetLoss(S), NormL1(λ, offdiag), ProxNewton())
  E = BitMatrix((X .!= 0) .| (X' .!= 0))
  E[diagind(E)] .= false
  meta = (λ, result.passes, result.converged, "ok")
  atomic_write(file, worker) do io
    println(io, join(("# quic", meta...), '\t'))
    println(io, "i\tj")
    for j in 1:p, i in j+1:p
      E[i, j] && println(io, i, '\t', j)
    end
  end
  E, meta
end

function reference_matrix(S, E)
  pattern = copy(E)
  pattern[diagind(pattern)] .= true
  sparse(Symmetric(S .* pattern))
end

function reference_check(S, E, fill_skip_p)
  p = size(S, 1)
  p < fill_skip_p && return NaN, "ok"
  try
    F = ChordalCholesky(reference_matrix(S, E))
    # Match sparse(F.L)'s triangular count without allocating its values.
    entries = (ndz(F) + p) ÷ 2 + nlz(F)
    fraction = entries / (p * (p + 1) ÷ 2)
    status = fraction > 0.5 ? "skipped: dense chordal extension (fill $fraction)" : "ok"
    fraction, status
  catch err
    NaN, failure(err)
  end
end

function reference_fit(S, E)
  A = reference_matrix(S, E)
  F = ChordalCholesky(A)
  complete!(F, Symmetric(A))
  L = sparse(F.L)
  Ω = Matrix(F.P' * (L * L') * F.P)
  off = 0.0
  for j in axes(E, 2), i in axes(E, 1)
    if i != j && !E[i, j]
      off = max(off, abs(Ω[i, j]))
      Ω[i, j] = 0
    end
  end
  Ω, off
end

function fit(S, E, method, a, m)
  if method == "cliquetrees"
    Ω, off = reference_fit(S, E)
    return Ω, nothing, off
  end
  solver = if method == "ChordalCompletion"
    ChordalCompletion()
  elseif method == "DualCoordinate"
    DualCoordinate()
  elseif method == "PrimalCoordinate"
    PrimalCoordinate()
  else
    RefitADMM(ADMMOptions(; abstol=a, reltol=a, maxiter=m))
  end
  info = RefitInfo()
  Ω = covsel_refit(S, E, solver; info)
  Ω, info, 0.0
end

function diagnostics(Ω, S, E, need_kkt)
  Ω = Matrix(Ω)
  F = cholesky(Symmetric(Ω); check=false)
  isposdef(F) || return NaN, Inf
  objective = dot(S, Ω) - logdet(F)
  need_kkt || return objective, NaN
  W = inv(F)
  kkt = 0.0
  for j in axes(E, 2), i in j:size(E, 1)
    if i == j || E[i, j]
      kkt = max(kkt, abs(W[i, j] - S[i, j]) / sqrt(S[i, i] * S[j, j]))
    end
  end
  objective, kkt
end

function calibrate(S, E)
  target, reference_status = 1e-9, "ok"
  try
    _, info, _ = fit(S, E, "PrimalCoordinate", NaN, 0)
    if isfinite(info.kkt)
      target = max(info.kkt, 1e-9)
    else
      reference_status = "nonfinite_primal_kkt"
    end
  catch err
    reference_status = failure(err)
  end
  kkt, a, status = Inf, 1e-10, reference_status
  for tolerance in (1e-6, 1e-7, 1e-8, 1e-9, 1e-10)
    a = tolerance
    try
      _, info, _ = fit(S, E, "RefitADMM", a, 20000)
      kkt = info.kkt
      status = reference_status
      if kkt <= target && reference_status == "ok"
        return (a=a, target=target, kkt=kkt, matched=true, status=status)
      end
    catch err
      kkt, status = Inf, failure(err)
    end
  end
  (a=a, target=target, kkt=kkt, matched=false, status=status == "ok" ? "target_not_reached" : status)
end

function calibration(S, E, p, n, graph, outdir, worker, quick)
  quick && return (a=1e-8, target=NaN, kkt=NaN, matched=missing, status="quick_uncalibrated", source="none")
  q = min(p, 500)
  nq = n < p ? q ÷ 2 : 2q
  source = case_id(q, nq, graph)
  file = joinpath(outdir, "calibration", source * ".tsv")
  if isfile(file)
    values = split(readlines(file)[2], '\t')
    return (a=parse(Float64, values[1]), target=parse(Float64, values[2]),
      kkt=parse(Float64, values[3]), matched=values[4] == "true", status=String(values[5]), source=source)
  end
  if q == p
    result = calibrate(S, E)
  else
    Sq = sample_covariance(q, nq, graph)
    Eq, _ = graph_data(Sq, nq, graph, outdir, worker)
    result = calibrate(Sq, Eq)
  end
  atomic_write(file, worker) do io
    println(io, "a\ttarget\tkkt\tmatched\tstatus")
    println(io, join(clean.(Tuple(result)), '\t'))
  end
  merge(result, (source=source,))
end

function timed_fit(S, E, method, a, m)
  GC.gc()
  timed = @timed try
    (value=fit(S, E, method, a, m), status="ok")
  catch err
    (value=nothing, status=failure(err))
  end
  result = timed.value
  if result.status != "ok"
    return (seconds=timed.time, bytes=timed.bytes, iterations=NaN, converged=missing,
      kkt=NaN, objective=NaN, off=NaN, status=result.status)
  end
  Ω, info, off = result.value
  objective, kkt = diagnostics(Ω, S, E, info === nothing || method == "ChordalCompletion")
  (seconds=timed.time, bytes=timed.bytes, iterations=info === nothing ? NaN : info.iterations,
    converged=info === nothing ? missing : info.converged, kkt=info === nothing || method == "ChordalCompletion" ? kkt : info.kkt,
    objective=objective, off=off, status="ok")
end

function run_case(p, n, graph, outdir, worker, quick, fill_skip_p)
  id = case_id(p, n, graph)
  file = joinpath(outdir, "cases", id * ".tsv")
  if isfile(file)
    println("skip ", id)
    return
  end
  println("run ", id)
  flush(stdout)
  S = sample_covariance(p, n, graph)
  E, quic = graph_data(S, n, graph, outdir, worker)
  chordal = is_chordal(E)
  degree = maximum(vec(sum(E; dims=1)))
  cal = calibration(S, E, p, n, graph, outdir, worker, quick)
  m = quick ? 1000 : p <= 500 ? 20000 : p == 1000 ? 600 : 300
  embedding_fill, reference_status = reference_check(S, E, fill_skip_p)
  prefix = (id, p, n, graph, count(E) ÷ 2, degree, chordal, case_seed(p, n, graph), quic...)
  rows = Vector{Any}()
  counts, targets = zeros(Int, 5), ones(Int, 5)
  elapsed = zeros(5)
  first_times = zeros(5)
  r = 1
  while any(counts .< targets)
    for offset in 0:4
      j = mod(r - 1 + offset, 5) + 1
      counts[j] >= targets[j] && continue
      method = METHODS[j]
      if j == 1 && !chordal
        targets[j] = 0
        continue
      end
      result = if j == 5 && reference_status != "ok"
        (seconds=NaN, bytes=NaN, iterations=NaN, converged=missing,
          kkt=NaN, objective=NaN, off=NaN, status=reference_status)
      else
        timed_fit(S, E, method, cal.a, m)
      end
      counts[j] += 1
      elapsed[j] += result.seconds
      settings = j == 4 ? "abstol=$(cal.a);reltol=$(cal.a);maxiter=$m" :
        j in (2, 3) ? "tol=1e-10;maxiter=1000" : j == 5 ? "complete! defaults;project=true" : "none"
      push!(rows, (prefix..., method, settings, j == 4 ? cal.a : NaN, j == 4 ? m : NaN,
        cal.source, cal.target, cal.kkt, cal.matched, cal.status, counts[j], result.seconds,
        result.bytes, result.iterations, result.converged, result.kkt, result.objective, result.off,
        j == 5 ? embedding_fill : NaN, result.status))
      if counts[j] == 1
        t = first_times[j] = result.seconds
        targets[j] = quick || (j == 4 && p >= 1000) || result.status != "ok" ? 1 : t < 0.1 ? 30 : t < 1 ? 10 :
          t < 60 ? max(3, floor(Int, 30 / t)) : 1
      elseif result.status != "ok"
        targets[j] = counts[j]
      end
      if 1 <= first_times[j] < 60 && counts[j] >= 3
        expected = elapsed[j] / counts[j]
        elapsed[j] + expected > 30 && (targets[j] = counts[j])
      end
    end
    r += 1
  end
  atomic_write(file, worker) do io
    println(io, join(COLUMNS, '\t'))
    for row in rows
      println(io, join(clean.(row), '\t'))
    end
  end
end

function metadata(outdir, worker, workers, explicit_worker, quick, fill_skip_p)
  commit = try
    readchomp(pipeline(`git -C $(joinpath(@__DIR__, "../..")) rev-parse HEAD`; stderr=devnull))
  catch
    "unavailable"
  end
  utc = try
    readchomp(`date -u +%Y-%m-%dT%H:%M:%SZ`)
  catch
    "unix_seconds=$(time())"
  end
  data = [("date_utc", utc), ("julia", string(VERSION)), ("package_commit", commit),
    ("cpu_model", Sys.cpu_info()[1].model), ("blas_vendor", string(BLAS.vendor())),
    ("blas_config", sprint(show, BLAS.get_config())), ("julia_threads", string(Threads.nthreads())),
    ("blas_threads", string(BLAS.get_num_threads())), ("worker", "$worker/$workers"), ("quick", string(quick)),
    ("fill_skip_p", string(fill_skip_p))]
  for package in (CovarianceSelection, ProximalBase, ProximalSolvers, CliqueTrees, Krylov)
    push!(data, (string(nameof(package)), string(pkgversion(package))))
  end
  file = joinpath(outdir, explicit_worker ? "meta-$worker.tsv" : "meta.tsv")
  isfile(file) && return
  atomic_write(file, worker) do io
    println(io, "key\tvalue")
    for (key, value) in data
      println(io, key, '\t', clean(value))
    end
  end
end

function main(args)
  isempty(args) && error("usage: run.jl OUTDIR [--quick] [--worker k/K]")
  outdir = args[1]
  quick, explicit_worker = false, false
  worker, workers = 1, 1
  i = 2
  while i <= length(args)
    if args[i] == "--quick"
      quick = true
    elseif args[i] == "--worker" && i < length(args)
      worker, workers = parse.(Int, split(args[i+1], '/'))
      explicit_worker = true
      i += 1
    else
      error("unknown or incomplete option")
    end
    i += 1
  end
  1 <= worker <= workers || error("worker must satisfy 1 <= k <= K")
  Threads.nthreads() == 1 || error("start Julia with one thread")
  fill_skip_p = parse(Int, get(ENV, "REFIT_FILL_SKIP_P", "1000"))
  fill_skip_p >= 1 || error("REFIT_FILL_SKIP_P must be positive")
  mkpath(outdir)
  metadata(outdir, worker, workers, explicit_worker, quick, fill_skip_p)
  dimensions = quick ? (30, 60) : (100, 500, 1000, 2000)
  graphs = quick ? GRAPHS[1:3] : GRAPHS
  cases = [(p, n, g) for p in dimensions for n in (quick ? (2p,) : (p ÷ 2, 2p)) for g in graphs]
  sort!(cases; by=c -> (-c[1], c[2], findfirst(==(c[3]), GRAPHS)))
  assigned = cases[worker:workers:end]
  pending = any(!isfile(joinpath(outdir, "cases", case_id(c...) * ".tsv")) for c in assigned)
  if pending
    S = sample_covariance(8, 16, "band")
    E = fixed_graph(8, "band")
    for method in METHODS
      try
        fit(S, E, method, 1e-8, 1000)
      catch err
        println("warmup ", method, ": ", failure(err))
      end
    end
  end
  for (p, n, graph) in assigned
    run_case(p, n, graph, outdir, worker, quick, fill_skip_p)
  end
end

main(ARGS)
