# Summarize timings, objective agreement, hypotheses, and failures.
#
#   julia --project=benchmark/refit benchmark/refit/summarize.jl OUTDIR

using Printf, Statistics

const METHODS = ("ChordalCompletion", "DualCoordinate", "PrimalCoordinate", "RefitADMM", "cliquetrees")
number(row, key) = parse(Float64, row[key])
ok(row) = row["status"] == "ok"
skipped(row) = startswith(row["status"], "skipped:")
failed(row) = !ok(row) && !skipped(row)
finite_objective(row) = ok(row) && isfinite(number(row, "objective"))
exception_type(row) = first(split(row["status"], ':'; limit=2))
fmt(x) = @sprintf("%.6g", x)
md(x) = replace(string(x), '|' => "\\|")

function read_rows(outdir)
  directory = joinpath(outdir, "cases")
  isdir(directory) || error("no cases directory")
  rows = Dict{String,String}[]
  for file in sort(filter(f -> endswith(f, ".tsv"), readdir(directory; join=true)))
    lines = readlines(file)
    isempty(lines) && error("empty case file")
    header = split(lines[1], '\t')
    for line in lines[2:end]
      fields = split(line, '\t'; keepempty=true)
      length(fields) == length(header) || error("malformed case row")
      push!(rows, Dict(String(k) => String(v) for (k, v) in zip(header, fields)))
    end
  end
  isempty(rows) && error("no completed cases")
  rows
end

method_rows(rows, method) = filter(r -> r["method"] == method, rows)
function method_skipped(rows, method)
  subset = method_rows(rows, method)
  !isempty(subset) && all(skipped, subset)
end

lower_bound(rows, method) = any(r -> ok(r) && r["converged"] == "false", method_rows(rows, method))

function median_time(rows, method)
  subset = filter(ok, method_rows(rows, method))
  isempty(subset) ? nothing : median(number(r, "time_s") for r in subset)
end

function timing_cell(rows, method)
  subset = method_rows(rows, method)
  isempty(subset) && return "—"
  successes = filter(ok, subset)
  failures = filter(failed, subset)
  skips = filter(skipped, subset)
  if isempty(successes) && isempty(failures) && !isempty(skips)
    return "skipped (fill $(fmt(number(skips[1], "embedding_fill"))))"
  end
  errors = join(sort(unique(exception_type.(failures))), ", ")
  isempty(successes) && return "fails: $errors"
  times = [number(r, "time_s") for r in successes]
  iterations = [number(r, "iterations") for r in successes]
  k = all(isfinite, iterations) ? fmt(median(iterations)) : "NaN"
  cell = "$(fmt(median(times))) ($(fmt(minimum(times)))) / $k"
  lower_bound(rows, method) && (cell *= " [lower bound]")
  isempty(failures) ? cell : cell * "; fails: $errors"
end

function hypothesis(io, label, description, cases, eligible, left, right; fastest=false)
  deciding, counterexamples, unavailable, bounds, omitted = String[], String[], String[], String[], String[]
  for id in sort(collect(keys(cases)))
    rows = cases[id]
    eligible(rows) || continue
    opponents = fastest ? [m for m in METHODS if m != left && !isempty(method_rows(rows, m))] : [right]
    filter!(m -> !method_skipped(rows, m), opponents)
    if method_skipped(rows, left) || isempty(opponents)
      push!(omitted, id)
      continue
    end
    a = median_time(rows, left)
    if a === nothing
      push!(unavailable, id)
      continue
    end
    capped_subject = lower_bound(rows, left)
    outcomes = Symbol[]
    for opponent in opponents
      b = median_time(rows, opponent)
      outcome = if b === nothing
        :unavailable
      elseif capped_subject && lower_bound(rows, opponent)
        :lower_bound
      elseif capped_subject
        a >= b ? :counterexample : :lower_bound
      elseif lower_bound(rows, opponent)
        a < b ? :support : :lower_bound
      else
        b <= a ? :counterexample : :support
      end
      push!(outcomes, outcome)
    end
    if :counterexample in outcomes
      push!(counterexamples, id)
    elseif :lower_bound in outcomes
      push!(bounds, id)
    elseif :unavailable in outcomes
      push!(unavailable, id)
    else
      push!(deciding, id)
    end
  end
  verdict = !isempty(counterexamples) ? "refuted" :
    isempty(deciding) || !isempty(unavailable) || !isempty(bounds) ? "not evaluated completely" : "confirmed"
  println(io, "- **$label: $verdict.** $description")
  !isempty(deciding) && println(io, "  Supporting cases: ", join(deciding, ", "), ".")
  !isempty(counterexamples) && println(io, "  Counterexamples (including ties against uncapped opponents): ", join(counterexamples, ", "), ".")
  !isempty(unavailable) && println(io, "  Unavailable comparisons: ", join(unavailable, ", "), ".")
  println(io, "  Unavailable (lower bound): $(length(bounds)) cases", isempty(bounds) ? "." : ": $(join(bounds, ", ")).")
  !isempty(omitted) && println(io, "  Excluded by skips: ", join(omitted, ", "), ".")
  isempty(deciding) && isempty(counterexamples) && isempty(unavailable) && isempty(bounds) && isempty(omitted) && println(io, "  No eligible cases.")
end

function summary(rows)
  cases = Dict{String,Vector{Dict{String,String}}}()
  for row in rows
    push!(get!(cases, row["case_id"], Dict{String,String}[]), row)
  end
  io = IOBuffer()
  println(io, "# Precision refit benchmark\n")
  println(io, "Completed cases: $(length(cases)); repetition rows: $(length(rows)). Conclusions apply to these cases only.\n")
  println(io, "Entries are median s (minimum s) / median iterations. NaN denotes unavailable iterations; — denotes an inapplicable method.\n")
  for p in sort(unique(parse(Int, r["p"]) for r in rows))
    println(io, "## p = $p\n")
    println(io, "| graph, n | ", join(METHODS, " | "), " |")
    println(io, "|---|", join(fill("---:", 5), '|'), "|")
    ids = sort([id for (id, rs) in cases if parse(Int, rs[1]["p"]) == p];
      by=id -> (cases[id][1]["graph"], parse(Int, cases[id][1]["n"])))
    for id in ids
      rs = cases[id]
      println(io, "| $(rs[1]["graph"]), $(rs[1]["n"]) | ",
        join([md(timing_cell(rs, m)) for m in METHODS], " | "), " |")
    end
    println(io)
  end
  println(io, "## Accuracy gate\n")
  println(io, "Qualifying fits have status ok and kkt ≤ 1e-6. Objective differences use unnormalized tr(SΩ) − logdet(Ω).\n")
  println(io, "| case | qualifying methods | largest pairwise difference | ≤ 1e-8 |")
  println(io, "|---|---:|---:|---|")
  agreement, convergence = true, true
  violations, no_pairs = String[], String[]
  for id in sort(collect(keys(cases)))
    rs = cases[id]
    qualified = filter(r -> finite_objective(r) && number(r, "kkt") <= 1e-6, rs)
    methods = unique(r["method"] for r in qualified)
    difference = 0.0
    for (i, a) in enumerate(methods), b in methods[i+1:end]
      oa = [number(r, "objective") for r in qualified if r["method"] == a]
      ob = [number(r, "objective") for r in qualified if r["method"] == b]
      difference = max(difference, abs(maximum(oa) - minimum(ob)), abs(maximum(ob) - minimum(oa)))
    end
    pass = difference <= 1e-8
    agreement &= pass
    length(methods) < 2 && push!(no_pairs, id)
    println(io, "| $id | $(length(methods)) | $(fmt(difference)) | $(pass ? "pass" : "FAIL") |")
    finite = filter(finite_objective, rs)
    best = isempty(finite) ? Inf : minimum(number(r, "objective") for r in finite)
    for r in rs
      if ok(r) && r["converged"] == "true"
        objective = number(r, "objective")
        if !isfinite(objective) || objective - best > 1e-8
          convergence = false
          push!(violations, "$id/$(r["method"])/rep$(r["repetition"]) (excess $(fmt(objective - best)))")
        end
      end
    end
  end
  println(io, "\nPairwise gate: **$(agreement ? "PASS" : "FAIL")**. Converged-objective gate: **$(convergence ? "PASS" : "FAIL")**. Overall: **$(agreement && convergence ? "PASS" : "FAIL")**.\n")
  !isempty(no_pairs) && println(io, "No qualifying method pair (difference defined as zero): ", join(no_pairs, ", "), ".\n")
  for violation in violations
    println(io, "- Converged-objective violation: $violation")
  end
  println(io, "\n## Hypotheses\n")
  println(io, "Strict comparisons use median times of successful repetitions. Skipped opponents are excluded. Nonconverged iterative times are lower bounds. With subject time a and opponent time b: a capped opponent supports only when a < b and cannot refute; a capped subject refutes only when a ≥ b against an uncapped opponent. Other capped comparisons, including both capped, are unavailable (lower bound). With neither capped, b ≤ a refutes and a < b supports. Missing comparisons cannot confirm a hypothesis.\n")
  hypothesis(io, "H1", "PrimalCoordinate is faster than RefitADMM at maximum degree ≤ 50.", cases,
    rs -> parse(Int, rs[1]["max_degree"]) <= 50, "PrimalCoordinate", "RefitADMM")
  hypothesis(io, "H2", "DualCoordinate is faster than PrimalCoordinate wherever DualCoordinate succeeds, at maximum degree ≤ 10 and p ≥ 500.", cases,
    rs -> parse(Int, rs[1]["max_degree"]) <= 10 && parse(Int, rs[1]["p"]) >= 500 && median_time(rs, "DualCoordinate") !== nothing,
    "DualCoordinate", "PrimalCoordinate")
  hypothesis(io, "H3", "ChordalCompletion is fastest on band cases.", cases,
    rs -> rs[1]["graph"] == "band", "ChordalCompletion", ""; fastest=true)
  println(io, "\n## Failures and diagnostics\n")
  for method in METHODS
    failures = filter(r -> failed(r) && r["method"] == method, rows)
    ids = sort(unique(r["case_id"] for r in failures))
    types = sort(unique(exception_type.(failures)))
    println(io, "- $method: $(length(failures)) failed repetitions in $(length(ids)) cases",
      isempty(ids) ? "." : "; cases: $(join(ids, ", ")); types: $(join(types, ", ")).")
    skips = filter(r -> skipped(r) && r["method"] == method, rows)
    details = sort(unique("$(r["case_id"]) (fill $(fmt(number(r, "embedding_fill"))))" for r in skips))
    println(io, "- $method: $(length(skips)) skipped rows",
      isempty(details) ? "." : "; cases: $(join(details, ", ")).")
  end
  outside = [number(r, "off_graph_max") for r in rows if r["method"] == "cliquetrees" && ok(r)]
  filter!(isfinite, outside)
  println(io, "- cliquetrees largest |Ω| off E before projection: ", isempty(outside) ? "unavailable" : fmt(maximum(outside)), ".")
  nonconverged = sort(unique("$(r["case_id"])/$(r["method"])" for r in rows if ok(r) && r["converged"] == "false"))
  println(io, "- Nonconverged successful returns: ", isempty(nonconverged) ? "none" : join(nonconverged, ", "), ".")
  quic = sort(unique(r["case_id"] for r in rows if startswith(r["graph"], "glasso") && r["quic_converged"] != "true"))
  println(io, "- Nonconverged QUIC supports: ", isempty(quic) ? "none" : join(quic, ", "), ".")
  unmatched = sort(unique("$(r["case_id"]) ($(r["calibration_status"]))" for r in rows if r["method"] == "RefitADMM" && r["calibration_matched"] != "true"))
  println(io, "- ADMM without a matched calibration: ", isempty(unmatched) ? "none" : join(unmatched, ", "), ".")
  String(take!(io))
end

length(ARGS) == 1 || error("usage: summarize.jl OUTDIR")
report = summary(read_rows(ARGS[1]))
print(report)
file = joinpath(ARGS[1], "summary.md")
open(file * ".tmp", "w") do io
  write(io, report)
end
mv(file * ".tmp", file; force=true)
