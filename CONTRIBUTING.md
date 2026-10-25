# Development conventions

Use UpperCamelCase for types and modules and snake_case for descriptive functions,
fields and keyword arguments. Mutating functions end in `!`. Keep established
numerical abbreviations, mathematical identifiers and external interface names.
Old public names live only in explicit compatibility code.

Develop and test against local master checkouts of workspace dependencies. Run
`julia scripts/instantiate.jl` to prepare the runtime environment, and
`julia scripts/test_master.jl` for the tests.
Version numbering is deferred; retain existing package versions and UUIDs.

Run the complete tests with `julia scripts/test_master.jl`. This uses an isolated
master-based environment on Julia 1.10 and newer without changing package deps.
