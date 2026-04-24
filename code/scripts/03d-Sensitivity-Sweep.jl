# =============================================================================
# 03d-Sensitivity-Sweep.jl
#
# Sensitivity sweeps for the hybrid composer's two hyperparameters
# (peer-review P2, R2): the R²-preserve branch threshold and the
# idiosyncratic-variance floor f. Because these knobs only affect the
# hybrid composer, each inner call restricts `include_composers` to
# {"hybrid"} to keep runtime modest (≈6× faster than a full six-composer
# run per grid point).
#
# Grids:
#   R²_threshold ∈ {0.70, 0.75, 0.80, 0.85, 0.90}
#   f            ∈ {0.05, 0.10, 0.15, 0.20, 0.30}
#
# Output: data/results-thresh-{thresh}-f-{f}.jld2 per grid cell. The
# consolidated sensitivity plot is built in scripts/05-Figures.jl (new
# panel "Figure 4: threshold × floor sensitivity").
# =============================================================================

include(joinpath(@__DIR__, "..", "Include.jl"))

cfg = load_config()

R²_grid = [0.70, 0.75, 0.80, 0.85, 0.90]
f_grid  = [0.05, 0.10, 0.15, 0.20, 0.30]

@info "Sensitivity sweep" R²_grid = R²_grid f_grid = f_grid n_cells = length(R²_grid) * length(f_grid)

for thresh in R²_grid
    for floor_f in f_grid
        suffix = "-thresh-$(thresh)-f-$(floor_f)"
        @info ">>> Sensitivity cell" thresh=thresh f=floor_f
        run_composer_experiment(cfg;
            R²_threshold = thresh,
            f            = floor_f,
            include_composers = Set(["hybrid"]),
            output_suffix = suffix,
        )
    end
end

@info "Sensitivity sweep complete."
