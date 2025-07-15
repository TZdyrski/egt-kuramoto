# 2024_EGT_Kuramoto
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/TZdyrski/egt-kuramoto/HEAD?urlpath=pluto%2Fopen%3Furl%3Dhttps%253A%252F%252Fraw.githubusercontent.com%252FTZdyrski%252Fegt-kuramoto%252Frefs%252Fheads%252Fmain%252Fnotebooks%252FEKT_Plots.jl)

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> 2024_EGT_Kuramoto

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "2024_EGT_Kuramoto"
```
which auto-activate the project and enable local path handling from DrWatson.
