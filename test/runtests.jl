# SPDX-License-Identifier: GPL-3.0-or-later

using DrWatson, Test
quickactivate("..", "Chimera_EGT_Kuramoto")

# Run test suite
println("Starting tests")
ti = time()

@testset "Chimera_EGT_Kuramoto tests" begin
    @testset "Utils tests" begin
	include("utils_tests.jl")
    end

    @testset "Moran tests" begin
	include("moran_tests.jl")
    end

    @testset "Postprocess tests" begin
	include("postprocess_tests.jl")
    end
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
