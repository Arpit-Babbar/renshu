using DelimitedFiles
using Plots

files = readdir("output", join = true)

n_time_levels = length(files)
energies = zeros(n_time_levels, 2)

for (i, file) in enumerate(files)
    gamma = 1.4
    epsilon = 0.1
    open(file) do f
        t = parse(Float64,readline(f))
        data = readdlm(f)
        rho = data[:, 1]
        v = data[:, 2]
        p = rho .^ gamma
        e = p / (epsilon^2 * (gamma - 1)) .+ 0.5 * rho .* v.^2
        dx = 1.0 / length(rho)
        total_energy = sum(e * dx)
        # total_energy = sqrt(sum(e.^2 * dx))
        # @show file
        energies[i, 1] = t
        energies[i, 2] = total_energy
    end
end

energies[:, 2] ./= energies[1, 2]

scatter(energies[:,1], energies[:,2], label = "Total energy", xlabel = "Time", ylabel = "Total energy", legend = :topleft)
