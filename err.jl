using LinearAlgebra
using Plots
using ProgressMeter

include("utils.jl")
include("fd.jl")
include("fem.jl")
include("analytic.jl")

# ρ= [100, 1, 100];
# h= [1000, 4e4-1000];

ρ= [1, 10000, 10];
h= [100, 1000];

n= length(h)+1;
f= 10 .^range(-4, stop= 0, length= 41);
ω= 2π.*f;
μ= 4π*1e-7;

function get_err_fd(ρ, h, ω, dzp, mzp)
    Z_fd, nz= fd_fwd1d(ρ,h,ω, dzp, mzp);
    Z_analytical= fwd(ρ, h, ω);

    ρₐ_fd= inv.(ω).* inv(μ).* abs.(Z_fd).^2;
    ϕ_fd= 180/π .*atan.(imag.(Z_fd)./real.(Z_fd));

    ρₐ_analytical= inv.(ω).* inv(μ).* abs.(Z_analytical).^2;
    ϕ_analytical= 180/π .*atan.(imag.(Z_analytical)./real.(Z_analytical));

    err_ρₐ= 100*sum(abs.(ρₐ_analytical.- ρₐ_fd)./abs.(ρₐ_fd))/length(ω)
    err_ϕ= 100*sum(abs.(ϕ_analytical.- ϕ_fd)./abs.(ϕ_fd))/length(ω)
    
    return [err_ρₐ, err_ϕ, nz];
end

dzp_grid= reverse(10 .^collect(0:-0.1:-3));
mzp_grid= collect(2:0.5:10);

err_fd= zeros(length(dzp_grid), length(mzp_grid), 3);

@showprogress for idzp in 1:length(dzp_grid), imzp in 1:length(mzp_grid)
    err_fd[idzp, imzp,:] .= get_err_fd(ρ, h, ω, dzp_grid[idzp], mzp_grid[imzp]);
end

function get_err_fe(ρ, h, ω, mzp; n=16)
    Z_fd, nz= fem_fwd1d(ρ,h,ω, mzp,n=n);
    Z_analytical= fwd(ρ, h, ω);

    ρₐ_fd= inv.(ω).* inv(μ).* abs.(Z_fd).^2;
    ϕ_fd= 180/π .*atan.(imag.(Z_fd)./real.(Z_fd));

    ρₐ_analytical= inv.(ω).* inv(μ).* abs.(Z_analytical).^2;
    ϕ_analytical= 180/π .*atan.(imag.(Z_analytical)./real.(Z_analytical));

    err_ρₐ= 100*sum(abs.(ρₐ_analytical.- ρₐ_fd)./abs.(ρₐ_fd))/length(ω)
    err_ϕ= 100*sum(abs.(ϕ_analytical.- ϕ_fd)./abs.(ϕ_fd))/length(ω)
    
    return [err_ρₐ, err_ϕ, nz];
end

mzp_grid= collect(2:0.5:10);
ngrid= collect(5:2:500);

err_fe= zeros(length(mzp_grid), length(ngrid), 3);

@showprogress for imzp in 1:length(mzp_grid), inn in 1:length(ngrid)
    err_fe[imzp, inn,:] .= get_err_fe(ρ, h, ω, mzp_grid[imzp], n= ngrid[inn]);
end

p1= plot(err_fe[1,:,3], err_fe[1,:,1], scale=:log10, label= "fem") #
plot!(p1, err_fd[:,1,3], err_fd[:,1,1], scale=:log10, label= "fd") #
plot!(p1, title= "2 x δₘₐₓ", ylims= (0.01, 1e4), yticks= [0.01, 0.1, 1, 10, 100, 1000, 1e4], xticks= [1, 10, 100, 1000, 1e4, 1e5, 1e6])
plot!(p1, xlabel= " number of nodes", ylabel= "% error in ρₐ")

p2= plot(err_fe[end,:,3], err_fe[end,:,1], scale=:log10, label= "fem") #
plot!(p2, err_fd[:,end,3], err_fd[:,end,1], scale=:log10, label= "fd") #
plot!(p2, title= "10 x δₘₐₓ", ylims= (0.01, 1e4), yticks= [0.01, 0.1, 1, 10, 100, 1000, 1e4], xticks= [1, 10, 100, 1000, 1e4, 1e5, 1e6])
plot!(p2, xlabel= " number of nodes", ylabel= "% error in ρₐ")

p3= plot(err_fe[1,:,3], err_fe[1,:,1], scale=:log10, label= "2 x δₘₐₓ") #
plot!(p3, err_fe[end,:,3], err_fe[end,:,1], scale=:log10, label= "10 x δₘₐₓ") #
plot!(p3, title="fe:δₘₐₓ comparison", ylims= (0.01, 1e4), yticks= [0.01, 0.1, 1, 10, 100, 1000, 1e4], xticks= [1, 10, 100, 1000, 1e4, 1e5, 1e6])
plot!(p3, xlabel= " number of nodes", ylabel= "% error in ρₐ")

plot(p1,p2,p3, layout= (1,3), size= (1000, 400), margin=5Plots.mm)
savefig("appres_err.png")


p4= plot(err_fe[1,:,3], err_fe[1,:,2], scale=:log10, label= "fem") #
plot!(p4, err_fd[:,1,3], err_fd[:,1,2], scale=:log10, label= "fd") #
plot!(p4, title= "2 x δₘₐₓ", ylims= (0.01, 1e4), yticks= [0.01, 0.1, 1, 10, 100, 1000, 1e4], xticks= [1, 10, 100, 1000, 1e4, 1e5, 1e6])
plot!(p4, xlabel= " number of nodes", ylabel= "% error in ϕ")

p5= plot(err_fe[end,:,3], err_fe[end,:,2], scale=:log10, label= "fem") #
plot!(p5, err_fd[:,end,3], err_fd[:,end,2], scale=:log10, label= "fd") #
plot!(p5, title= "10 x δₘₐₓ", ylims= (0.01, 1e4), yticks= [0.01, 0.1, 1, 10, 100, 1000, 1e4], xticks= [1, 10, 100, 1000, 1e4, 1e5, 1e6])
plot!(p5, xlabel= " number of nodes", ylabel= "% error in ϕ")

p6= plot(err_fe[1,:,3], err_fe[1,:,2], scale=:log10, label= "2 x δₘₐₓ") #
plot!(p6, err_fe[end,:,3], err_fe[end,:,2], scale=:log10, label= "10 x δₘₐₓ") #
plot!(p6, title="fe:δₘₐₓ comparison", ylims= (0.01, 1e4), yticks= [0.01, 0.1, 1, 10, 100, 1000, 1e4], xticks= [1, 10, 100, 1000])
plot!(p6, xlabel= " number of nodes", ylabel= "% error in ϕ")

plot(p4,p5,p6, layout= (1,3), size= (1000, 400), margin=5Plots.mm)
savefig("phase_err.png")
