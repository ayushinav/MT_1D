"""
`solve_Thomas_compact_fast!(a,d,b,y)`
in-place modifies `y` with the solution of the linear system `Ax=y` with `A` being a tridiagonal matrix with `a`,`d`,`b` being the elemets of the diagonal band.
"""
function solve_Thomas_compact_fast!(a,d,b,y) # for the inverse of the matrix
    a[1], y[1]= [a[1], y[1]]./d[1]; 
    for i in 2:length(d) a[i], y[i]= [a[i], (y[i]- b[i]*y[i-1])]./(d[i]- b[i]*a[i-1]) end;
    for i in length(d)-1:-1:1 y[i]= y[i]- a[i]*y[i+1] end
end

"""
`get_ρ_at_z(model, zs)` returns the resistivities at points specified by `zs` according to the `model`. Note that `model`= [`rho`..., `'h`...]
"""
function get_ρ_at_z(model, zs)
    n= 1+ length(model) ÷2;
    h= model[n+1:end];
    res= zeros(length(zs));
    idx= zs.<= h[1];
    res[idx].= model[1];
    
    for (iz, z) in enumerate(zs)
        for i in 1:length(h)
            if zs[iz] > h[i] 
                res[iz]= model[i+1];
            end
        end
    end
    idx= zs.> h[end];
    res[idx].= model[n];
    return res;
end

"""
`plot_fig(Z; color= "blue", lable= false, line= :solid)` returns a tuple of plots of ρₐ, phase, both together
"""
function plot_fig(Z; color= "blue", lable= false, line= :solid)
    ρₐ= inv.(ω).* inv(μ).* abs.(Z).^2;
    ϕ= 180/π .*atan.(imag.(Z)./real.(Z));
    
    plt1= plot(1 ./f, ρₐ, scale= :log10,
        xlabel= "period (s)", ylabel= "ρₐ", label= lable, color= color, linestyle= line
    )

    plt2= plot(1 ./f, ϕ, xscale= :log10,
        xlabel= "period (s)", ylabel= "phase (deg)", label= lable, color= color, linestyle= line, ylim= (0,90)
    )

    plt3= plot(plt1, plt2, layout= (1,2), size= (900, 350))
    return plt1,plt2,plt3;
end

"""
`plot_fig!(plt1, plt2, plt3, Z; color= "blue", lable= false, line= :solid)` plots on a tuple of plots of ρₐ, phase, both together
"""
function plot_fig!(plt1, plt2, plt3, Z; color= "blue", lable= false, line= :solid)
    ρₐ= inv.(ω).* inv(μ).* abs.(Z).^2;
    ϕ= 180/π .*atan.(imag.(Z)./real.(Z));
    
    plot!(plt1, 1 ./f, ρₐ, scale= :log10,
        xlabel= "period (s)", ylabel= "ρₐ", label= lable, color= color, linestyle= line
    )

    plot!(plt2, 1 ./f, ϕ, xscale= :log10,
        xlabel= "period (s)", ylabel= "phase (deg)", label= lable, color= color, linestyle= line, ylim= (0,90)
    )
end