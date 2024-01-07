"""
`solve_Thomas_compact_fast!(a,d,b,y)`
in-place modifies `y` with the solution of the linear system `Ax=y` with `A` being a tridiagonal matrix with `a`,`d`,`b` being the elemets of the diagonal band.
"""
function solve_Thomas_compact_fast!(a,d,b,y) # for the inverse of the matrix
    a[1], y[1]= [a[1], y[1]]./d[1]; 
    for i in 2:length(d) a[i], y[i]= [a[i], (y[i]- b[i]*y[i-1])]./(d[i]- b[i]*a[i-1]) end;
    for i in length(d)-1:-1:1 y[i]= y[i]- a[i]*y[i+1] end
end

function get_ρ_at_z(pred, zs)
    n= 1+ length(pred) ÷2;
    h= pred[n+1:end];
    res= zeros(length(zs));
    idx= zs.<= h[1];
    res[idx].= pred[1];
    
    for (iz, z) in enumerate(zs)
        for i in 1:length(h)
            if zs[iz] > h[i] 
                res[iz]= pred[i+1];
            end
        end
    end
    idx= zs.> h[end];
    res[idx].= pred[n];
    return res;
end
"""
`fwd(ρ,h,ω)`
1D MT recursion response for a model with resistivity distribution `ρ` with layer thickness  as `h` for frequency `ω`.
"""
function fwd(ρ,h,ω)
    k= sqrt.(-im * ω* μ ./ ρ);
    R= 0 .*k;
    r= 0 .*k;
    for i in n-1:-1:1
        r[i]= (k[i]- k[i+1])/(k[i]+ k[i+1]);
        if i==n-1
            # r[i]= (k[i]- k[i+1])/(k[i]+ k[i+1]);
            
            R[i]= r[i]* exp(-im*k[i]*h[i]);
        else
            R[i]= (r[i]+ R[i+1]*exp(-im*k[i+1]*h[i+1]))/
            (1+ r[i]* R[i+1]*exp(-im*k[i+1]*h[i+1]))*
            exp(-im*k[i]*h[i]);
        end 
    end
    
    Z1= ω*μ/k[1];
    Z= -Z1*(R[1]*exp(-im*k[1]*h[1])+1)/(R[1]*exp(-im*k[1]*h[1])-1);
end

"""
`fwds(ρ,h,ω)`
1D recursion solution for a model with resistivity distribution `ρ` with layer thickness  as `h` for a vector of frequency `ωs`.
"""
fwds(ρ,h,ωs)= [fwd(ρ,h,iω) for iω in ωs];