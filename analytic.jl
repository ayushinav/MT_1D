"""
`fwd(ρ,h,ω)`
1D MT recursion response for a model with resistivity distribution `ρ` with layer thickness  as `h` for frequency `ω`.
"""
function fwd(ρ,h,ω::Real)
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
fwd(ρ,h,ω::AbstractVector)= [fwd(ρ,h,iω) for iω in ω];