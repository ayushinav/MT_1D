"""
`fd_fwd1d(ρ,h,ω, dz_p, max_z_p)`
finite difference for 1D MT for resistivity distribution `ρ` with thickness `h` at a frequency `ω`, 
`dz_p` is the percentage of smallest skin depth to be used for discretization of model,
`max_z_p` is the multiplication factor for the maximum skin depth for the deepest extent of the model space.
returns a tuple of impedance at the surface `Z` and the number of node points used `nz`.
"""
function fd_fwd1d(ρ,h,ω::Real, dz_p, max_z_p)
    dz= 500*sqrt(minimum(ρ)*2π/ω)*dz_p; # automate this, and then use variable grid
    zgrid= 1:dz:500*sqrt(maximum(ρ)*2π/ω)*max_z_p; # 10 times the max skin depth for that frq
    nz= length(zgrid);
    ρgrid= get_ρ_at_z([ρ..., h...], zgrid);
    
    E₀= 1;
    E= 0im .+ zeros(nz);
    E[1]= -E₀; 
    a= 0im .+ ones(nz);
    a[end]= 0;
    b= 0im .+ ones(nz);
    b[1]= 0;
    d= -2 .- im*μ.*ω.*inv.(ρgrid).*(dz)^2;

    solve_Thomas_compact_fast!(a,d,b,E);
    H= -1/(im*μ*ω) .*(E[2]-E₀)/(2dz);
    Z= E[1] *inv.(H);

    return Z,nz;
end

"""
`fd_fwd1d(ρ,h,ω, dz_p, max_z_p)`
finite difference response for 1D MT for resistivity distribution `ρ` with thieckness `h` at frequencies given by vector `ω`, 
`dz_p` is the percentage of smallest skin depth to be used for discretization of model,
`max_z_p` is the multiplication factor for the maximum skin depth for the deepest extent of the model space.
returns a tuple of impedance at the surface impdenances `Z` and the number of node points used `nz` for all the frequencies.
"""
function fd_fwd1d(ρ,h,ω::AbstractVector,dzp, mzp)
    arr= [fd_fwd1d(ρ,h,iω, dzp, mzp) for iω in ω];
    a1= [arr[i][1] for i in 1:length(arr)];
    nz= arr[1][2];
    return a1,nz;
end