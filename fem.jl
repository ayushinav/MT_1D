"""
xrange= vector to be updated, of the same length as Δx
Δx= interval to the left of the grid point. The first is 0, making it the same length as the number of elements.
"""
function build_grid!(xrange, Δx)
    for i in 2:length(Δx)
        xrange[i]= xrange[i-1]+ Δx[i];
    end
end

"""
`K1(dx)` returns the tuple of vector of the tridiagonal of the matrix `A` for discretization `dx` as `a`,`d`,`b`
where `d` is the vector of diagonal elements, `a` and `b` are the elements above and below the diagonal respectively
"""
function K1(dx)
    N= length(dx)

    a= zeros(N) .+0im; 
    d= zeros(N) .+0im;
    b= zeros(N) .+0im;

    for i in 2:N
        i<N && (d[i]= 1/dx[i]+ 1/dx[i+1]);
        a[i-1]= -1/dx[i]; # implicitly makes a[end]= 0
        b[i]= -1/dx[i]; # implicitly makes b[1]= 0
    end
    d[1]= 1/dx[2];
    d[end]= 1/dx[end];

    return a,d,b;
end

# k2 has zero at the first index, like dx. Even if it does not, it won't matter because we never access k2[1] unless calculating d[1] which also gets eliminated because of Dirichlet BC.
"""
`K2(dx)` returns the tuple of vector of the tridiagonal of the matrix `B` for discretization `dx` as `a`,`d`,`b`
where `d` is the vector of diagonal elements, `a` and `b` are the elements above and below the diagonal respectively
"""
function K2(dx, k2)
    N= length(dx);

    a= zeros(N) .+0im;
    d= zeros(N) .+0im;
    b= zeros(N) .+0im;

    for i in 2:N
        i<N && (d[i]= k2[i]*dx[i]/3+ k2[i+1]*dx[i+1]/3);
        a[i-1]= k2[i]*dx[i]/6;
        b[i]= k2[i]*dx[i]/6;
    end
    d[1]= k2[2]*dx[2]/3;
    d[end]= k2[end]*dx[end]/3;
    
    return a,d,b;
end

"""
`get_grid_for_layer(L; n= 16)` returns a set of grid points in a layer such that the separation between the points increases arithmatically towards the middle and then decreases after that with the same slope where `L` is the thickness of the layer constructing 2`n`-1 grid points.
"""
function get_grid_for_layer(L; n= 16)
    d= L/(n-1)/n;
    arr= zeros(n);
    an= 0;
    for i in 1:n
        arr[i]= an;
        an+=d;
    end
    b= arr[2:end]
    return [b..., reverse(b)...];
end

"""
`make_dz(ρ, h, ω, mzp; n= 16)` discretizes the resistivity distribution `ρ` with thickness `h` where `mzp` is used to get the maximum extent of the grid space for the frequency `ω` constructing 2`n`-1 grid points.
"""
function make_dz(ρ, h, ω, mzp; n= 16)
    n2= 2n-2;
    max_depth= 500*sqrt(maximum(ρ)*2π/ω)*mzp;
    ngrid= (n2)*(length(h)+1);
    arr= zeros(n2*length(h)+n2);
    for i in 1:length(h)
        arr[n2*(i-1)+1:n2*(i)].= get_grid_for_layer(h[i], n=n);
    end
    arr[end- n2+1: end].= get_grid_for_layer(max_depth- sum(h), n=n);

    return arr;
end

"""
`fem_fwd1d(ρ, h, ω::Real, mzp; n= 16)`
finite element response for 1D MT for resistivity distribution `ρ` with thickness `h` at a frequency `ω`, 
`mz_p` is the multiplication factor for the maximum skin depth for the deepest extent of the model space,
2`n`-1 grid points per layer will be constructed.
returns a tuple of impedance at the surface `Z` and the number of node points used `nz`.
"""
function fem_fwd1d(ρ, h, ω::Real, mzp; n= 16)
    max_depth= 500*sqrt(maximum(ρ)*2π/ω)*mzp;

    dz= make_dz(ρ, h, ω, mzp, n= n)
    zgrid= zero(dz);
    build_grid!(zgrid, dz);
    nz= length(zgrid);
    N= length(zgrid);
    
    ρgrid= get_ρ_at_z([ρ..., h...], zgrid);
    k2= im* ω* μ* inv.(ρgrid);
    
    a1,d1,b1= K1(dz)
    k22= [0. + 0im, k2...]; # because the element before node 0 will not have any value
    a2,d2,b2= K2(dz, k22);
    
    a= a1.+a2;
    d= d1.+d2;
    b= b1.+b2;
    
    E₀= 1e4 + 0im;
    E= zeros(N-2) .+0im;
    E[1]= -b[2]*E₀;
    E[end]= -a[end-1]*0;
    
    # Dirichlet BC
    
    deleteat!(a, [1, N-1])
    deleteat!(d, [1, N])
    deleteat!(b, [2, N]);
    
    solve_Thomas_compact_fast!(a,d,b,E);

    H= -1/(im*μ*ω) *(E[1]-E₀)/(dz[2]);
    Z= E₀ *inv.(H);
    return Z, nz;
end

"""
`fem_fwd1d(ρ, h, ω::Real, mzp; n= 16)`
finite element response for 1D MT for resistivity distribution `ρ` with thickness `h` at frequencies given by `ω`, 
`mz_p` is the multiplication factor for the maximum skin depth for the deepest extent of the model space,
2`n`-1 grid points per layer will be constructed.
returns a tuple of impedances at the surface `Z` and the number of node points used `nz` for all the frequencies.
"""
function fem_fwd1d(ρ,h,ω::AbstractVector, mzp; n= 16)
    arr= [fem_fwd1d(ρ,h,iω, mzp, n=n) for iω in ω];
    a1= [arr[i][1] for i in 1:length(arr)];
    nz= arr[1][2];
    return a1,nz;
end