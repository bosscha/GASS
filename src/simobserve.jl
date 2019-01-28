### observation simulation and beam fitting

## baseline coordinates
## arr: DataFrames with at least :X, :Y, :Z
## check at https://web.njit.edu/~gary/728/Lecture6.html

function calc_baselines(arr::AbstractDataFrame , lat=-23.0262015)
    
    coordinates= vcat(hcat(convert(Vector{Float64},arr[:X]))' ,
        hcat(convert(Vector{Float64},arr[:Y]))' ,
        hcat(convert(Vector{Float64},arr[:Z]))')
    rot= [0 -sind(lat) cosd(lat) ; 1 0 0 ; 0  cosd(lat)  sind(lat)]
    xyz= *(rot,coordinates)
    
    Nc= size(xyz)
    Nbl= convert(Int64,(Nc[2]*(Nc[2]-1))/2)   
    bl= zeros(Float64,3,2Nbl)   ## symmetry of the bl
    
    ibl= 1
    for i in 1:Nc[2]
        for j in i+1:Nc[2]
            bl[1,ibl]= xyz[1,i]-xyz[1,j] ; bl[1,ibl+Nbl]= xyz[1,j]-xyz[1,i]
            bl[2,ibl]= xyz[2,i]-xyz[2,j] ; bl[2,ibl+Nbl]= xyz[2,j]-xyz[2,i]
            bl[3,ibl]= xyz[3,i]-xyz[3,j] ; bl[3,ibl+Nbl]= xyz[3,j]-xyz[3,i]
            ibl += 1
        end
    end
    return(bl)
end

####
### bl: baselines computed in calc_baselines
### h0: hour angle in degrees
### d0: source declination in degrees
### only snapshot still...
function calc_uv(bl, h0=0, d0=-50)
    
    h0 *= 15   #hours to degree
    proj= [sind(h0) cosd(h0) 0 ; -sind(d0)cosd(h0) sind(d0)sind(h0) cosd(d0) ; cosd(d0)cosd(h0) -cosd(d0)sind(h0) sind(d0)]
    uv= *(proj, bl)
    return(uv)
end

function calc_uv_coverage
    ## will compute the uv_coverage with an observing and sampling time ..
    ## for now only a snapshot.
    ## ...
end

################################
## setup the npix,fov, freq, etc for the beam extraction from the uv
## coverage. It is based on the uniform weighting 2/FOV for the uv-cell
##
function calc_dirtybeam(uv, npix=511, sizefit=255 ; robust=2)
    c= 299792458.0  ; ν= 100*1e9 ; λ= c / ν
    D= 12 ; PB= 1.13 * λ / D
    duv= 2/(PB*npix) ; ds= duv/λ
    dr= 180*3600/(π*ds*npix)
    lmax= maximum(abs.(uv))

    if (duv*npix/2) < lmax
        println("##UV-gridding error, npix too small...")
    end
    
    binuv= range(-duv*npix/2, stop=duv*npix/2, length=npix)
    Hh= fit(Histogram, (uv[1,:] , uv[2,:]) , (binuv , binuv))
    H= Hh.weights

    ## Robust weighting R (Briggs)
    ## R=-2 : uniform
    ## R=2  : natural
    ## R=0.5: "briggs"
    Hsum= sum(H)
    H2sum= sum(1 ./ H[H .> 0].^2)/Hsum
    f2= (5 * 10.0^(-robust))^2 * H2sum 
    H2 = H ./ (1 .+ f2 .* H)
 
    h= fft(H2)
    hshift= fftshift(abs.(real.(h)))

    center= npix/2
    coord_a= convert(Int64,center-sizefit/2)
    coord_b= convert(Int64,center+sizefit/2)
    
    hshiftcut= hshift'[coord_a:coord_b, coord_a:coord_b]
    return(hshiftcut , dr)
end

## fitting of the Gaussian 2D for the beam
## assuming offset=0 for the beam

function gaussian2D(x,y, amplitude, xo, yo, σx, σy, θ)
    
    a = (cosd(θ)^2)/(2σx^2) + (sind(θ)^2)/(2*σy^2)
    b = -(sind(2θ))/(2σx^2) - (sind(2θ))/(2σy^2)
    c = (sind(θ)^2)/(2σx^2) + (cosd(θ)^2)/(2σy^2)
    gauss2d= amplitude*exp( -( a*((x-xo)^2) + 2b*(x-xo)*(y-yo) + c*((y-yo)^2)) )
    return(gauss2d)
end

#################################
function fit_beam(beam , dr)
    bs= size(beam)
    h= reshape(beam,(bs[1]*bs[2],))
    
    xy= zeros(Float64,bs[1]*bs[2], 2)
    
    for i in 1:bs[1]
        for j in 1:bs[2]
            xy[(i-1)*bs[1]+j,1]= i
        end
    end
    
    for i in 1:bs[1]
        for j in 1:bs[2]
            xy[(i-1)*bs[1]+j,2]= bs[2]+1-j
        end
    end        
    
    @. multimodel(x, p) = gaussian2D(x[:,1], x[:,2], p[1], p[2] , p[3] , p[4] , p[5], p[6])
    p0 = [maximum(beam) , bs[1]/2 , bs[2]/2 , 1, 1, 0]
    f= curve_fit(multimodel, xy, h, p0)
    
    gaussian_sigma_to_fwhm=2.3548200450309493
    
    bx= f.param[4]*dr*gaussian_sigma_to_fwhm
    by= f.param[5]*dr*gaussian_sigma_to_fwhm  
    e= max(bx/by, by/bx)
    ar= sqrt(bx*by)
    sidelobe= maximum(f.resid * 100 / f.param[1])
    
    b = synthbeam(bx, by,ar,e, sidelobe)
    
    return(b)
end

#################
### MRS in arcsec at 100 GHz
function calc_mrs(uv)
    minbl = minimum(uv[:,1].^2 .+ uv[:,2].^2)
    lmin= sqrt(minbl)
    ν= 100
    
    mrs= 37100 / (lmin*ν)
    return(mrs)
end
