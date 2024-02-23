# Functions for plotting the SRG of an LTI transfer function and a matrix, in the complex plane and on the Poincare disc.
# Tom Chaffey, 23/02/2024, adapted from earlier code.
using Plots
using LinearAlgebra
using Printf
using Interpolations
using ControlSystems
using NumericalRange

include("./Graham_scan.jl")

######## SRG plotting functions ################

# Plots the SRG of an LTI transfer function by computing the convex hull in hyperbolic space, then transforming back.  Uses "Graham_scan.jl".
function SRG_LTI(G)
        r, i, w = nyquist(G)
        nyq = vec(r) .+ im.*abs.(vec(i))
        nyq_p = BeltramiKlein.(nyq)
        srg_p = Graham_scan(nyq_p)
        push!(srg_p, last(srg_p))
        srg =  BeltramiKleinInverse.(srg_p)
        p = plot(real(nyq), imag(nyq), aspect_ratio = 1, framestyle = :origin, color = :cornflowerblue)
        for i = 1:length(srg) - 1
                am = arc_min!(srg[i], srg[i+1], p)
        end
        return p

end

# Plots the SRG of an LTI transfer function in the Poincare plane.  Uses "Graham_scan.jl".
function SRG_LTI_Poincare(G)
        r, i, w = nyquist(G)
        nyq = vec(r) .+ im.*abs.(vec(i))
        nyq_p = BeltramiKlein.(nyq)
        p = scatter(real(nyq_p), imag(nyq_p), aspect_ratio = 1, framestyle = :origin, color = :red, legend = false)
        srg_p = Graham_scan(nyq_p)
        push!(srg_p, last(srg_p))
        plot!(p, real.(srg_p), imag.(srg_p), aspect_ratio = 1, framestyle = :origin, color = :blue, legend = false)
        plot_circle!(p)
        return p
end

# Plot the SRG of a matrix, by taking the convex hull of the numerical range under the BK mapping.  See Pates, "The SRG of a linear operator".
function SRG_matrix(A)
        (r, e) = NumericalRange.nrange(BeltramiKlein(A), noplot=true, thmax = 128)
        plot(real(BeltramiKleinInverse.(r)), imag(BeltramiKleinInverse.(r)), framestyle=:origin, aspect_ratio=1, legend=false, color= :cornflowerblue)
end
############ Utility functions ##################

# plot unit circle
function plot_circle!(p)
    x = LinRange(-1.0, 1.0, 50)
    y = (1.0 .- x.^2).^(0.5)
    plot!(p, x, y, legend=false, color = :cornflowerblue)
    plot!(p, x, -y, legend=false, color = :cornflowerblue)
end

# plot arc_min between two complex numbers
function arc_min!(z1, z2, p)
        arcmin = arc_min(z1, z2)
        plot!(p, real(arcmin), imag(arcmin), legend=false, color = :cornflowerblue)
        display(p)
end

function arc_min(z1, z2)
        x1 = real(z1)
        x2 = real(z2)
        y1 = imag(z1)
        y2 = imag(z2)

        # edge and corner case
        if z1 == z2
                return 0 # arcmin is single point
        elseif x1 == x2 # arcmin is line
                return arcmin = x1 .+ im .* range(y1, length=5, stop=y2)
                #return plot!(p, real(arcmin), imag(arcmin))
        end
        
        # centre of circle on real axis 
        xc = (y1 - y2)*(y1 + y2)/2/(x1 - x2) + (x1 + x2)/2

        # radius
        r = sqrt((x1 - xc)^2 + y1^2)

        # start at the point with lowest argument
        θ1 = 0
        θ2 = 0
        θ1 = angle(z1 - xc + im*0)
        θ2 = angle(z2 - xc + im*0)

        arcmin = []
        # generate arcmin
        for ϕ = range(θ1, length=50, stop=θ2)#θ1:0.1:θ2
                append!(arcmin, xc + r*exp(im*ϕ))
        end

        return arcmin
end

function BeltramiKlein(A)
        return (I + A'*A)^(-0.5)*(A' - im*I)*(A - im*I)*(I + A'*A)^(-0.5)
end

function BeltramiKleinInverse(z)
        return (imag(z) - im*sqrt(1 - z*conj(z)))/(real(z) - 1)

end

