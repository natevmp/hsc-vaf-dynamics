function completeParams!(params::Dict)
	if !haskey(params, "p")
		params["p"] = params["ϕ"] / (params["ϕ"] + params["ρ"])
	end
	if !haskey(params, "λ")
		params["λ"] = params["ϕ"] + params["ρ"]
	end
	return params
end

function freqToInd(freq, df)
    1 + Int(round(freq/df))
end


""" Finite differences """

# === first order ===
function fd1(f0, f1, dx::Float64)
    (f1 .- f0)/dx
end

function fd1(f::Array{T} where T<:Number, dx::Float64)
    fd1(f[1:end-1], f[2:end], dx)
end

# === second order forward ===
function fd2(f0, f1, f2, dx::Float64)
    ( f2 .- 2f1 .+ f1 )/dx^2
end

function fd2(f::Array{T} where T<:Number, dx::Float64)
    fd2(f[1:end-2], f[2:end-1], f[3:end], dx)
end

# === second order central ===
function cd2(f0, fu, fd, dx::Float64)
    ( fu .- 2*f0 .+ fd )/dx^2
end

function cd2(f::Array{T} where T<:Number, dx::Float64)
    cd2(f[2:end-1], f[3:end], f[1:end-2], dx)
end