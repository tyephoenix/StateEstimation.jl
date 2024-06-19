

abstract type Control end

struct Valve <: Control
    α::Float64
    
    function control(self::Valve, setPoints::Array{Float64}, du, u, p, t)


        return self.α * u
    end
    Valve(α) => new(α)
end

struct Pump <: Control
    
end 