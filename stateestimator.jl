using DifferentialEquations

abstract type StateObserver end

abstract type OptimizableObserver <: StateObserver end

struct NullObserver <: StateObserver end

struct StateEstimator
    trueProblem::DEProblem
    measurement::Function
    tSpan::Tuple
    estimatorProblem::ODEProblem
    observer::StateObserver
    StateEstimator(trueProblem, measurement, tSpan) = new(trueProblem, measurement, tSpan, ODEProblem(trueProblem.f.f, zeros(size(trueProblem.u0)), tSpan), NullObserver())
    StateEstimator(trueProblem, estimatorProblem, measurement, tSpan) = new(trueProblem, measurement, tSpan, estimatorProblem, NullObserver())
    StateEstimator(trueProblem, estimatorProblem, measurement, tSpan, observer) = new(trueProblem, measurement, tSpan, estimatorProblem, observer)
end

mutable struct Luenberger <: OptimizableObserver
    L::Matrix{Float64}
    accuracy::Float64

    observe::Function
    set::Function
    bounds::Function
    format::Function

    function observe(self::Luenberger, stateEstimator::StateEstimator, ym::Vector{Float64}, y_hat::Vector{Float64}, du, u, p, t) 
        mat = self.L * (ym - y_hat)
        i = 1
        while (i <= size(du)[1])
            du[i] = du[i] + mat[i]
            i = i + 1
        end
    end

    function set(self::Luenberger, stateEstimator::StateEstimator, newValue::Matrix{Float64}, accuracy::Float64)
        self.L = newValue
        self.accuracy = accuracy
    end

    function bounds(self::Luenberger, stateEstimator::StateEstimator)
        testArr = stateEstimator.measurement(stateEstimator.trueProblem.u0)
        s = size(testArr)[1] * size(stateEstimator.trueProblem.u0)[1]
        return [fill(0,s),fill(1,s)]
    end

    function formatArr(self::Luenberger, stateEstimator::StateEstimator, x::Array)
        testArr = stateEstimator.measurement(stateEstimator.trueProblem.u0)
        formatted = reshape(x, (size(stateEstimator.trueProblem.u0)[1], size(testArr)[1]))
        return formatted
    end

    function Luenberger(L) 
        self = new(L, 0.0)
        function observe0(stateEstimator::StateEstimator, ym::Vector{Float64}, y_hat::Vector{Float64}, du, u, p, t)
            observe(self, stateEstimator, ym, y_hat, du, u, p, t)
        end
        self.observe = observe0
        function set0(stateEstimator::StateEstimator, newValue::Matrix{Float64}, accuracy::Float64)
            set(self, stateEstimator, newValue, accuracy)
        end
        self.set = set0
        function bounds0(stateEstimator::StateEstimator)
            bounds(self, stateEstimator)
        end
        self.bounds = bounds0
        function format0(stateEstimator::StateEstimator, x::Array)
            formatArr(self, stateEstimator, x)
        end
        self.format = format0
        return self
    end
    function Luenberger()
        self = new(zeros(Float64,1,1), 0.0)
        function observe0(stateEstimator::StateEstimator, ym::Vector{Float64}, y_hat::Vector{Float64}, du, u, p, t)
            observe(self, stateEstimator, ym, y_hat, du, u, p, t)
        end
        self.observe = observe0
        function set0(stateEstimator::StateEstimator, newValue::Matrix{Float64}, accuracy::Float64)
            set(self, stateEstimator, newValue, accuracy)
        end
        self.set = set0
        function bounds0(stateEstimator::StateEstimator)
            bounds(self, stateEstimator)
        end
        self.bounds = bounds0
        function format0(stateEstimator::StateEstimator, x::Array)
            formatArr(self, stateEstimator, x)
        end
        self.format = format0
        return self
    end
end

function solveEstimatedState(stateEstimator::StateEstimator) 
    trueSolution = solve(stateEstimator.trueProblem, Tsit5())

    function stateFunction(du, u, p, t)
        tv = trueSolution(t)

        du[1] = -k * u[1] * u[2]
        du[2] = -k * u[1] * u[2]
        du[3] = k * u[1] * u[2]

        ym = stateEstimator.measurement(tv)
        y_hat = stateEstimator.measurement(u)
        if (!(stateEstimator.observer isa NullObserver))
            stateEstimator.observer.observe(stateEstimator, ym, y_hat, du, u, p, t)
        end
    end
    stateProblem = ODEProblem(stateFunction, stateEstimator.estimatorProblem.u0, stateEstimator.tSpan)
    
    return solve(stateProblem, Tsit5())
end