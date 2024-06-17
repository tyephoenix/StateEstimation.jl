using DifferentialEquations
using Optimization
using OptimizationNLopt
using PolyChaos

include("stateestimator.jl")

abstract type SE_Optimization end

mutable struct SE_NelderMead <: SE_Optimization

    interval::Int
    maxTime::Float64

    optimize::Function

    function optimize(self::SE_NelderMead, stateEstimator::StateEstimator)
        if (stateEstimator.observer isa OptimizableObserver)
            pF = [[],-1]
            function optFunction(du, u, p, t)
                stateEstimator.estimatorProblem.f.f(du,u,p,t)

                y_hat = stateEstimator.measurement(u)
                stateEstimator.observer.set(stateEstimator, p[1], 0.0)
                stateEstimator.observer.observe(stateEstimator, p[2], y_hat, du, u, p, t)
            end
            optProblem = ODEProblem(optFunction, stateEstimator.estimatorProblem.u0, stateEstimator.tSpan)

            cbs = []
            tStops = []
            i = 1
            while i*self.interval <= stateEstimator.tSpan[2]
                timeIndex = i*self.interval
                preTimeIndex = (i-1)*self.interval
                function c(u,t,intregrator) 
                    if (t == timeIndex)
                        println("Optimizing at $(timeIndex)...")
                    end
                    t == timeIndex
                end
                function a!(intregrator)
                    ti = time()
                    map = Dict()
                    function opt(x, p) 
                        ym = stateEstimator.measurement(intregrator.u)
                        tF = stateEstimator.observer.format(stateEstimator, x)
                        tv  = intregrator.u
                        eSol = solve(optProblem, Tsit5(), p=[tF, ym], verbose=false)
                        ev = eSol(timeIndex)
                        err = tv-ev
                        sq = sum(e -> e^2, err)
                        if (pF[2] >= 0)
                            peSol = solve(optProblem, Tsit5(), p=[pF[1], ym], verbose=false)
                            pev = peSol(timeIndex)
                            peErr = tv-pev
                            peSq =  sum(e -> e^2, peErr)
                            map[sq] = [sq, peSq]
                        else 
                            map[sq] = [sq, 0]
                        end
                        return sq
                    end
                    bounds = stateEstimator.observer.bounds(stateEstimator)
                    x0 = bounds[1]
                    p = [1.0]
                    funca = OptimizationFunction(opt)
                    probl = Optimization.OptimizationProblem(funca, x0, p, lb = bounds[1], ub = bounds[2])
                    sol = solve(probl, NLopt.LN_NELDERMEAD(), maxtime=self.maxTime)
                    nF = stateEstimator.observer.format(stateEstimator, sol.u)
                    s = map[sol.objective]
                    if (pF[2] < 0) 
                        pF[1] = stateEstimator.observer.format(stateEstimator, sol.u)
                        pF[2] = s[1]
                    else 
                        tw = s[1] + s[2]
                        itw = (tw/s[2]) + (tw/(s[1]))
                        v = ((tw/s[2])/itw)*pF[1] + ((tw/(s[1]))/itw)*(nF)
                        pF[1] = v
                        acc = ((tw/s[2])/itw)*pF[2] + ((tw/(s[1]))/itw)*(s[1])
                        pF[2] = acc
                    end
                    println("SEOptimization: Discrete Event at $(timeIndex). Elapsed time = $(time() - ti).")
                end
                cb = DiscreteCallback(c,a!) 
                push!(cbs,cb)
                push!(tStops, timeIndex)
                i = i + 1
            end

            trueSolution = solve(stateEstimator.trueProblem, Tsit5(), callback=CallbackSet(cbs...), tstops=tStops)

            println("SEOptimization: Finished with $(pF[1]) with accuracy $(pF[2]).")
            stateEstimator.observer.set(stateEstimator, pF[1], pF[2])
        end
    end

    function SE_NelderMead(interval, maxTime)
        self = new(interval, maxTime)
        function optimize0(stateEstimator::StateEstimator)
            optimize(self, stateEstimator)
        end
        self.optimize = optimize0
        return self
    end
end

mutable struct SE_PolyChaos <: SE_Optimization

    interval::Int
    n::Int
    d::Int

    optimize::Function

    function optimize(self::SE_PolyChaos, stateEstimator::StateEstimator)
        if (stateEstimator.observer isa OptimizableObserver)
            pF = [[],-1]
            function optFunction(du, u, p, t)
                stateEstimator.estimatorProblem.f.f(du,u,p,t)

                y_hat = stateEstimator.measurement(u)
                stateEstimator.observer.set(stateEstimator, p[1], 0.0)
                stateEstimator.observer.observe(stateEstimator, p[2], y_hat, du, u, p, t)
            end
            optProblem = ODEProblem(optFunction, stateEstimator.estimatorProblem.u0, stateEstimator.tSpan)

            cbs = []
            tStops = []
            i = 1
            while i*self.interval <= stateEstimator.tSpan[2]
                timeIndex = i*self.interval
                preTimeIndex = (i-1)*self.interval
                function c(u,t,intregrator) 
                    if (t == timeIndex)
                        println("Optimizing at $(timeIndex)...")
                    end
                    t == timeIndex
                end
                function a!(intregrator)
                    ti = time()
                    map = Dict()
                    function opt(x) 
                        ym = stateEstimator.measurement(intregrator.u)
                        tF = stateEstimator.observer.format(stateEstimator, [x...])
                        tv  = intregrator.u
                        eSol = solve(optProblem, Tsit5(), p=[tF, ym], verbose=false)
                        ev = eSol(timeIndex)
                        err = tv-ev
                        sq = sum(e -> e^2, err)
                        if (pF[2] >= 0)
                            peSol = solve(optProblem, Tsit5(), p=[pF[1], ym], verbose=false)
                            pev = peSol(timeIndex)
                            peErr = tv-pev
                            peSq =  sum(e -> e^2, peErr)
                            map[sq] = [sq, peSq]
                        else 
                            map[sq] = [sq, 0]
                        end
                        return sq
                    end
                    bounds = stateEstimator.observer.bounds(stateEstimator)
                    uniform01 = Uniform01OrthoPoly(self.d)
                    samples = []
                    j = 1
                    while j <= size(bounds[1])[1]
                        push!(samples, evaluatePCE(convert2affinePCE(bounds[1][j],bounds[2][j],uniform01), sampleMeasure(self.n,uniform01),uniform01))
                        j = j + 1
                    end
                    lErr = Inf
                    lEle = []
                    z = 1
                    for elements in Iterators.product(samples...)
                        sq = opt(elements)
                        if (sq < lErr) 
                            lErr = sq
                            lEle = elements
                        end
                        z = z + 1
                        println(z)
                    end
                    nF = stateEstimator.observer.format(stateEstimator, [lEle...])
                    s = map[lErr]
                    if (pF[2] < 0) 
                        pF[1] = stateEstimator.observer.format(stateEstimator, sol.u)
                        pF[2] = s[1]
                    else 
                        tw = s[1] + s[2]
                        itw = (tw/s[2]) + (tw/(s[1]))
                        v = ((tw/s[2])/itw)*pF[1] + ((tw/(s[1]))/itw)*(nF)
                        pF[1] = v
                        acc = ((tw/s[2])/itw)*pF[2] + ((tw/(s[1]))/itw)*(s[1])
                        pF[2] = acc
                    end
                    println("SEOptimization: Discrete Event at $(timeIndex). Elapsed time = $(time() - ti).")
                end
                cb = DiscreteCallback(c,a!) 
                push!(cbs,cb)
                push!(tStops, timeIndex)
                i = i + 1
            end

            trueSolution = solve(stateEstimator.trueProblem, Tsit5(), callback=CallbackSet(cbs...), tstops=tStops)

            println("SEOptimization: Finished with $(pF[1]) with accuracy $(pF[2]).")
            stateEstimator.observer.set(stateEstimator, pF[1], pF[2])
        end
    end

    function SE_PolyChaos(interval, n, d)
        self = new(interval, n, d)
        function optimize0(stateEstimator::StateEstimator)
            optimize(self, stateEstimator)
        end
        self.optimize = optimize0
        return self
    end
end