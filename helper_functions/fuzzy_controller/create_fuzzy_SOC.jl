using FuzzyLogic

fis = @mamfis function fuzzy_control(p_now, p_trend, p_volatility, SOC_now)::SOC_target
    p_now := begin
        domain = -50:950 # TODO change to min max 
        cheap = TrapezoidalMF(-100, -50, 0, 450)
        average = TrapezoidalMF(-50, 300, 450, 950)
        expensive = TriangularMF(450, 950, 1500)
        trapezoidalMF
    end

    p_trend := begin
        domain = -50:950 # TODO change to min max 
        cheap = TriangularMF(-100, -50, 450)
        average = TriangularMF(-50, 450, 950)
        expensive = TriangularMF(450, 950, 1500)
    end

    p_volatility := begin
        domain = 0:10 # TODO change to min max 
        low = TriangularMF(-5, 0, 5)
        average = TriangularMF(0, 5, 10)
        high = TriangularMF(5, 10, 15)
    end

    SOC_now := begin
        domain = 0:1
        low = TriangularMF(-1.0, 0.0, 0.5)
        middle = TriangularMF(0.0, 0.5, 1.0)
        high = TriangularMF(0.5, 1.0, 2.0)
    end

    SOC_target := begin
        domain = 0:1
        low = TriangularMF(-1.0, 0.0, 0.5)
        middle = TriangularMF(0.0, 0.5, 1.0)
        high = TriangularMF(0.5, 1.0, 2.0)
    end

    p_now == cheap || p_trend == expensive --> SOC_target == low
    p_now == expensive || p_volatility == high --> SOC_target == high
    p_trend == cheap || p_volatility == low --> SOC_target == middle
end

# output of compilefis should be used in the controller because it's a lot faster and 
# doesn't need FuzzyLogic as a dependecy
compilefis(fis)