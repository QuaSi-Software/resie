using FuzzyLogic

fis = @mamfis function fuzzy_control(p_now, p_future, p_stability, SOC_now)::SOC_target
    p_now := begin
        domain = -50:950 # TODO change to min max 
        cheap = TriangularMF(-100, -50, 450)
        average = TriangularMF(-50, 450, 950)
        expensive = TriangularMF(450, 950, 1500)
        trapezoidalMF
    end

    p_future := begin
        domain = -50:950 # TODO change to min max 
        cheap = TriangularMF(-100, -50, 450)
        average = TriangularMF(-50, 450, 950)
        expensive = TriangularMF(450, 950, 1500)
    end

    p_stability := begin
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

    p_now == cheap || p_future == high --> SOC_target == low
    p_now == expensive || p_stability == high --> SOC_target == high
    p_future == low || p_stability == low --> SOC_target == middle
end

# output of compilefis should be used in the controller because it's a lot faster and 
# doesn't need FuzzyLogic as a dependecy
compilefis(fis)