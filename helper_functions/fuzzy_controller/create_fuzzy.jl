using FuzzyLogic

fis = @mamfis function fuzzy_control(p_stock, p_reserve_pos, p_reserve_neg)::power_reserve
    p_stock := begin
        domain = -50:950 # TODO change to min max 
        cheap = TriangularMF(-100, -50, 450)
        average = TriangularMF(-50, 450, 950)
        expensive = TriangularMF(450, 950, 1500)
    end

    p_reserve_pos := begin
        domain = -50:950 # TODO change to min max 
        low = TriangularMF(-100, -50, 450)
        average = TriangularMF(-50, 450, 950)
        high = TriangularMF(450, 950, 1500)
    end

    p_reserve_neg := begin
        domain = -50:950 # TODO change to min max 
        low = TriangularMF(-100, -50, 450)
        average = TriangularMF(-50, 450, 950)
        high = TriangularMF(450, 950, 1500)
    end

    power_reserve := begin
        domain = -1:1
        strong_neg = TriangularMF(-2.0, -1.0, -0.5)
        weak_neg = TriangularMF(-1.0, -0.5, 0.0)
        zero = TriangularMF(-0.5, 0.0, 0.5)
        weak_pos = TriangularMF(0.0, 0.5, 1.0)
        strong_pos = TriangularMF(0.5, 1.0, 2.0)
    end

    p_stock == cheap || p_reserve_neg == high --> power_reserve == strong_neg
    p_stock == expensive || p_reserve_pos == high --> power_reserve == strong_pos
    p_reserve_neg == low || p_reserve_pos == low --> power_reserve == zero
end

# output of compilefis should be used in the controller because it's a lot faster and 
# doesn't need FuzzyLogic as a dependecy
compilefis(fis)