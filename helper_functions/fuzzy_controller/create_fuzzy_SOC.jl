using FuzzyLogic

fis = @mamfis function fuzzy_control(p_now, p_trend, p_volatility, SOC_now)::{P2H, SOC_target}
    p_now := begin
        domain = -136:936 # TODO change to min max 
        low = TrapezoidalMF(-137, -136, 55, 80)
        med = TriangularMF(55, 80, 100)
        high = TrapezoidalMF(80, 100, 1000, 1001)
        
    end

    p_trend := begin
        domain = -257:246 # TODO change to min max 
        falling = TrapezoidalMF(-258, -257, -7, 6)
        rising = TrapezoidalMF(-7, 6, 246, 247)
    end

    p_volatility := begin
        domain = 0:265 # TODO change to min max 
        low = TrapezoidalMF(-1, 0, 4, 15)
        high = TrapezoidalMF(4, 15, 265, 266)
    end

    SOC_now := begin
        domain = 0:1
        low = TriangularMF(-0.5, 0.0, 0.5)
        high = TriangularMF(0.5, 1.0, 1.5)
    end

    P2H := begin
        domain = -1:1
        lower = TriangularMF(-1.5, -1.0, 0.0)
        keep = TriangularMF(-1.0, 0.0, 1.0)
        rise = TriangularMF(0.0, 1.0, 1.5)
    end

    SOC_target := begin
        domain = 0:1
        low = TriangularMF(-0.5, 0.0, 0.5)
        med = TriangularMF(0.0, 0.5, 1.0)
        high = TriangularMF(0.5, 1.0, 1.5)
    end

    p_now == low && p_trend == falling && p_volatility == low  && SOC == low  --> (P2H == keep , SOC_target == med)
    p_now == low && p_trend == falling && p_volatility == low  && SOC == high --> (P2H == keep , SOC_target == med)
    p_now == low && p_trend == falling && p_volatility == high && SOC == low  --> (P2H == rise , SOC_target == high)
    p_now == low && p_trend == falling && p_volatility == high && SOC == high --> (P2H == rise , SOC_target == high)

    p_now == low && p_trend == rising  && p_volatility == low  && SOC == low  --> (P2H == rise , SOC_target == high)
    p_now == low && p_trend == rising  && p_volatility == low  && SOC == high --> (P2H == rise , SOC_target == high)
    p_now == low && p_trend == rising  && p_volatility == high && SOC == low  --> (P2H == rise , SOC_target == high)
    p_now == low && p_trend == rising  && p_volatility == high && SOC == high --> (P2H == rise , SOC_target == high)

    p_now == med && p_trend == falling && p_volatility == low  && SOC == low  --> (P2H == keep  , SOC_target == low)
    p_now == med && p_trend == falling && p_volatility == low  && SOC == high --> (P2H == lower , SOC_target == low)
    p_now == med && p_trend == falling && p_volatility == high && SOC == low  --> (P2H == keep  , SOC_target == med)
    p_now == med && p_trend == falling && p_volatility == high && SOC == high --> (P2H == keep  , SOC_target == med)

    p_now == med && p_trend == rising  && p_volatility == low  && SOC == low  --> (P2H == rise  , SOC_target == high)
    p_now == med && p_trend == rising  && p_volatility == low  && SOC == high --> (P2H == rise  , SOC_target == high)
    p_now == med && p_trend == rising  && p_volatility == high && SOC == low  --> (P2H == rise  , SOC_target == high)
    p_now == med && p_trend == rising  && p_volatility == high && SOC == high --> (P2H == rise  , SOC_target == high)

    p_now == high && p_trend == falling && p_volatility == low  && SOC == low  --> (P2H == lower , SOC_target == low)
    p_now == high && p_trend == falling && p_volatility == low  && SOC == high --> (P2H == lower , SOC_target == low)
    p_now == high && p_trend == falling && p_volatility == high && SOC == low  --> (P2H == lower , SOC_target == low)
    p_now == high && p_trend == falling && p_volatility == high && SOC == high --> (P2H == lower , SOC_target == low)

    p_now == high && p_trend == rising  && p_volatility == low  && SOC == low  --> (P2H == keep  , SOC_target == med)
    p_now == high && p_trend == rising  && p_volatility == low  && SOC == high --> (P2H == rise  , SOC_target == high)
    p_now == high && p_trend == rising  && p_volatility == high && SOC == low  --> (P2H == keep  , SOC_target == high)
    p_now == high && p_trend == rising  && p_volatility == high && SOC == high --> (P2H == keep  , SOC_target == med)
  

end

# output of compilefis should be used in the controller because it's a lot faster and 
# doesn't need FuzzyLogic as a dependecy
compilefis(fis)