using FuzzyLogic

fis = @mamfis function fuzzy_control(p_now, p_trend, p_volatility, SOC_now)::{P2H, SOC_target}
    p_now := begin
        domain = -136:1000 # TODO change to min max 
        low = TrapezoidalMF(-200, -136, 55, 80)
        med = TriangularMF(55, 80, 100)
        high = TrapezoidalMF(80, 100, 1000, 1500)
        
    end

    p_trend := begin
        domain = -260:250 # TODO change to min max 
        falling = TrapezoidalMF(-270, -260, -7, -1)
        stable = TriangularMF(-7, 1, 6)
        rising = TrapezoidalMF(1, 6, 250, 270)
    end

    p_volatility := begin
        domain = 0:265 # TODO change to min max 
        low = TrapezoidalMF(-10, 0, 4, 8)
        med = TriangularMF(4, 8, 15)
        high = TrapezoidalMF(8, 15, 265, 275)
    end

    SOC_now := begin
        domain = 0:1
        low = TriangularMF(-0.5, 0.0, 0.5)
        med = TriangularMF(0.0, 0.5, 1.0)
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

    p_now == low    &&  p_trend	== falling	&&	p_volatility == low     &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == low    &&	p_trend	== falling	&&	p_volatility == low     &&	SOC_now	== med	  -->	(P2H == lower	,	SOC_target == med)
    p_now == low    &&	p_trend	== falling	&&	p_volatility == low     &&	SOC_now	== high	  -->	(P2H == lower	,	SOC_target == high)
    p_now == low    &&	p_trend	== falling	&&	p_volatility == med     &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== falling	&&	p_volatility == med     &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== falling	&&	p_volatility == med     &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== falling	&&	p_volatility == high    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== falling	&&	p_volatility == high    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== falling	&&	p_volatility == high    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== stable	&&	p_volatility == low	    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== stable	&&	p_volatility == low	    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== stable	&&	p_volatility == low	    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== stable	&&	p_volatility == med	    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== stable	&&	p_volatility == med	    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== stable	&&	p_volatility == med	    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== stable	&&	p_volatility == high    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== stable	&&	p_volatility == high    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== stable	&&	p_volatility == high    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== rising	&&	p_volatility == low	    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== rising	&&	p_volatility == low	    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== rising	&&	p_volatility == low	    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== rising	&&	p_volatility == med	    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== rising	&&	p_volatility == med	    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== rising	&&	p_volatility == med	    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== rising	&&	p_volatility == high    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== rising	&&	p_volatility == high    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == low    &&	p_trend	== rising	&&	p_volatility == high    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== falling	&&	p_volatility == low	    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == med    &&	p_trend	== falling	&&	p_volatility == low	    &&	SOC_now	== med	  -->	(P2H == lower	,	SOC_target == med)
    p_now == med    &&	p_trend	== falling	&&	p_volatility == low	    &&	SOC_now	== high	  -->	(P2H == lower	,	SOC_target == high)
    p_now == med    &&	p_trend	== falling	&&	p_volatility == med	    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == med    &&	p_trend	== falling	&&	p_volatility == med	    &&	SOC_now	== med	  -->	(P2H == lower	,	SOC_target == med)
    p_now == med    &&	p_trend	== falling	&&	p_volatility == med	    &&	SOC_now	== high	  -->	(P2H == lower	,	SOC_target == high)
    p_now == med    &&	p_trend	== falling	&&	p_volatility == high    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == med    &&	p_trend	== falling	&&	p_volatility == high    &&	SOC_now	== med	  -->	(P2H == keep	    ,	SOC_target == med)
    p_now == med    &&	p_trend	== falling	&&	p_volatility == high    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== stable	&&	p_volatility == low	    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == med    &&	p_trend	== stable	&&	p_volatility == low	    &&	SOC_now	== med	  -->	(P2H == keep	    ,	SOC_target == med)
    p_now == med    &&	p_trend	== stable	&&	p_volatility == low	    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== stable	&&	p_volatility == med	    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == med    &&	p_trend	== stable	&&	p_volatility == med	    &&	SOC_now	== med    -->	(P2H == keep	    ,	SOC_target == med)
    p_now == med    &&	p_trend	== stable	&&	p_volatility == med	    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== stable	&&	p_volatility == high    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == med    &&	p_trend	== stable	&&	p_volatility == high    &&	SOC_now	== med	  -->	(P2H == keep	    ,	SOC_target == med)
    p_now == med    &&	p_trend	== stable	&&	p_volatility == high    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== rising	&&	p_volatility == low	    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== rising	&&	p_volatility == low	    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== rising	&&	p_volatility == low	    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== rising	&&	p_volatility == med	    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== rising	&&	p_volatility == med	    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== rising	&&	p_volatility == med	    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== rising	&&	p_volatility == high    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== rising	&&	p_volatility == high    &&	SOC_now	== med	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == med    &&	p_trend	== rising	&&	p_volatility == high    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == high   &&	p_trend	== falling	&&	p_volatility == low	    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == high   &&	p_trend	== falling	&&	p_volatility == low	    &&	SOC_now	== med	  -->	(P2H == lower	,	SOC_target == low)
    p_now == high   &&	p_trend	== falling	&&	p_volatility == low	    &&	SOC_now	== high	  -->	(P2H == lower	,	SOC_target == low)
    p_now == high   &&	p_trend	== falling	&&	p_volatility == med	    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == high   &&	p_trend	== falling	&&	p_volatility == med	    &&	SOC_now	== med	  -->	(P2H == lower	,	SOC_target == low)
    p_now == high   &&	p_trend	== falling	&&	p_volatility == med	    &&	SOC_now	== high	  -->	(P2H == lower	,	SOC_target == low)
    p_now == high   &&	p_trend	== falling	&&	p_volatility == high    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == high   &&	p_trend	== falling	&&	p_volatility == high    &&	SOC_now	== med	  -->	(P2H == lower	,	SOC_target == low)
    p_now == high   &&	p_trend	== falling	&&	p_volatility == high    &&	SOC_now	== high	  -->	(P2H == lower	,	SOC_target == low)
    p_now == high   &&	p_trend	== stable	&&	p_volatility == low	    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == high   &&	p_trend	== stable	&&	p_volatility == low	    &&	SOC_now	== med	  -->	(P2H == keep	    ,	SOC_target == med)
    p_now == high   &&	p_trend	== stable	&&	p_volatility == low	    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == high   &&	p_trend	== stable	&&	p_volatility == med	    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == high   &&	p_trend	== stable	&&	p_volatility == med	    &&	SOC_now	== med	  -->	(P2H == keep	    ,	SOC_target == med)
    p_now == high   &&	p_trend	== stable	&&	p_volatility == med	    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == high   &&	p_trend	== stable	&&	p_volatility == high    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == high   &&	p_trend	== stable	&&	p_volatility == high    &&	SOC_now	== med	  -->	(P2H == keep	    ,	SOC_target == med)
    p_now == high   &&	p_trend	== stable	&&	p_volatility == high    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == high   &&	p_trend	== rising	&&	p_volatility == low	    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == med)
    p_now == high   &&	p_trend	== rising	&&	p_volatility == low	    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == high   &&	p_trend	== rising	&&	p_volatility == low	    &&	SOC_now	== high	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == high   &&	p_trend	== rising	&&	p_volatility == med	    &&	SOC_now	== low	  -->	(P2H == rise	    ,	SOC_target == med)
    p_now == high   &&	p_trend	== rising	&&	p_volatility == med	    &&	SOC_now	== med	  -->	(P2H == rise	    ,	SOC_target == high)
    p_now == high   &&	p_trend	== rising	&&	p_volatility == med	    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high)
    p_now == high   &&	p_trend	== rising	&&	p_volatility == high    &&	SOC_now	== low	  -->	(P2H == keep	    ,	SOC_target == low)
    p_now == high   &&	p_trend	== rising	&&	p_volatility == high    &&	SOC_now	== med	  -->	(P2H == keep	    ,	SOC_target == med)
    p_now == high   &&	p_trend	== rising	&&	p_volatility == high    &&	SOC_now	== high	  -->	(P2H == keep	    ,	SOC_target == high) 
    

end

# output of compilefis should be used in the controller because it's a lot faster and 
# doesn't need FuzzyLogic as a dependecy
compilefis(fis)