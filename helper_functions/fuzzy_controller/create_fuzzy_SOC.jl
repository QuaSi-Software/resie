using FuzzyLogic

fis = @mamfis function fuzzy_control_ems(p_now, p_trend, p_volatility, SOC_now)::{P2H, SOC_target}
    p_now := begin
        domain = -150:1100 # TODO change to min max 
        cheap = TrapezoidalMF(-137, -136, 10, 75)
        average = TriangularMF(55, 80, 105)
        expensive = TrapezoidalMF(100, 120, 1000, 1001)
        
    end

    p_trend := begin
        domain = -300:300 # TODO change to min max 
        falling = TrapezoidalMF(-258, -257, -7, 1)
        rising = TrapezoidalMF(-1, 6, 246, 247)
    end

    p_volatility := begin
        domain = -10:300 # TODO change to min max 
        stable = TriangularMF(-1, 0, 10)
        volatile = TrapezoidalMF(8, 15, 265, 266)
    end

    SOC_now := begin
        domain = -1.0:2.0
        empty = TrapezoidalMF(-0.5, 0.0, 0.15, 0.3)
        middle = TriangularMF(0.2, 0.5, 0.8)
        full = TrapezoidalMF(0.7, 0.85, 1.0, 1.5)
    end

    P2H := begin
        domain = -2.0:2.0
        lower = TriangularMF(-2.0, -1.0, -0.2)
        keep = TriangularMF(-0.3, 0.0, 0.3)
        rise = TriangularMF(0.2, 1.0, 2.0)
    end

    SOC_target := begin
        domain = -2.0:2.0
        low = TriangularMF(-1.0, -0.5, 0.5)
        med = TriangularMF(0.0, 0.5, 1.0)
        high = TriangularMF(0.5, 1.0, 1.5)
    end

p_now	==	cheap	    &&	p_trend	==	falling	&&	p_volatility	==	stable	    &&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	med	)
p_now	==	cheap	    &&	p_trend	==	falling	&&	p_volatility	==	stable	    &&	SOC_now	==	middle	-->	(	P2H	==	keep	,	SOC_target	==	med	)
p_now	==	cheap	    &&	p_trend	==	falling	&&	p_volatility	==	stable	    &&	SOC_now	==	full	-->	(	P2H	==	keep	,	SOC_target	==	high)
p_now	==	cheap	    &&	p_trend	==	falling	&&	p_volatility	==	volatile	&&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	cheap	    &&	p_trend	==	falling	&&	p_volatility	==	volatile	&&	SOC_now	==	middle	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	cheap	    &&	p_trend	==	falling	&&	p_volatility	==	volatile	&&	SOC_now	==	full	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	cheap	    &&	p_trend	==	rising	&&	p_volatility	==	stable	    &&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	cheap	    &&	p_trend	==	rising	&&	p_volatility	==	stable	    &&	SOC_now	==	middle	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	cheap	    &&	p_trend	==	rising	&&	p_volatility	==	stable	    &&	SOC_now	==	full	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	cheap	    &&	p_trend	==	rising	&&	p_volatility	==	volatile	&&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	cheap	    &&	p_trend	==	rising	&&	p_volatility	==	volatile	&&	SOC_now	==	middle	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	cheap	    &&	p_trend	==	rising	&&	p_volatility	==	volatile	&&	SOC_now	==	full	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	average	    &&	p_trend	==	falling	&&	p_volatility	==	stable	    &&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	med	)
p_now	==	average	    &&	p_trend	==	falling	&&	p_volatility	==	stable	    &&	SOC_now	==	middle	-->	(	P2H	==	lower	,	SOC_target	==	low	)
p_now	==	average	    &&	p_trend	==	falling	&&	p_volatility	==	stable	    &&	SOC_now	==	full	-->	(	P2H	==	lower	,	SOC_target	==	low	)
p_now	==	average	    &&	p_trend	==	falling	&&	p_volatility	==	volatile	&&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	med	)
p_now	==	average	    &&	p_trend	==	falling	&&	p_volatility	==	volatile	&&	SOC_now	==	middle	-->	(	P2H	==	keep	,	SOC_target	==	med	)
p_now	==	average	    &&	p_trend	==	falling	&&	p_volatility	==	volatile	&&	SOC_now	==	full	-->	(	P2H	==	keep	,	SOC_target	==	high)
p_now	==	average	    &&	p_trend	==	rising	&&	p_volatility	==	stable	    &&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	med	)
p_now	==	average	    &&	p_trend	==	rising	&&	p_volatility	==	stable	    &&	SOC_now	==	middle	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	average	    &&	p_trend	==	rising	&&	p_volatility	==	stable	    &&	SOC_now	==	full	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	average	    &&	p_trend	==	rising	&&	p_volatility	==	volatile	&&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	med	)
p_now	==	average	    &&	p_trend	==	rising	&&	p_volatility	==	volatile	&&	SOC_now	==	middle	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	average	    &&	p_trend	==	rising	&&	p_volatility	==	volatile	&&	SOC_now	==	full	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	expensive	&&	p_trend	==	falling	&&	p_volatility	==	stable	    &&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	med	)
p_now	==	expensive	&&	p_trend	==	falling	&&	p_volatility	==	stable	    &&	SOC_now	==	middle	-->	(	P2H	==	lower	,	SOC_target	==	low	)
p_now	==	expensive	&&	p_trend	==	falling	&&	p_volatility	==	stable	    &&	SOC_now	==	full	-->	(	P2H	==	lower	,	SOC_target	==	low	)
p_now	==	expensive	&&	p_trend	==	falling	&&	p_volatility	==	volatile	&&	SOC_now	==	empty	-->	(	P2H	==	keep	,	SOC_target	==	med	)
p_now	==	expensive	&&	p_trend	==	falling	&&	p_volatility	==	volatile	&&	SOC_now	==	middle	-->	(	P2H	==	keep	,	SOC_target	==	low	)
p_now	==	expensive	&&	p_trend	==	falling	&&	p_volatility	==	volatile	&&	SOC_now	==	full	-->	(	P2H	==	keep	,	SOC_target	==	low	)
p_now	==	expensive	&&	p_trend	==	rising	&&	p_volatility	==	stable	    &&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	med	)
p_now	==	expensive	&&	p_trend	==	rising	&&	p_volatility	==	stable	    &&	SOC_now	==	middle	-->	(	P2H	==	rise	,	SOC_target	==	med	)
p_now	==	expensive	&&	p_trend	==	rising	&&	p_volatility	==	stable	    &&	SOC_now	==	full	-->	(	P2H	==	keep	,	SOC_target	==	high)
p_now	==	expensive	&&	p_trend	==	rising	&&	p_volatility	==	volatile	&&	SOC_now	==	empty	-->	(	P2H	==	rise	,	SOC_target	==	high)
p_now	==	expensive	&&	p_trend	==	rising	&&	p_volatility	==	volatile	&&	SOC_now	==	middle	-->	(	P2H	==	keep	,	SOC_target	==	med	)
p_now	==	expensive	&&	p_trend	==	rising	&&	p_volatility	==	volatile	&&	SOC_now	==	full	-->	(	P2H	==	keep	,	SOC_target	==	med	)



end

# output of compilefis should be used in the controller because it's a lot faster and 
# doesn't need FuzzyLogic as a dependecy
compilefis(fis)
