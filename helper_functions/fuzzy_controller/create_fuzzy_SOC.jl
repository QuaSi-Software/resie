using FuzzyLogic

fis = @mamfis function fuzzy_control_ems(p_now, p_trend)::{plr_max, chargemode}
    p_now := begin
        domain = -150:1100 
        cheap = TrapezoidalMF(-137, -136, 10, 75)
        average = TriangularMF(55, 80, 105)
        expensive = TrapezoidalMF(100, 120, 1000, 1001)
        
    end

    p_trend := begin
        domain = -300:300
        falling = TrapezoidalMF(-258, -257, -7, 1)
        stable = TriangularMF(-7, 0, 6)
        rising = TrapezoidalMF(-1, 6, 246, 247)
    end

    plr_max := begin
        domain = -1.0:2.0
        low = TriangularMF(-1.0, 0.0, 0.5)
        mid = TriangularMF(0.0, 0.5, 1.0)
        high = TriangularMF(0.5, 1.0, 2.0)
    end

    chargemode := begin
        domain = -1.0:101.0
        off = TriangularMF(-1.0, 0.0, 51.0)
        on = TriangularMF(49.0, 100.0, 101.0)
    end

    
    #p_volatility := begin
     #   domain = -10:300  
     #   stable = TriangularMF(-1, 0, 10)
     #   volatile = TrapezoidalMF(8, 15, 265, 266)
    #end

    #SOC_now := begin
     #   domain = -1.0:2.0
      #  empty = TrapezoidalMF(-0.5, 0.0, 0.15, 0.3)
      #  middle = TriangularMF(0.2, 0.5, 0.8)
      #  full = TrapezoidalMF(0.7, 0.85, 1.0, 1.5)
    #end


p_now	==	cheap	    &&	p_trend	==	falling	-->	(	plr_max	==	high	,	chargemode	==	off	)
p_now	==	cheap	    &&	p_trend	==	stable	-->	(	plr_max	==	high	,	chargemode	==	on	)
p_now	==	cheap	    &&	p_trend	==	rising	-->	(	plr_max	==	high	,	chargemode	==	on	)
p_now	==	average	    &&	p_trend	==	falling	-->	(	plr_max	==	mid	    ,	chargemode	==	off	)
p_now	==	average	    &&	p_trend	==	stable	-->	(	plr_max	==	mid	    ,	chargemode	==	on	)
p_now	==	average	    &&	p_trend	==	rising	-->	(	plr_max	==	mid	    ,	chargemode	==	on	)
p_now	==	expensive	&&	p_trend	==	falling	-->	(	plr_max	==	low	    ,	chargemode	==	off	)
p_now	==	expensive	&&	p_trend	==	stable	-->	(	plr_max	==	low	    ,	chargemode	==	off	)
p_now	==	expensive	&&	p_trend	==	rising	-->	(	plr_max	==	low	    ,	chargemode	==	on	)



end

# output of compilefis should be used in the controller because it's a lot faster and 
# doesn't need FuzzyLogic as a dependecy
compilefis(fis)
