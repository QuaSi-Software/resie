using LinearAlgebra

"""
I_1, I_2
V_full (P1)
V_2, Q_2 in exponential zone (P2)
V_3, Q_3 in exponential zone Q_3 = 2*Q_2 (P3)
V_4, Q_4 near end of nominal zone (P4)
V_cut, Q_full_1 end of curve for I_1 (P5)
Q_full_1, Q_full_2 capacities for I_1 and I_2 at end of curve
n_1, n_2 for n_2 > n_1 > 1
Q_full_1, Q_full_n_2 capacities at n_1 and n_2 cycles at reference Temperature at I_n and n_2 > n_1 > 1
V_full, V_full_n_2 Voltages at n_1 and n_2 cycles at reference Temperature at I_n and n_2 > n_1 > 1
V_nom_n_2, Q_nom_n_2 in nominal zone of curve for n_2 (P8)
T_1, T_2, T_ref for T_1 != T_2 != T_ref
Q_full_T_1, Q_full_T_2 capacities at T_1 and T_2 temperature at 1 cycle at I_T
V_full, V_full_T_2 voltages at T_1 and T_2 temperature at 1 cycle at I_T
V_nom_T_2, Q_nom_T_2 in nominal zone of curve for T_2 (P11)
r_i, 
V_6, Q_6 in nominal zone of curve I_2 only necessary if r_i not given

example Paper SP-LFP1000AHA:
calc_cell_values(300, 3.4, 3.318, 7.56, 3.301, 15.12, 3.172, 820, 2.0, 1070.4,
calc_cell_values(300, 3.4, 3.301, 15, 3.284, 30, 3.172, 820, 2.0, 1070.4,
calc_cell_values(100, 3.4, 3.346, 10, 3.332, 20, 3.216, 820, 2.0, 1087.5,
                 1000, 1061.5,
                 nothing, 3.141, 550,
                 500, 1000, 1026,
                 3000, 3.2, 3.15, 400, 806,
                 500, -20, 987,         
                 55, 3.31, 3.277, 696, 1035,
                 25)

paper reversed:
V_full = 3.4
I_1 = 100
Q_full_1 = 1090
AB_5 = 1.31835
V_cut = 2.0
V_3 = 3.3293
Q_3 = 20
V_2 = 3.3483
Q_2 = 10
V_4 = 3.2295
Q_4 = 800
I_2 = 1000
Q_full_2 = 1060
calc_cell_values(100, 3.4, 3.3483, 10, 3.3293, 20, 3.2295, 800, 2.0, 1090,
                 1000, 1060,
                 nothing, 3.141, 550,
                 500, 1000, 1027,
                 3000, 3.2, 3.15, 400, 809,
                 500, -20, 987,
                 55, 3.4, 3.2782, 757.57, 1102, 
                 25)
calc_cell_values(10, 3.4, 3.3483, 1.0, 3.3293, 2.0, 3.2295, 80, 2.0, 109,
                 100, 106.0,
                 nothing, 3.141, 55.0,
                 50, 1000, 102.7,
                 3000, 3.2, 3.15, 40.0, 80.9,
                 50, -20, 98.7,
                 55, 3.4, 3.2782, 75.757, 110.2, 
                 25)
"""
function calc_cell_values(I_1, V_full, V_2, Q_2, V_3, Q_3, V_4, Q_4, V_cut, Q_full_1,
                          I_2, Q_full_2,
                          r_i::Union{Number,Nothing}=nothing, V_6=0.0, Q_6=0.0,
                          I_n=0.0, n_1=0.0, Q_full_n_1=0.0,
                          n_2::Union{Number,Nothing}=nothing, V_full_n_2=0.0, V_nom_n_2=0.0, Q_nom_n_2=0.0, Q_full_n_2=0.0,
                          I_T=0.0, T_1=0.0, Q_full_T_1=0.0,
                          T_2::Union{Number,Nothing}=nothing, V_full_T_2=0.0, V_nom_T_2=0.0, Q_nom_T_2=0.0,Q_full_T_2=0.0,
                          T_ref=25)

    # basic values new
    if Q_2 == Q_full_1
        A=0
        B=1
    else
        B = -1 / Q_2 * log((V_full - V_3) / (V_full - V_2) - 1)
        A = (V_full - V_3) / (1 - exp(-B * Q_3))
    end

    AB_4 = V_full - V_4 - A * (1 - exp(-B * Q_4)) # exp(-B*Q_4) -> 0
    AB_5 = V_full - V_cut - A * (1 - exp(-B * Q_full_1)) # exp(-B*Q_full) -> 0

    m = (1 - AB_5 / AB_4) / (1 - (AB_5 * Q_4) / (AB_4 * Q_full_1)) * (Q_4 / Q_full_1)
    #TODO test m = 1
    K = AB_5 * ((m * Q_full_1 / Q_full_1) - 1)
    if isnothing(r_i)
        r_i = 1 / (I_1 - I_2) * (V_6 + K * (m * Q_full_2 / (m * Q_full_2 - Q_6)) - A * exp(-B * Q_6) -
                                 V_3 - K * (m * Q_full_1 / (m * Q_full_1 - Q_3)) + A * exp(-B * Q_3))
    end
    V_0 = V_full + K + r_i * I_1 - A
    alpha = log(Q_full_2 / Q_full_1) / log(I_2 / I_1)
    I_ref = I_1
    Q_ref = Q_full_1

    # basic values old
    # A = V_full - V_exp
    # B = 3/Q_exp
    # K = ((V_full - V_nom + A*(exp(-B*Q_nom) - 1))*(Q_full_1 - Q_nom))/(Q_nom)
    # V_0 = V_full + K + r_i*I_1 - A

    # for aging if the needed values are provided
    k_qn = [0.0, 0.0]
    k_n = [0.0, 0.0, 0.0, 0.0]
    if !isnothing(n_2)
        N = [(n_1 - 1) * n_1/2 (n_1 - 1) * n_1 * (2 * n_1 - 1)/2
             (n_2 - 1) * n_2/2 (n_2 - 1) * n_2 * (2 * n_2 - 1)/2]
        q = [Q_full_n_1 / (Q_ref * (I_n / I_ref)^alpha) - 1
             Q_full_n_2 / (Q_ref * (I_n / I_ref)^alpha) - 1]
        k_qn = N \ q

        Q_n = Q_ref * (I_n / I_ref)^alpha *
              (1 + k_qn[1] * (n_2 - 1) * n_2 / 2 + k_qn[2] * (n_2 - 1) * n_2 * (2 * n_2 - 1) / 2)
        b_n = [V_full - (V_0 - r_i * I_1 - K + A)
               V_full_n_2 - (V_0 - r_i * I_n - K + A)
               V_nom_n_2 - (V_0 - r_i * I_n - K * (m * Q_n / (m * Q_n - Q_nom_n_2)) + A * exp(-B * Q_nom_n_2))
               V_cut - (V_0 - r_i * I_n - K * (m * Q_n / (m * Q_n - Q_full_n_2)) + A * exp(-B * Q_full_n_2))]
        G_n = [V_0 -r_i*I_1 -K A
               V_0*(n_2 - 1) -r_i*I_n*(n_2-1) -K*(n_2 - 1) A*(n_2 - 1)
               V_0*(n_2 - 1) -r_i*I_n*(n_2-1) -K*(n_2-1)*(m * Q_n/(m * Q_n - Q_nom_n_2)) A*(n_2-1)*exp(-B * Q_nom_n_2)
               V_0*(n_2 - 1) -r_i*I_n*(n_2-1) -K*(n_2-1)*(m * Q_n/(m * Q_n - Q_full_n_2)) A*(n_2-1)*exp(-B * Q_full_n_2)]
        try k_n = G_n \ b_n
        catch
            try k_n = svd(G_n) \ b_n
            catch 
                k_n = pinv(G_n) * b_n
            end
        end
    end

    k_qT = [0.0, 0.0]
    k_T = [0.0, 0.0, 0.0, 0.0]
    if !isnothing(T_2)
        T = [T_1-T_ref (T_1 - T_ref)^2
             T_2-T_ref (T_2 - T_ref)^2]
        q = [Q_full_T_1 / (Q_ref * (I_T / I_ref)^alpha) - 1
             Q_full_T_2 / (Q_ref * (I_T / I_ref)^alpha) - 1]
        k_qT = T \ q

        Q_T = Q_ref * (I_T / I_ref)^alpha *
              (1 + k_qT[1] * (T_2 - T_ref) + k_qT[2] * (T_2 - T_ref)^2)
        b_T = [V_full - (V_0 - r_i * I_1 - K + A)
               V_full_T_2 - (V_0 - r_i * I_T - K + A)
               V_nom_T_2 - (V_0 - r_i * I_T - K * (m * Q_T / (m * Q_T - Q_nom_T_2)) + A * exp(-B * Q_nom_T_2))
               V_cut - (V_0 - r_i * I_T - K * (m * Q_T / (m * Q_T - Q_full_T_2)) + A * exp(-B * Q_full_T_2))]
        G_T = [V_0 -r_i*I_1 -K A
               V_0*(T_2 - T_ref) -r_i*I_T*(T_2-T_ref) -K*(T_2 - T_ref) A*(T_2 - T_ref)
               V_0*(T_2 - T_ref) -r_i*I_T*(T_2-T_ref) -K*(T_2-T_ref)*(m * Q_T/(m * Q_T - Q_nom_T_2)) A*(T_2-T_ref)*exp(-B * Q_nom_T_2)
               V_0*(T_2 - T_ref) -r_i*I_T*(T_2-T_ref) -K*(T_2-T_ref)*(m * Q_T/(m * Q_T - Q_full_T_2)) A*(T_2-T_ref)*exp(-B * Q_full_T_2)]
        try k_T = G_T \ b_T
        catch
            try k_T = svd(G_T) \ b_T
            catch 
                k_T = pinv(G_T) * b_T
            end
        end
    end

    return V_0, K, A, B, r_i, m, alpha, k_qn, k_qT, k_n, k_T, I_ref, T_ref
end

# NiMH
# calc_cell_values(0.03,1.4,1.2238,0.03057,1.2069,0.06114,1.1296,0.141093,1.0121,0.15627,
#                  0.15,0.0993345,
#                  0.8,0,0,
#                  0.03,2,0.15626,
#                  nothing,1.399,1.2033,0.075743,0.15625,
#                  0.03,0,0.134967,
#                  50,1.305,1.2019,0.0757,0.15504,
#                  23)

# NMC
calc_cell_values(
1.500,
4.137,
3.688,
1.205,
3.422,
2.409,
3.296,
2.805,
2.760,
3.034,

3.000,
3.005,

0.030,
0.000,
0.000,

1.500,
200.0,
1.455,

2.000,
4.136,
3.507,
1.950,
3.033,

1.500,
10.00,
2.782,

55.00,
4.193,
3.546,
2.137,
3.122,

25.00
)