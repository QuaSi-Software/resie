using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles
import Plots as p
import Roots
using CSV
using DataFrames
p.plotlyjs()

include("../test/test_util.jl")

function get_config_battery_detailed()
    return Dict{String,Any}( 
        "TST_GRI_EL_01"=> Dict{String,Any}(
            "type"=> "GridConnection",
            "medium"=> "m_e_ac_230v",
            "output_refs"=> ["TST_BUS_EL_01"],
            "is_source"=> true,
        ),
        "TST_GRI_EL_01_out"=> Dict{String,Any}(
            "type"=> "GridConnection",
            "medium"=> "m_e_ac_230v",
            "output_refs"=> [],
            "is_source"=> false,
        ),
        "TST_BUS_EL_01"=> Dict{String,Any}(
            "type"=> "Bus",
            "medium"=> "m_e_ac_230v",
            "connections"=> Dict{String,Any}(
                "input_order"=> [
                    "SP-LFP1000AHA",
                    "TST_GRI_EL_01",
                ],
                "output_order"=> [
                    "SP-LFP1000AHA",
                    "TST_GRI_EL_01_out",
                ],
                "energy_flow"=>[
                    [0,1],
                    [1,0],
                ]
            )
        ),
        "SP-LFP1000AHA"=> Dict{String,Any}(
            "type"=> "Battery",
            "output_refs"=> ["TST_BUS_EL_01"],
            "model_type"=> "detailed",
            "capacity"=> 3200,
            "load"=> 3200,
            "self_discharge_rate"=> 0.05,
            "max_charge_C_rate"=> 1,
            "max_discharge_C_rate"=> 1,

            "V_n_bat"=> 4,
            "SOC_min"=> 10,
            "SOC_max"=> 100,
         
            "V_n"=> 4,
            "V_cell_min"=> 2.0,
            "V_cell_max"=> 3.4,
            "capacity_cell_Ah"=> 0.15627,
        ),
    )  
end

function plot_curves(params, V_min, V_max, nominal_cell_capacity, c_rates, cycles, c_ref_cyc, Temps, c_ref_temp)
    components_config = get_config_battery_detailed()
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    setup_mock_run!(components, simulation_parameters)

    battery = components["SP-LFP1000AHA"]
    battery.V_0, battery.K, battery.A, battery.B, battery.r_i, battery.m, battery.alpha, 
    battery.k_qn, battery.k_qT, battery.k_n, battery.k_T, battery.I_ref, battery.T_ref = params
    battery.extracted_charge = 0
    battery.V_cell_min = V_min
    battery.V_cell_max = V_max
    battery.capacity_cell_Ah = nominal_cell_capacity

    # plot parameters
    minorticks=10
    size=(1500, 950)
    titlefontsize=24
    guidefontsize=18
    tickfontsize=18
    legendfontsize=18
    gridlinewidth=1
    margin=15p.mm

    # discharge curves of the cell chemistry at different discharge rates
    Q_max = battery.capacity_cell_Ah
    try Q_max = Roots.find_zero(Q -> f_V_cell(Q, minimum(c_rates) * battery.capacity_cell_Ah, battery, 1, maximum(Temps)) - battery.V_cell_min,
                                battery.capacity_cell_Ah) * 1.001
    catch    
    end
    Q_max = max(Q_max, battery.capacity_cell_Ah)
    x = 0:(Q_max / 500):Q_max
    y1 = fill(NaN, length(x))
    y2 = fill(NaN, length(x))
    y3 = fill(NaN, length(x))
    for (idx, Q) in enumerate(x)
        y1[idx] = f_V_cell(Q, c_rates[1] * battery.capacity_cell_Ah, battery, 1, 25)
        if y1[idx] <= battery.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y2[idx] = f_V_cell(Q, c_rates[2] * battery.capacity_cell_Ah, battery, 1, 25)
        if y2[idx] <= battery.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y3[idx] = f_V_cell(Q, c_rates[3] * battery.capacity_cell_Ah, battery, 1, 25)
        if y3[idx] <= battery.V_cell_min
            break
        end
    end
    fig = p.plot(x, [y1, y2, y3], 
                   labels=["$(c_rates[1])C" "$(c_rates[2])C" "$(c_rates[3])C"], 
                   lw=3)
    p.plot!(fig, ; title="Battery discharge curve at different C-rates",
                xlabel="Removed Cell Charge / Ah",
                ylabel="Cell Voltage / V",
                ylims=(battery.V_cell_min, battery.V_cell_max),
                xlims=(0, Q_max),
                minorticks=minorticks,
                size=size,
                titlefontsize=titlefontsize,
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                grid=true,
                minorgrid=true,
                gridlinewidth=gridlinewidth,
                margin=margin)
    display(fig)

    # discharge curves of the cell chemistry at different Temperatures
    y1 = fill(NaN, length(x))
    y2 = fill(NaN, length(x))
    y3 = fill(NaN, length(x))
    for (idx, Q) in enumerate(x)
        y1[idx] = f_V_cell(Q, c_ref_temp * battery.capacity_cell_Ah, battery, 1, Temps[1])
        if y1[idx] <= battery.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y2[idx] = f_V_cell(Q, c_ref_temp * battery.capacity_cell_Ah, battery, 1, Temps[2])
        if y2[idx] <= battery.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y3[idx] = f_V_cell(Q, c_ref_temp * battery.capacity_cell_Ah, battery, 1, Temps[3])
        if y3[idx] <= battery.V_cell_min
            break
        end
    end
    fig = p.plot(x, [y1, y2, y3], 
                   labels=["$(Temps[1]) °C" "$(Temps[2]) °C" "$(Temps[3]) °C"], 
                   lw=3)
    p.plot!(fig, ; title="Battery discharge curve at different Temperatures",
                xlabel="Removed Cell Charge / Ah",
                ylabel="Cell Voltage / V",
                ylims=(battery.V_cell_min, battery.V_cell_max),
                xlims=(0, Q_max),
                minorticks=minorticks,
                size=size,
                titlefontsize=titlefontsize,
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                grid=true,
                minorgrid=true,
                gridlinewidth=gridlinewidth,
                margin=margin)
    display(fig)

    # discharge curves of the cell chemistry after different number of cycles
    y1 = fill(NaN, length(x))
    y2 = fill(NaN, length(x))
    y3 = fill(NaN, length(x))
    for (idx, Q) in enumerate(x)
        y1[idx] = f_V_cell(Q, c_ref_cyc * battery.capacity_cell_Ah, battery, cycles[1], 25)
        if y1[idx] <= battery.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y2[idx] = f_V_cell(Q, c_ref_cyc * battery.capacity_cell_Ah, battery, cycles[2], 25)
        if y2[idx] <= battery.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y3[idx] = f_V_cell(Q, c_ref_cyc * battery.capacity_cell_Ah, battery, cycles[3], 25)
        if y3[idx] <= battery.V_cell_min
            break
        end
    end
    fig = p.plot(x, [y1, y2, y3], 
                   labels=["$(cycles[1]) cycles" "$(cycles[2]) cycles" "$(cycles[3]) cycles"], 
                   lw=3)
    p.plot!(fig, ; title="Battery discharge curve after different number of cycles",
                xlabel="Removed Cell Charge / Ah",
                ylabel="Cell Voltage / V",
                ylims=(battery.V_cell_min, battery.V_cell_max),
                xlims=(0, Q_max),
                minorticks=minorticks,
                size=size,
                titlefontsize=titlefontsize,
                guidefontsize=guidefontsize,
                tickfontsize=tickfontsize,
                legendfontsize=legendfontsize,
                grid=true,
                minorgrid=true,
                gridlinewidth=gridlinewidth,
                margin=margin)
    display(fig)
end


# V_0, K, A, B, r_i, m, alpha, 
# k_qn, k_qT, 
# k_n, 
# k_T, i_ref, T_ref

#LIIon LFP
# parameters = 3.36964, 0.03546, 0.08165, 0.1003, 0.00016, 1.0269, -0.01212, 
#              [-1.27571e-7, 1.22095e-11], [1.32729e-3, -7.9763e-6], 
#              [9.71249e-6, 7.51635e-4, -8.59363e-5, -2.92489e-4], [1.05135e-3, 1.83721e-2, -7.72438e-3, -4.31833e-2], 100, 25
# V_min = 2.0
# V_max = 3.4
# Q_full_1 = 1090

# parameters = 3.36964, 0.03546, 0.08165, 10.0102, 0.0157, 1.0269, -0.01212, 
#              [-1.27571e-7, 1.22095e-11], [1.32729e-3, -7.9763e-6], 
#              [9.71249e-6, 7.51635e-4, -8.59363e-5, -2.92489e-4], [1.05135e-3, 1.83721e-2, -7.72438e-3, -4.31833e-2], 1, 25
# parameters = 3.3694, 0.03543, 0.08174, 10.0102, 0.015747, 1.02687, -0.01212, 
#              [-9.079645e-8, 6.1242e-12], [1.2999e-3, -8.96904e-6], 
#              [9.0743e-6, 7.2532e-4, -3.51154e-5, -2.49544e-4], [1.84665e-4, -3.3333e-2, 8.474e-4, -1.3666e-2], 10, 25
Q_full_1 = 10.9
parameters = 3.3694, 0.03543, 0.08174, 1.00102, 0.0015747, 1.02687 -0.0121, 
             [-9.079639e-8, 6.1242e-12], [1.2999e-3, -8.96904e-6], 
             [9.0743e-6, 7.2532e-4, -3.51154e-5, -2.49544e-4], [1.84665e-4, -3.3333e-2, 8.474e-4, -1.3666e-2], 10, 25
Q_full_1 = 109
V_min = 2.0
V_max = 3.4

# NiMH
# parameters = 1.2434404981985387, 0.014333404664326461, 0.1948929064657876, 76.68649301841862, 0.8, 1.0742631476488769, -0.28152209798835337, 
#              [0.0, 0.0], [0.003066499236889949, -0.00012437101770780085], 
#              [0.0, 0.0, 0.0, 0.0], [0.00031992282665392217, -6.1749218003040444e-6, 0.0019188297799150977, -0.07791049727450391], 0.03, 23
# V_min = 1.0121
# V_max = 1.4
# Q_full_1 = 0.15627
   
#LiIon NMC
# parameters = 4.322089864650338, 0.22465804646851967, 0.08456818181818074, 47.320760004338, 0.03, 1.1738258400234733, -0.013856094648485657, 
#              [0.0, 0.0], [0.0006298981908737985, 1.1230742938060141e-5], 
#              [0.0, 0.0, 0.0, 0.0], [-0.0027225713491162828, 2.8346405221964717e-5, -0.00775726863401747, -0.4310834811137324], 1.5, 25
# V_min = 2.76
# V_max = 4.137
# Q_full_1 = 3.034

c_rates = [0.1, 0.5, 1]
cycles = [1, 500, 1000]
c_rate_for_cycles = 0.5
Temps = [40, 25, -10]
c_rate_for_Temps = 0.5
plot_curves(parameters, V_min, V_max, Q_full_1, c_rates, cycles, c_rate_for_cycles, Temps, c_rate_for_Temps)

# datasheet = CSV.read("profiles/SP_LFP1000AHA/0.1C.csv", DataFrame, header=["Q","V"])
# p.scatter!(datasheet.Q, datasheet.V, label="datasheet 0.1C", ms=2 , markerstrokeopacity=0)
# datasheet = CSV.read("profiles/SP_LFP1000AHA/0.3C.csv", DataFrame, header=["Q","V"])
# p.scatter!(datasheet.Q, datasheet.V, label="datasheet 0.3C", ms=2 , markerstrokeopacity=0)
# datasheet = CSV.read("profiles/SP_LFP1000AHA/1C.csv", DataFrame, header=["Q","V"])
# p.scatter!(datasheet.Q, datasheet.V, label="datasheet 1C", ms=2 , markerstrokeopacity=0)
