using Debugger
using Test
using Resie
using Resie.EnergySystems
using Resie.Profiles
import Plots as p
using CSV
using DataFrames
p.plotlyjs()

include("../test_util.jl")

function get_config_battery_detailed()
    return Dict{String,Any}( 
        "TST_DEM_01"=> Dict{String,Any}(
            "type"=> "Demand",
            "medium"=> "m_e_ac_230v",
            "output_refs"=> [],
            "energy_profile_file_path"=> "./profiles/tests/demand_electricity.prf",
            "scale"=> 1500,
           ),
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
                    "TST_PV_01",
                    "BYD_HVM_8.3",
                    "SP-LFP1000AHA",
                    "TST_GRI_EL_01",
                ],
                "output_order"=> [
                    "TST_DEM_01",
                    "BYD_HVM_8.3",
                    "SP-LFP1000AHA",
                    "TST_GRI_EL_01_out",
                ],
                "energy_flow"=>[
                    [1,1,1,1],
                    [1,0,0,0],
                    [1,0,0,0],
                    [1,0,0,0],
                ]
            )
        ),
        "TST_PV_01"=> Dict{String,Any}(
            "type"=> "PVPlant",
            "output_refs"=> ["TST_BUS_EL_01"],
            "energy_profile_file_path"=> "./profiles/examples/district/PV_Stuttgart_30_south_Wh_per_square_meter.prf",
            "scale"=> 1000,
        ),
        "BYD_HVM_8.3"=> Dict{String,Any}(
            "type"=> "Battery",
            "output_refs"=> ["TST_BUS_EL_01"],
            "model_type"=> "with_aging",
            "capacity"=> 8.590881628999776,
            "load"=> 8.590881628999776,
            "self_discharge_rate"=> 0.03,
            "max_charge_power"=> 2.25,
            "max_discharge_power"=> 2.25,

            "V_n_bat"=> 3.6,
            "SOC_min"=> 10,
            "SOC_max"=> 100,
         
            "V_n"=> 3.6,
            "r_i"=> 0.002,
            "V_0"=> 4.13215,
            "K"=> 0.08125,
            "A"=> 0.05,
            "B"=> 75,
            "V_cell_min"=> 2.706,
            "V_cell_max"=> 4.1,
            "capacity_cell_Ah"=> 2.25,
        ),        
        "SP-LFP1000AHA"=> Dict{String,Any}(
            "type"=> "Battery",
            "output_refs"=> ["TST_BUS_EL_01"],
            "model_type"=> "with_aging",
            "capacity"=> 3200,
            "load"=> 3200,
            "self_discharge_rate"=> 0.05,
            "max_charge_power"=> 3000,
            "max_discharge_power"=> 3000,

            "V_n_bat"=> 3.2,
            "SOC_min"=> 10,
            "SOC_max"=> 100,
         
            "V_n"=> 3.2,
            "V_cell_min"=> 2.0,
            "V_cell_max"=> 3.4,
            "capacity_cell_Ah"=> 1090,
        ),
    )  
end

function calc_and_plot_discharge_charge_curve(c_rate, delta_t, aging_params, n=1, T=25)
    components_config = get_config_battery_detailed()
    simulation_parameters = get_default_sim_params()

    components = Resie.load_components(components_config, simulation_parameters)
    setup_mock_run!(components, simulation_parameters)
    simulation_parameters["time_step_seconds"] = delta_t * 3600

    battery = components["SP-LFP1000AHA"]
    battery.V_0, battery.K, battery.A, battery.B, battery.r_i, battery.m, battery.alpha, 
    battery.k_qn, battery.k_qT, battery.k_n, battery.k_T, battery.i_ref, battery.T_ref = aging_params
    battery.current_charge = battery.capacity_cell_Ah

    # battery.capacity_cell_Ah_cutoff = battery.capacity_cell_Ah

    I = battery.capacity_cell_Ah * c_rate * 1000/1090
    V_cell_arr_sum = [battery.V_cell_max]
    V_cell_arr = [battery.V_cell_max]
    t_arr = [0.0]
    SOC_arr = [100.0]
    charge_arr = [0.0]

    energy_discharge = 0
    while (V_cell_arr_sum[end] > battery.V_cell_min) && (t_arr[end] < (1/c_rate * 5)) 
    # while SOC_arr[end] > 0 && t_arr[end] < (1/c_rate * 5)
        push!(t_arr, t_arr[end] + delta_t)
        discharge_charge = I*t_arr[end]

        V = calc_V_cell_cc_aging(I*delta_t, I, n, T, battery)
        battery.current_charge -= I*delta_t
        V_sum = f_V_cell(discharge_charge, I, n, T, battery)
        energy_discharge += delta_t * I * (V_cell_arr_sum[end]+V_sum)/2

        push!(V_cell_arr_sum, V_sum)
        push!(V_cell_arr, V)
        push!(SOC_arr, battery.current_charge/battery.capacity_cell_Ah*100)
        push!(charge_arr, discharge_charge)
    end
    end_idx = length(t_arr)
    charge_end_discharge = I*t_arr[end]
    
    I = -I
    energy_charge = 0
    # while (V_cell_arr_sum[end] < 3.8) && (t_arr[end] < (1/c_rate * 5))
    while charge_arr[end] > 0 && t_arr[end] < (1/c_rate * 5)
        push!(t_arr, t_arr[end]+delta_t)
        charge_charge = I*(t_arr[end]-t_arr[end_idx]) + charge_end_discharge

        V = calc_V_cell_cc_aging(I*delta_t, I, n, T, battery)
        battery.current_charge -= I*delta_t

        V_sum = f_V_cell(charge_charge, I, n, T, battery)
        energy_charge -= delta_t * I * (V_cell_arr_sum[end]+V_sum)/2

        push!(V_cell_arr, V)
        push!(V_cell_arr_sum, V_sum)
        push!(SOC_arr, battery.current_charge/battery.capacity_cell_Ah*100) 
        push!(charge_arr, charge_charge)        
    end
    V = calc_V_cell_cc_aging(0, 0, n, T, battery)
    V_sum = f_V_cell(charge_arr[end], 0, n, T, battery)
    push!(V_cell_arr, V)
    push!(V_cell_arr_sum, V_sum)
    push!(SOC_arr, battery.current_charge/battery.capacity_cell_Ah*100) 
    push!(charge_arr, charge_arr[end])     


    # p.plot(size=(1200,600), legend=:right)
    # p.plot(legend=:bottom, xlabel="Charge / Ah", ylabel="Zellspannung / V", 
    #      ylims=(1.8,3.7), yticks=1.8:0.1:3.7, 
    #      xlims=(0,1150), xticks=0:100:1100, size=(600,400))
    # p.plot!(charge_arr[1:end_idx], V_cell_arr[1:end_idx], label="Entladen:$c_rate C; Energie: $(round(charge_arr[1]-charge_arr[end_idx], digits=4)) Ah")
    # p.plot!(charge_arr[end_idx:end], V_cell_arr[end_idx:end], label="Laden:$c_rate C; Energie: $(round(charge_arr[end]-charge_arr[end_idx], digits=4)) Ah")
    p.plot!(charge_arr[1:end_idx], V_cell_arr_sum[1:end_idx], label="Entladen_sum:$(c_rate)C,$(n)n,$(T)°C; Energie: $(round(energy_discharge, digits=4)) Wh")
    p.plot!(charge_arr[end_idx:end], V_cell_arr_sum[end_idx:end], label="Laden_sum:$(c_rate)C,$(n)n,$(T)°C; Energie: $(round(energy_charge, digits=4)) Wh")
    # v_hystersis_arr = reverse(V_cell_arr_sum[end_idx:end]) - V_cell_arr_sum[(end_idx+1-length(V_cell_arr_sum[end_idx:end])):end_idx]
    # p.plot(SOC_arr[end_idx:end], v_hystersis_arr)
    # p.plot!(t_arr, charge_arr)

end


# calc_and_plot_discharge_charge_curve_comp(5, 1/5/3600)
# calc_and_plot_discharge_charge_curve_comp(1, 1/3600)
# calc_and_plot_discharge_charge_curve_comp(0.1, 10/3600)

# V_0, K, A, B, r_i, m, alpha, 
# k_qn, k_qT, 
# k_n, 
# k_T, i_ref, T_ref
ap_paper = 3.36964, 0.03546, 0.08165, 0.1003, 0.00016, 1.0269, -0.01212, 
           [-1.27571e-7, 1.22095e-11], [1.32729e-3, -7.9763e-6], 
           [9.71249e-6, 7.51635e-4, -8.59363e-5, -2.92489e-4], [1.05135e-3, 1.83721e-2, -7.72438e-3, -4.31833e-2], 100, 25

p.plot(legend=:bottomleft, xlabel="Charge / Ah", ylabel="Zellspannung / V", 
     ylims=(1.8,3.7), ticks=:native,
     xlims=(0,1150), size=(1450,850))
calc_and_plot_discharge_charge_curve(0.5, 0.1/3600, ap_paper, 1)
calc_and_plot_discharge_charge_curve(0.5, 0.1/3600, ap_paper, 1, -20)
calc_and_plot_discharge_charge_curve(0.5, 0.1/3600, ap_paper, 1, 55)

# datasheet = CSV.read("profiles/SP_LFP1000AHA/0.1C.csv", DataFrame, header=["Q","V"])
# scatter!(datasheet.Q, datasheet.V, label="datasheet 0.1C", ms=2 , markerstrokeopacity=0)
# datasheet = CSV.read("profiles/SP_LFP1000AHA/0.3C.csv", DataFrame, header=["Q","V"])
# scatter!(datasheet.Q, datasheet.V, label="datasheet 0.3C", ms=2 , markerstrokeopacity=0)
# datasheet = CSV.read("profiles/SP_LFP1000AHA/1C.csv", DataFrame, header=["Q","V"])
# scatter!(datasheet.Q, datasheet.V, label="datasheet 1C", ms=2 , markerstrokeopacity=0)

# @testset "calc_and_plot_discharge_curve" begin
#     calc_and_plot_discharge_curve()
# end
