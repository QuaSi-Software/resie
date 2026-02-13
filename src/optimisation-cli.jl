using JSON: JSON
using OrderedCollections: OrderedDict
using Random
using Plots
using Measures

const BOUNDS = Dict{String,Tuple}(
    "TST_HP_01" => ("power_th", 3000, 12000),
    "TST_PV_01" => ("scale", 2000, 4500),
    "TST_BAT_01" => ("capacity", 0, 30000),
    "TST_BFT_01" => ("capacity", 0, 60000),
)
const NBH_SCALE = 0.5
const NR_TRIES = 200
const OPEX_YEARS = 20
const HEAT_DEMAND = 12000000
const POWER_DEMAND = 5000000
const SUNSHINE_HOURS = 2000

function read_JSON(filepath::String)::OrderedDict{AbstractString,Any}
    open(filepath, "r") do file_handle
        content = read(file_handle, String)
        return JSON.parse(content; dicttype=OrderedDict)
    end
end

function run_simulation(inputs::OrderedDict{AbstractString,Any})::NamedTuple
    # naive capex, opex and emissions calculations depending purely on parameters as we
    # don't perform a simulation
    capex = inputs["components"]["TST_HP_01"]["power_th"] * 0.001 *
            inputs["components"]["TST_HP_01"]["capex_per_kW"] +
            inputs["components"]["TST_BAT_01"]["capacity"] * 0.001 *
            inputs["components"]["TST_BAT_01"]["capex_per_kWh"] +
            inputs["components"]["TST_BFT_01"]["capacity"] * 0.001 *
            inputs["components"]["TST_BFT_01"]["capex_per_kWh"] +
            inputs["components"]["TST_PV_01"]["scale"] * 0.001 *
            inputs["components"]["TST_PV_01"]["capex_per_kW"] +
            10000.0

    avg_cop = inputs["components"]["TST_HP_01"]["base_cop"] +
              (inputs["components"]["TST_BFT_01"]["capacity"] / BOUNDS["TST_BFT_01"][3])
    hp_load_hours = HEAT_DEMAND / inputs["components"]["TST_HP_01"]["power_th"]
    power_demand = POWER_DEMAND + (inputs["components"]["TST_HP_01"]["power_th"] / avg_cop) * hp_load_hours
    grid_power_demand = power_demand -
                        (inputs["components"]["TST_PV_01"]["scale"] * SUNSHINE_HOURS) *
                        (inputs["components"]["TST_BAT_01"]["capacity"] / BOUNDS["TST_BAT_01"][3])
    grid_power_demand = max(grid_power_demand, 0.0)
    price = inputs["components"]["TST_BAT_01"]["base_grid_price"] -
            0.1 * (inputs["components"]["TST_BAT_01"]["capacity"] / BOUNDS["TST_BAT_01"][3])
    opex = grid_power_demand * 0.001 * price
    opex = max(opex, 0.0)

    return (capex=capex,
            opex=opex,
            emissions=(0.35 - 0.1 * (inputs["components"]["TST_BAT_01"]["capacity"] /
                                     BOUNDS["TST_BAT_01"][3])) * grid_power_demand * 0.001)
end

function objective_function(results::NamedTuple)::NamedTuple
    return (lcc=(results.capex + OPEX_YEARS * results.opex), emissions=results.emissions)
end

optim



function main()
    inputs = read_JSON("./examples/optimisation_example.json")

    print("Starting N=$NR_TRIES simulation runs: ")

    all_results = []
    min_lcc = Inf
    min_emissions = Inf

    # combined monte carlo and simulated annealing
    # the temperature determines if a completely random or existing sample is used as starting
    # point, determines the size of the neighborhood and the number of results (sorted by)
    # global measure, from which a new sample is drawn
    for idx in range(1, NR_TRIES)
        # temperature schedule is simple inverse logistic curve
        temperature = 1.0 - 1.0 / (1.0 + exp(-8.0 * (idx / NR_TRIES - 0.5)))

        if length(all_results) == 0 || rand() < temperature
            # set parameters to equally distributed random values across whole parameter space
            hp_power = rand(BOUNDS["TST_HP_01"][2]:BOUNDS["TST_HP_01"][3])
            pv_power = rand(BOUNDS["TST_PV_01"][2]:BOUNDS["TST_PV_01"][3])
            bat_cap = rand(BOUNDS["TST_BAT_01"][2]:BOUNDS["TST_BAT_01"][3])
            bt_cap = rand(BOUNDS["TST_BFT_01"][2]:BOUNDS["TST_BFT_01"][3])
            inputs["components"]["TST_HP_01"]["power_th"] = hp_power
            inputs["components"]["TST_PV_01"]["scale"] = pv_power
            inputs["components"]["TST_BAT_01"]["capacity"] = bat_cap
            inputs["components"]["TST_BFT_01"]["capacity"] = bt_cap
            print(". ")
        else
            # set parameters to neighborhood of existing result, drawn from the top results
            # by global measure, where temperature determines the results pool and size of
            # neighborhood
            sample_idx = rand(1:max(1, Int(round(length(all_results) * temperature))))
            print("$sample_idx ")
            sample = sample_idx >= 1 && sample_idx <= length(all_results) ? all_results[sample_idx] : all_results[1]

            range = NBH_SCALE * temperature * (BOUNDS["TST_HP_01"][3] - BOUNDS["TST_HP_01"][2])
            hp_power = sample["hp_power"] + rand((-0.5 * range):(0.5 * range))
            hp_power = max(BOUNDS["TST_HP_01"][2], min(BOUNDS["TST_HP_01"][3], hp_power))

            range = NBH_SCALE * temperature * (BOUNDS["TST_PV_01"][3] - BOUNDS["TST_PV_01"][2])
            pv_power = sample["pv_power"] + rand((-0.5 * range):(0.5 * range))
            pv_power = max(BOUNDS["TST_PV_01"][2], min(BOUNDS["TST_PV_01"][3], pv_power))

            range = NBH_SCALE * temperature * (BOUNDS["TST_BAT_01"][3] - BOUNDS["TST_BAT_01"][2])
            bat_cap = sample["bat_cap"] + rand((-0.5 * range):(0.5 * range))
            bat_cap = max(BOUNDS["TST_BAT_01"][2], min(BOUNDS["TST_BAT_01"][3], bat_cap))

            range = NBH_SCALE * temperature * (BOUNDS["TST_BFT_01"][3] - BOUNDS["TST_BFT_01"][2])
            bt_cap = sample["bt_cap"] + rand((-0.5 * range):(0.5 * range))
            bt_cap = max(BOUNDS["TST_BFT_01"][2], min(BOUNDS["TST_BFT_01"][3], bt_cap))

            inputs["components"]["TST_HP_01"]["power_th"] = hp_power
            inputs["components"]["TST_PV_01"]["scale"] = pv_power
            inputs["components"]["TST_BAT_01"]["capacity"] = bat_cap
            inputs["components"]["TST_BFT_01"]["capacity"] = bt_cap
        end

        # run sim and record results
        run_results = run_simulation(inputs)
        obj_results = objective_function(run_results)
        result = Dict(
            "gm" => 0.0,
            "cost" => obj_results.lcc,
            "emissions" => obj_results.emissions,
            "hp_power" => hp_power,
            "hp_capex_per_kw" => inputs["components"]["TST_HP_01"]["capex_per_kW"],
            "pv_power" => pv_power,
            "pv_capex_per_kw" => inputs["components"]["TST_PV_01"]["capex_per_kW"],
            "bat_cap" => bat_cap,
            "bat_capex_per_kwh" => inputs["components"]["TST_BAT_01"]["capex_per_kWh"],
            "bt_cap" => bt_cap,
            "bt_capex_per_kwh" => inputs["components"]["TST_BFT_01"]["capex_per_kWh"],
            "base_cop" => inputs["components"]["TST_HP_01"]["base_cop"],
            "base_grid_price" => inputs["components"]["TST_BAT_01"]["base_grid_price"],
        )
        push!(all_results, result)
        min_lcc = obj_results.lcc < min_lcc ? obj_results.lcc : min_lcc
        min_emissions = obj_results.emissions < min_emissions ? obj_results.emissions : min_emissions

        # calculate global measure and sort by it
        for res in all_results
            res["gm"] = sqrt(((res["cost"]) / min_lcc - 1.0)^2 + (res["emissions"] / min_emissions - 1.0)^2)
        end
        sort!(all_results; by=x -> x["gm"])

        if idx % 5 == 0
            print("| ")
        end
    end
    print("\nFinished simulation, now plotting\n")

    cost_values = [result["cost"] for result in all_results]
    em_values = [result["emissions"] for result in all_results]
    hp_powers = [result["hp_power"] for result in all_results]
    pv_powers = [result["pv_power"] for result in all_results]
    bat_caps = [result["bat_cap"] for result in all_results]
    bt_caps = [result["bt_cap"] for result in all_results]

    p1 = scatter(cost_values, em_values; xlabel="Life-cycle cost (capex + $OPEX_YEARS*opex) [€]",
                 ylabel="Operational emissions [kg/a]", title="Pareto Front", legend=false,
                 xlims=(4e4, 9e4), ylims=(0, 3500), left_margin=18mm, right_margin=1.8mm,
                 bottom_margin=8mm, top_margin=8mm)
    p2 = scatter([0.0], [0.0]; xlabel="-", ylabel="-", title="Placeholder", legend=false,
                 left_margin=13mm, right_margin=1.8mm, bottom_margin=8mm, top_margin=8mm)

    p3 = scatter(hp_powers, cost_values; xlabel="HP power [W]",
                 ylabel="Life-cycle cost [€]", title="HP Power vs Cost", legend=false,
                 left_margin=18mm, right_margin=1.8mm, bottom_margin=8mm, top_margin=8mm)
    p4 = scatter(hp_powers, em_values; xlabel="HP power [W]",
                 ylabel="Operational emissions [kg/a]", title="HP Power vs Emissions", legend=false,
                 left_margin=13mm, right_margin=1.8mm, bottom_margin=8mm, top_margin=8mm)

    p5 = scatter(pv_powers, cost_values; xlabel="PV power [W]",
                 ylabel="Life-cycle cost [€]", title="PV Power vs Cost", legend=false,
                 left_margin=18mm, right_margin=1.8mm, bottom_margin=8mm, top_margin=8mm)
    p6 = scatter(pv_powers, em_values; xlabel="PV power [W]",
                 ylabel="Operational emissions [kg/a]", title="PV Power vs Emissions", legend=false,
                 left_margin=13mm, right_margin=1.8mm, bottom_margin=8mm, top_margin=8mm)

    p7 = scatter(bat_caps, cost_values; xlabel="Battery capacity [Wh]",
                 ylabel="Life-cycle cost [€]", title="Battery Capacity vs Cost", legend=false,
                 left_margin=18mm, right_margin=1.8mm, bottom_margin=8mm, top_margin=8mm)
    p8 = scatter(bat_caps, em_values; xlabel="Battery capacity [Wh]",
                 ylabel="Operational emissions [kg/a]", title="Battery Capacity vs Emissions", legend=false,
                 left_margin=13mm, right_margin=1.8mm, bottom_margin=8mm, top_margin=8mm)

    p9 = scatter(bt_caps, cost_values; xlabel="BFT capacity [Wh]",
                 ylabel="Life-cycle cost [€]", title="BFT Capacity vs Cost", legend=false,
                 left_margin=18mm, right_margin=1.8mm, bottom_margin=8mm, top_margin=8mm)
    p10 = scatter(bt_caps, em_values; xlabel="BFT capacity [Wh]",
                  ylabel="Operational emissions [kg/a]", title="BFT Capacity vs Emissions", legend=false,
                  left_margin=13mm, right_margin=1.8mm, bottom_margin=8mm, top_margin=8mm)

    combined_plot = plot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10; layout=(5, 2), size=(1200, 2500))
    savefig(combined_plot, "./output/pareto_front.html")

    println("Plotting finished, now finding Pareto-optimal points")
    pareto_points = []
    for (i, result_i) in enumerate(all_results)
        is_dominated = false
        for (j, result_j) in enumerate(all_results)
            if i != j && result_j["cost"] <= result_i["cost"] && result_j["emissions"] <= result_i["emissions"]
                if result_j["cost"] < result_i["cost"] || result_j["emissions"] < result_i["emissions"]
                    is_dominated = true
                    break
                end
            end
        end
        if !is_dominated
            push!(pareto_points, result_i)
        end
    end

    println("Pareto Front Points:")
    for point in pareto_points
        println("Cost: $(point["cost"]), Emissions: $(point["emissions"])")
        println("  HP power: $(point["hp_power"]), PV scale: $(point["pv_power"]), " *
                "BAT capacity: $(point["bat_cap"]), BFT capacity: $(point["bt_cap"])")
    end

    println("Top 10 Global Measure Points:")
    for point in all_results[1:10]
        println("Global Measure: $(point["gm"]), Cost: $(point["cost"]), Emissions: $(point["emissions"])")
        println("  HP power: $(point["hp_power"]), PV scale: $(point["pv_power"]), " *
                "BAT capacity: $(point["bat_cap"]), BFT capacity: $(point["bt_cap"])")
    end
end

main()
