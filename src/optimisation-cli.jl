using JSON: JSON
using OrderedCollections: OrderedDict
using Random
using Plots

const BOUNDS = Dict{String,Tuple}(
    "TST_HP_01" => ("power_th", 3000, 12000),
    "TST_PV_01" => ("scale", 2000, 4500),
    "TST_BAT_01" => ("capacity", 0, 30000),
    "TST_BFT_01" => ("capacity", 0, 60000),
)
const NR_TRIES = 1000
const OPEX_YEARS = 20
const HEAT_DEMAND = 12000000
const POWER_DEMAND = 3000000

function read_JSON(filepath::String)::OrderedDict{AbstractString,Any}
    open(filepath, "r") do file_handle
        content = read(file_handle, String)
        return JSON.parse(content; dicttype=OrderedDict)
    end
end

function run_simulation(inputs::OrderedDict{AbstractString,Any})::Dict{String,Float64}
    # artifical delay, useful for debugging
    # sleep(0.05)

    # naive capex, opex and emissions calculations depending purely on parameters as we
    # don't perform a simulation
    capex = inputs["components"]["TST_HP_01"]["power_th"] * 0.001 * 1500.0 +
            inputs["components"]["TST_BAT_01"]["capacity"] * 0.001 * 350.0 +
            inputs["components"]["TST_BFT_01"]["capacity"] * 0.001 * 150.0 +
            inputs["components"]["TST_PV_01"]["scale"] * 0.001 * 1000.0 +
            10000.0

    avg_cop = 2.0 + (inputs["components"]["TST_BFT_01"]["capacity"] / BOUNDS["TST_BFT_01"][3])
    hp_load_hours = HEAT_DEMAND / inputs["components"]["TST_HP_01"]["power_th"]
    power_demand = POWER_DEMAND + (inputs["components"]["TST_HP_01"]["power_th"] / avg_cop) * hp_load_hours
    grid_power_demand = power_demand -
                        (inputs["components"]["TST_PV_01"]["scale"] * 2000) *
                        (inputs["components"]["TST_BAT_01"]["capacity"] / BOUNDS["TST_BAT_01"][3])
    grid_power_demand = max(grid_power_demand, 0.0)
    price = 0.3 - 0.1 * (inputs["components"]["TST_BAT_01"]["capacity"] / BOUNDS["TST_BAT_01"][3])
    opex = grid_power_demand * 0.001 * price
    opex = max(opex, 0.0)

    results = Dict{String,Float64}(
        "capex" => capex,
        "opex" => opex,
        "emissions" => (0.35 - 0.1 * (inputs["components"]["TST_BAT_01"]["capacity"] /
                                      BOUNDS["TST_BAT_01"][3])) * grid_power_demand * 0.001,
    )
    return results
end

function objective_function(results::Dict{String,Float64})::Tuple{Float64,Float64}
    return results["capex"] + OPEX_YEARS * results["opex"], results["emissions"]
end

function main()
    inputs = read_JSON("./examples/optimisation_example.json")

    # simple monte carlo
    all_results = []
    for _ in range(0, NR_TRIES)
        # set parameters to equally distributed random values
        hp_power = rand(BOUNDS["TST_HP_01"][2]:BOUNDS["TST_HP_01"][3])
        pv_power = rand(BOUNDS["TST_PV_01"][2]:BOUNDS["TST_PV_01"][3])
        bat_cap = rand(BOUNDS["TST_BAT_01"][2]:BOUNDS["TST_BAT_01"][3])
        bt_cap = rand(BOUNDS["TST_BFT_01"][2]:BOUNDS["TST_BFT_01"][3])
        inputs["components"]["TST_HP_01"]["power_th"] = hp_power
        inputs["components"]["TST_PV_01"]["scale"] = pv_power
        inputs["components"]["TST_BAT_01"]["capacity"] = bat_cap
        inputs["components"]["TST_BFT_01"]["capacity"] = bt_cap

        # run sim and record results
        run_results = run_simulation(inputs)
        obj_results = objective_function(run_results)
        push!(all_results, (obj_results[1], obj_results[2], hp_power, pv_power, bat_cap, bt_cap))

        # debug output
        # print("run with HP $(hp_power) PV $(pv_power) BAT $(bat_cap) BFT $(bt_cap) ")
        # print("results in $(run_results["capex"]) $(run_results["opex"]) $(run_results["emissions"])\n")
    end

    # create scatter plot with objective function values
    x_values = [result[1] for result in all_results]
    y_values = [result[2] for result in all_results]

    scatter(x_values, y_values; xlabel="Cost (capex + 20*opex)", ylabel="Emissions",
            title="Pareto Front", legend=false)
    savefig("./output/pareto_front.html")

    # find and print Pareto optimal points
    pareto_points = []
    for (i, result_i) in enumerate(all_results)
        is_dominated = false
        for (j, result_j) in enumerate(all_results)
            if i != j && result_j[1] <= result_i[1] && result_j[2] <= result_i[2]
                if result_j[1] < result_i[1] || result_j[2] < result_i[2]
                    is_dominated = true
                    break
                end
            end
        end
        if !is_dominated
            push!(pareto_points, result_i)
        end
    end

    println("\nPareto Front Points:")
    for point in pareto_points
        cost, emissions, hp_power, pv_power, bat_cap, bt_cap = point
        println("Cost: $(cost), Emissions: $(emissions)")
        println("  HP power: $(hp_power), PV scale: $(pv_power), BAT capacity: $(bat_cap), BFT capacity: $(bt_cap)")
    end
end

main()
