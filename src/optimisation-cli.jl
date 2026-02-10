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

    # create scatter plot with objective function values
    x_values = [result["cost"] for result in all_results]
    y_values = [result["emissions"] for result in all_results]

    scatter(x_values, y_values; xlabel="Life-cycle cost (capex + $OPEX_YEARS*opex) [€]",
            ylabel="Operational emissions [kg/a]", title="Pareto Front", legend=false,
            xlims=(4e4, 9e4), ylims=(0, 3500))
    savefig("./output/pareto_front.html")

    println("Plotting finished, now finding Pareto-optimal points")
    # find and print Pareto optimal points
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
