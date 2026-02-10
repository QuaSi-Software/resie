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

function read_JSON(filepath::String)::OrderedDict{AbstractString,Any}
    open(filepath, "r") do file_handle
        content = read(file_handle, String)
        return JSON.parse(content; dicttype=OrderedDict)
    end
end

function run_simulation(inputs::OrderedDict{AbstractString,Any})::NamedTuple
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
            hp_power = sample[4] + rand((-0.5 * range):(0.5 * range))
            hp_power = max(BOUNDS["TST_HP_01"][2], min(BOUNDS["TST_HP_01"][3], hp_power))

            range = NBH_SCALE * temperature * (BOUNDS["TST_PV_01"][3] - BOUNDS["TST_PV_01"][2])
            pv_power = sample[5] + rand((-0.5 * range):(0.5 * range))
            pv_power = max(BOUNDS["TST_PV_01"][2], min(BOUNDS["TST_PV_01"][3], pv_power))

            range = NBH_SCALE * temperature * (BOUNDS["TST_BAT_01"][3] - BOUNDS["TST_BAT_01"][2])
            bat_cap = sample[6] + rand((-0.5 * range):(0.5 * range))
            bat_cap = max(BOUNDS["TST_BAT_01"][2], min(BOUNDS["TST_BAT_01"][3], bat_cap))

            range = NBH_SCALE * temperature * (BOUNDS["TST_BFT_01"][3] - BOUNDS["TST_BFT_01"][2])
            bt_cap = sample[7] + rand((-0.5 * range):(0.5 * range))
            bt_cap = max(BOUNDS["TST_BFT_01"][2], min(BOUNDS["TST_BFT_01"][3], bt_cap))

            inputs["components"]["TST_HP_01"]["power_th"] = hp_power
            inputs["components"]["TST_PV_01"]["scale"] = pv_power
            inputs["components"]["TST_BAT_01"]["capacity"] = bat_cap
            inputs["components"]["TST_BFT_01"]["capacity"] = bt_cap
        end

        # run sim and record results
        run_results = run_simulation(inputs)
        obj_results = objective_function(run_results)
        push!(all_results, [0.0, obj_results.lcc, obj_results.emissions, hp_power, pv_power, bat_cap, bt_cap])
        min_lcc = obj_results.lcc < min_lcc ? obj_results.lcc : min_lcc
        min_emissions = obj_results.emissions < min_emissions ? obj_results.emissions : min_emissions

        # calculate global measure and sort by it
        for res in all_results
            res[1] = sqrt(((res[2]) / min_lcc - 1.0)^2 + (res[3] / min_emissions - 1.0)^2)
        end
        sort!(all_results; by=x -> x[1])

        if idx % 5 == 0
            print("| ")
        end
    end
    print("\nFinished simulation, now plotting\n")

    # create scatter plot with objective function values
    x_values = [result[2] for result in all_results]
    y_values = [result[3] for result in all_results]

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
            if i != j && result_j[2] <= result_i[2] && result_j[3] <= result_i[3]
                if result_j[2] < result_i[2] || result_j[3] < result_i[3]
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
        _, cost, emissions, hp_power, pv_power, bat_cap, bt_cap = point
        println("Cost: $(cost), Emissions: $(emissions)")
        println("  HP power: $(hp_power), PV scale: $(pv_power), BAT capacity: $(bat_cap), BFT capacity: $(bt_cap)")
    end

    println("Top 10 Global Measure Points:")
    for point in all_results[1:10]
        gb, cost, emissions, hp_power, pv_power, bat_cap, bt_cap = point
        println("Global Measure: $gb, Cost: $(cost), Emissions: $(emissions)")
        println("  HP power: $(hp_power), PV scale: $(pv_power), BAT capacity: $(bat_cap), BFT capacity: $(bt_cap)")
    end
end

main()
