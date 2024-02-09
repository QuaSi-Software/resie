function merge(first::Bus, second::Bus, uac::String)::Bus
    new_bus = deepcopy(first)
    new_bus.uac = uac

    # insert outputs of second to outputs of first in place of bus->bus interface
    output_idx = [
        idx for (idx,f) in pairs(first.output_interfaces)
            if occursin(f.target.uac, second.uac)
    ][1]
    splice!(new_bus.output_interfaces, output_idx, second.output_interfaces)
    splice!(
        new_bus.connectivity.output_order,
        output_idx,
        second.connectivity.output_order
    )

    # insert inputs of first to inputs of second in place of bus->bus interface
    input_idx = [
        idx for (idx,f) in pairs(second.input_interfaces)
            if occursin(f.source.uac, first.uac)
    ][1]
    new_bus.input_interfaces = Base.deepcopy(second.input_interfaces)
    splice!(new_bus.input_interfaces, input_idx, first.input_interfaces)
    new_bus.connectivity.input_order = Base.deepcopy(second.connectivity.input_order)
    splice!(
        new_bus.connectivity.input_order,
        input_idx,
        first.connectivity.input_order
    )

    # append balance table inputs and rewrite priorities
    merge!(new_bus.balance_table_inputs, Base.deepcopy(second.balance_table_inputs))
    delete!(new_bus.balance_table_inputs, first.uac)
    for (idx,inface) in pairs(new_bus.input_interfaces)
        new_bus.balance_table_inputs[inface.source.uac].priority = idx
        new_bus.balance_table_inputs[inface.source.uac].input_index = idx
    end

    # append balance table outputs and rewrite priorities
    merge!(new_bus.balance_table_outputs, Base.deepcopy(second.balance_table_outputs))
    delete!(new_bus.balance_table_outputs, second.uac)
    for (idx,outface) in pairs(new_bus.output_interfaces)
        new_bus.balance_table_outputs[outface.target.uac].priority = idx
        new_bus.balance_table_outputs[outface.target.uac].output_index = idx
    end

    # merge energy_flow
    new_bus.connectivity.energy_flow = []

    for input_row in sort(collect(values(new_bus.balance_table_inputs)), by=x->x.priority)
        is_allowed_flags = []
        is_input_from_first = input_row.source.uac in keys(first.balance_table_inputs)

        for output_row in sort(collect(values(new_bus.balance_table_outputs)), by=x->x.priority)
            is_output_from_first = output_row.target.uac in keys(first.balance_table_outputs)

            if is_input_from_first && is_output_from_first
                input_index = first.balance_table_inputs[input_row.source.uac].input_index
                output_index = first.balance_table_outputs[output_row.target.uac].output_index
                is_allowed = first.connectivity.energy_flow[input_index][output_index]

            elseif is_input_from_first && !is_output_from_first
                input_index = first.balance_table_inputs[input_row.source.uac].input_index
                output_index = first.balance_table_outputs[second.uac].output_index
                is_allowed = first.connectivity.energy_flow[input_index][output_index]
                input_index = second.balance_table_inputs[first.uac].input_index
                output_index = second.balance_table_outputs[output_row.target.uac].output_index
                is_allowed = is_allowed && second.connectivity.energy_flow[input_index][output_index]

            elseif !is_input_from_first && is_output_from_first
                # inputs from the second bus can't see outputs of the first
                is_allowed = false

            else
                input_index = second.balance_table_inputs[input_row.source.uac].input_index
                output_index = second.balance_table_outputs[output_row.target.uac].output_index
                is_allowed = second.connectivity.energy_flow[input_index][output_index]
            end

            push!(is_allowed_flags, is_allowed)
        end

        push!(new_bus.connectivity.energy_flow, is_allowed_flags)
    end

    # reset recreates the balance table
    reset(new_bus)

    return new_bus
end
