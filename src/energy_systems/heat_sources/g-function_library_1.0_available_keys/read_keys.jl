## standalone tool to read out the keys of the g-fuction library
using JSON

function write_json_keys_preserving_order(json_filepath::String, output_filepath::String, max_depth::Int)
    # Parse the JSON file
    json_content = JSON.parsefile(json_filepath)
    keys = []
    # Function to recursively write keys
    function _collect_keys(obj, depth=1, parent_key="")
        if depth > max_depth
            return
        end
        if isa(obj, Dict)
            for (key, value) in obj
                # Construct a composite key if nested
                composite_key = parent_key == "" ? key : parent_key * "." * key
                if depth == max_depth
                    push!(keys, composite_key)
                end
                _collect_keys(value, depth + 1, composite_key)
            end
        elseif isa(obj, Array)
            for (index, value) in enumerate(obj)
                # Handle array elements, appending the index to the parent key
                array_key = parent_key * "[$index]"
                _collect_keys(value, depth + 1, array_key)
            end
        end
    end

    # collekt and write keys to the output file
    _collect_keys(json_content)
    sort!(keys, lt=custom_sort_hierarchical)
    open(output_filepath, "w") do io
        for key in keys
            println(io, key)
        end
    end
end

# Custom sort function
function custom_sort_hierarchical(s1::String, s2::String)
    # Helper function to convert a segment into a tuple of integers
    function parse_segment(segment)
        return Tuple(parse(Int, sub) for sub in split(segment, "."))
    end
    
    # Split the strings by underscore, and convert each segment to a tuple of integers
    nums1 = map(parse_segment, split(s1, "_"))
    nums2 = map(parse_segment, split(s2, "_"))
    
    # Compare the sequences of tuples
    min_length = min(length(nums1), length(nums2))
    for i in 1:min_length
        if nums1[i] != nums2[i]
            return nums1[i] < nums2[i]  # Tuple comparison
        end
    end
    
    # If all compared tuples are equal, the shorter sequence comes first
    return length(nums1) < length(nums2)
end



probe_field_configurations = Dict("rectangle" => "rectangle_5m_v1.0.json",
                                  "open_rectangle"=> "Open_configurations_5m_v1.0.json",
                                  "zoned_rectangle" => "zoned_rectangle_5m_v1.0.json",
                                  "U_configurations" => "U_configurations_5m_v1.0.json",
                                  "lopsided_U_configuration" => "LopU_configurations_5m_v1.0.json",
                                  "C_configuration" => "C_configurations_5m_v1.0.json",
                                  "L_configuration" => "L_configurations_5m_v1.0.json"
                                  )

 depth = Dict("rectangle" => 1,
             "open_rectangle"=> 2,
             "zoned_rectangle" => 2,
             "U_configurations" => 2,
             "lopsided_U_configuration" => 2,
             "C_configuration" => 2,
             "L_configuration" => 1
             )

for library in keys(probe_field_configurations)
    probe_field_geometry = library
    libfile_path = "src/energy_systems/heat_sources/g-function_library_1.0/" * probe_field_configurations[probe_field_geometry]
    outputfile_path = "src/energy_systems/heat_sources/g-function_library_1.0_available_keys/" * probe_field_geometry * "_keys.txt"
    write_json_keys_preserving_order(libfile_path, outputfile_path, depth[probe_field_geometry])
end