# Security-related limits and helpers shared by project loading and output generation.
# These limits are intentionally conservative and prevent accidental or malicious
# configurations from allocating unbounded resources inside a ReSiE process.
const MAX_CONFIG_FILE_BYTES = 16 * 1024 * 1024
const MAX_JSON_NESTING_DEPTH = 200
const MAX_INPUT_FILE_BYTES = 128 * 1024 * 1024
const MAX_PROFILE_FILE_BYTES = 64 * 1024 * 1024
const MAX_PROFILE_LINES = 2_000_000
const MAX_LINE_BYTES = 1 * 1024 * 1024
const MAX_SIMULATION_STEPS = 100_000_000
const MAX_COMPONENTS = 10_000
const MAX_COMPONENT_IDENTIFIER_LENGTH = 128
const MAX_CONTROL_MODULES_PER_COMPONENT = 100
const MAX_ORDER_OPERATIONS = 100_000
const MAX_PARAMETER_STUDY_RUNS = 10_000
const MAX_PARAMETER_STUDY_PARAMETERS = 100
const MAX_OBJECTIVE_VALUES = 1_000
const MAX_PARAMETER_VALUES = 100_000
const MAX_GRAPH_PATHS = 100_000
const MAX_GRAPH_DEPTH = 500
const MAX_MESH_CELLS = 1_000_000
const MAX_STORED_VALUES = 1_000_000_000
const MAX_G_FUNCTION_ROWS = 1_000_000

const PATH_MODES = (:local, :confined)

function validate_path_mode(path_mode)::Symbol
    path_mode isa Symbol || throw(ArgumentError("Path mode must be `:local` or `:confined`."))
    path_mode in PATH_MODES || throw(ArgumentError("Path mode must be `:local` or `:confined`."))
    return path_mode
end

function validate_path_string(path)::String
    path isa AbstractString || throw(InputError("File paths must be strings."))

    raw_path = String(path)
    occursin(r"[\x00-\x1f\x7f]", raw_path) &&
        throw(InputError("File paths must not contain control characters."))
    ncodeunits(raw_path) <= 4096 || throw(InputError("File path exceeds the supported length."))
    isempty(raw_path) && throw(InputError("File path must not be empty."))

    return raw_path
end

function path_is_within(path::AbstractString, root::AbstractString)::Bool
    root_abs = abspath(root)
    path_abs = abspath(path)
    relative = relpath(path_abs, root_abs)
    parts = splitpath(relative)
    return relative == "." || (!isabspath(relative) && !isempty(parts) && first(parts) != "..")
end

function normalise_relative_path(path)::String
    raw_path = validate_path_string(path)
    normalised = replace(strip(raw_path), '\\' => '/')
    isempty(normalised) && throw(InputError("File path must not be empty."))

    if isabspath(normalised) || startswith(normalised, "//") || occursin(r"^[A-Za-z]:", normalised)
        throw(InputError("Absolute file paths are not allowed in confined mode."))
    end
    any(part -> part == "..", split(normalised, '/')) &&
        throw(InputError("File paths must not contain parent-directory components in confined mode."))

    return normpath(normalised)
end

function reject_symlink_components(root::AbstractString, path::AbstractString)
    root_abs = abspath(root)
    path_abs = abspath(path)
    path_is_within(path_abs, root_abs) ||
        throw(InputError("File path must remain inside the configured directory."))

    current = root_abs
    relative = relpath(path_abs, root_abs)
    relative == "." && return

    for part in splitpath(relative)
        current = joinpath(current, part)
        if ispath(current) && islink(current)
            throw(InputError("Symbolic links are not allowed in confined file paths."))
        end
    end
end

function resolve_project_file(path::AbstractString;
                              max_bytes::Integer=MAX_CONFIG_FILE_BYTES,
                              confined::Bool=false)::String
    requested = abspath(validate_path_string(path))
    isfile(requested) || throw(InputError("Project config must point to a regular file: $requested"))
    confined && islink(requested) &&
        throw(InputError("Project config must not be a symbolic link in confined mode."))

    resolved = realpath(requested)
    filesize(resolved) <= max_bytes ||
        throw(InputError("Project config exceeds the maximum allowed file size."))
    return resolved
end

function resolve_local_directory(path)::String
    requested = abspath(validate_path_string(path))
    isdir(requested) || throw(InputError("Directory path does not point to a directory: $requested"))
    return realpath(requested)
end

function resolve_project_directory(config_root::AbstractString, path)::String
    root = realpath(config_root)
    relative = normalise_relative_path(path)
    candidate = abspath(joinpath(root, relative))
    path_is_within(candidate, root) ||
        throw(InputError("Base path must remain inside the confined input directory."))
    reject_symlink_components(root, candidate)
    isdir(candidate) || throw(InputError("Configured base path does not point to a directory: $candidate"))
    return realpath(candidate)
end

function resolve_local_input_path(base_path::AbstractString, path;
                                  max_bytes::Integer=MAX_INPUT_FILE_BYTES)::String
    raw_path = validate_path_string(path)
    candidate = isabspath(raw_path) ? normpath(raw_path) : abspath(joinpath(base_path, raw_path))
    validate_regular_input(candidate; max_bytes=max_bytes, allow_symlink=true)
    return realpath(candidate)
end

function resolve_input_path(root::AbstractString, path;
                            max_bytes::Integer=MAX_INPUT_FILE_BYTES)::String
    root_real = realpath(root)
    relative = normalise_relative_path(path)
    candidate = abspath(joinpath(root_real, relative))
    path_is_within(candidate, root_real) ||
        throw(InputError("Input path must remain inside the confined input directory."))
    reject_symlink_components(root_real, candidate)

    isfile(candidate) || throw(InputError("Input path must point to a regular file: $candidate"))
    resolved = realpath(candidate)
    path_is_within(resolved, root_real) ||
        throw(InputError("Input path resolves outside the confined input directory."))
    filesize(resolved) <= max_bytes ||
        throw(InputError("Input file exceeds the maximum allowed file size: $resolved"))
    return resolved
end

function remove_output_prefix(path::String)::String
    parts = splitpath(path)
    if !isempty(parts) && lowercase(first(parts)) == "output"
        return length(parts) == 1 ? "." : joinpath(parts[2:end]...)
    end
    return path
end

function resolve_local_output_path(base_path::AbstractString, path)::String
    raw_path = validate_path_string(path)
    candidate = isabspath(raw_path) ? normpath(raw_path) : abspath(joinpath(base_path, raw_path))
    mkpath(dirname(candidate))
    return candidate
end

function resolve_output_path(root::AbstractString, path)::String
    root_abs = abspath(root)
    mkpath(root_abs)
    relative = remove_output_prefix(normalise_relative_path(path))
    candidate = abspath(joinpath(root_abs, relative))
    path_is_within(candidate, root_abs) ||
        throw(InputError("Output path must remain inside the confined output directory."))
    reject_symlink_components(root_abs, candidate)
    mkpath(dirname(candidate))
    reject_symlink_components(root_abs, candidate)
    return candidate
end

function validate_regular_input(path::AbstractString;
                                max_bytes::Integer=MAX_INPUT_FILE_BYTES,
                                allow_symlink::Bool=false)::String
    isfile(path) || throw(InputError("Input path must point to a regular file: $path"))
    !allow_symlink && islink(path) &&
        throw(InputError("Symbolic links are not allowed for input files in confined mode."))
    filesize(path) <= max_bytes || throw(InputError("Input file exceeds the maximum allowed file size: $path"))
    return String(path)
end

function validate_json_nesting(content::AbstractString)::Nothing
    depth = 0
    in_string = false
    escaped = false

    for character in content
        if in_string
            if escaped
                escaped = false
            elseif character == '\\'
                escaped = true
            elseif character == '"'
                in_string = false
            end
        elseif character == '"'
            in_string = true
        elseif character == '{' || character == '['
            depth += 1
            depth <= MAX_JSON_NESTING_DEPTH ||
                throw(InputError("Project config exceeds the maximum JSON nesting depth."))
        elseif character == '}' || character == ']'
            depth -= 1
            depth >= 0 || throw(InputError("Project config contains unbalanced JSON delimiters."))
        end
    end

    depth == 0 || throw(InputError("Project config contains unbalanced JSON delimiters."))
    return nothing
end

function bounded_read_string(path::AbstractString, max_bytes::Integer)::String
    validate_regular_input(path; max_bytes=max_bytes)
    open(path, "r") do io
        data = read(io, max_bytes + 1)
        length(data) <= max_bytes || throw(InputError("Input file exceeds the maximum allowed file size."))
        return String(data)
    end
end

function bounded_eachline(path::AbstractString;
                          max_bytes::Integer=MAX_INPUT_FILE_BYTES,
                          max_lines::Integer=MAX_PROFILE_LINES,
                          max_line_bytes::Integer=MAX_LINE_BYTES,
                          callback::Function)
    validate_regular_input(path; max_bytes=max_bytes)
    open(path, "r") do io
        line_count = 0
        for line in eachline(io)
            line_count += 1
            line_count <= max_lines || throw(InputError("Input file contains too many lines."))
            ncodeunits(line) <= max_line_bytes || throw(InputError("Input file contains an excessively long line."))
            callback(line, line_count)
        end
    end
end

function checked_product(values::AbstractVector{<:Integer}; limit::Integer, label::String)::Int
    product = 1
    for value in values
        value >= 0 || throw(InputError("$label contains a negative dimension."))
        if value != 0 && product > div(limit, value)
            throw(InputError("$label exceeds the supported size."))
        end
        product *= value
    end
    return product
end

# Unknown JSON fields are permitted for annotations, metadata, and forward compatibility.
function report_unknown_keys(config::AbstractDict, allowed_keys, context::AbstractString;
                             extra_keys=String[])
    allowed = Set(String.(collect(allowed_keys)))
    union!(allowed, String.(extra_keys))
    unsupported = sort([String(key) for key in keys(config) if String(key) ∉ allowed]; by=lowercase)
    isempty(unsupported) ||
        @debug "Ignoring additional JSON fields" context=context fields=unsupported
    return unsupported
end

function validate_finite(value, name::AbstractString)
    if value isa AbstractFloat && !isfinite(value)
        throw(InputError("Parameter `$name` must be finite."))
    elseif value isa AbstractArray
        for item in value
            validate_finite(item, name)
        end
    elseif value isa AbstractDict
        for item in values(value)
            validate_finite(item, name)
        end
    end
    return value
end

function safe_filename(value; fallback::String="item", max_length::Integer=80)::String
    text = replace(strip(string(value)), r"[^A-Za-z0-9_.-]+" => "_")
    text = replace(text, r"^[._-]+|[._-]+$" => "")
    isempty(text) && return fallback
    characters = collect(text)
    return join(characters[1:min(length(characters), max_length)])
end

function sanitise_log_text(value)::String
    return replace(string(value), '\0' => "", '\r' => "\\r", '\n' => "\\n")
end

function csv_cell(value; decimal_comma::Bool=false)::String
    text = string(value)
    stripped = lstrip(text)
    numeric_value = value isa Number || tryparse(Float64, stripped) !== nothing

    if decimal_comma && numeric_value
        text = replace(text, '.' => ',')
        stripped = lstrip(text)
    end

    if !numeric_value && !isempty(stripped) && first(stripped) in ('=', '+', '-', '@', '\t', '\r')
        text = "'" * text
    end

    text = replace(text, '"' => "\"\"")
    return "\"" * text * "\""
end
