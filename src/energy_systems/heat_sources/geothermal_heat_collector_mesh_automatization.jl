# Mesh automatization

# TODO: Modify "Special" Node(s) between increasing and decreasing mesh width. Should not be smaller than previous mesh width or min. mesh width.
function create_mesh_y(min_mesh_width,
                       max_mesh_width,
                       expansion_factor,
                       pipe_laying_depth,
                       pipe_diameter_outer,
                       total_depth_simulation_domain)
    # function to create varying mesh-width array in y-direction.

    # define segments of the computing grid in y-direction 
    sy1 = 0                                                 # surface
    sy2 = pipe_laying_depth - 0.5 * pipe_diameter_outer     # node above fluid node
    sy3 = pipe_laying_depth + 0.5 * pipe_diameter_outer     # node below fluid node
    sy4 = total_depth_simulation_domain                     # lower simulation boundary
    segment_array_y = [sy1, sy2, sy3, sy4]                  # array containing segment points

    # Mesh in y-direction - initializations
    dy_1 = zeros(300)
    number_of_nodes_dy_1 = 0        # to get number of nodes in segment later            
    dy_2 = zeros(2)
    dy_3 = zeros(300)
    number_of_nodes_dy_3 = 0        # to get number of nodes in segment later

    # turning points: Mesh width changes from increasing to decreasing 
    turning_point = 0           # set for if-condition in loop
    turning_point_index = 0     # set for if-condition in loop

    for i in 1:Int(floor((segment_array_y[2] - segment_array_y[1]) / min_mesh_width))

        # When Sum of Elements in dy_1 is equal to the size of this segment, the loop will be stoped.
        # Get number of nodes to cut zeros later.
        if round(sum(dy_1[1:i]); digits=2) >= round(segment_array_y[2] - segment_array_y[1]; digits=2)
            number_of_nodes_dy_1 = i - 1

            # reset variables
            turning_point = 0
            turning_point_index = 0
            break

            # first Element of dy_1 (directly under the surface) starts with minimal mesh width
        elseif i == 1
            dy_1[i] .= min_mesh_width
            difference = copy(dy_1[1])

            # as long as the depth (sum of dy_1) in the loop iteration stays smaller than the half segment size,
            # the next element of dy increases (according to the defined expansion factor)
        elseif sum(dy_1[1:i]) + min(dy_1[i - 1] * expansion_factor, max_mesh_width) <
               (segment_array_y[2] - segment_array_y[1]) / 2
            difference *= expansion_factor
            dy_1[i] .= min(difference, max_mesh_width)
            turning_point = 0

        elseif sum(dy_1[1:i]) + min(dy_1[i - 1] * expansion_factor, max_mesh_width) >=
               (segment_array_y[2] - segment_array_y[1]) / 2

            # if sum of dy_1 would be greater than the half segment size for the first time, the turning point is reached.
            if turning_point == 0
                # special dy_1-elements are calculated in between the increasing and decreasing mesh widths
                dy_1[i] .= ((segment_array_y[2] - segment_array_y[1]) - sum(dy_1[1:(i - 1)]) * 2) / node_arrangement
                dy_1[i + 1] .= copy(dy_1[i - 1])
                turning_point_index = copy(i)
                turning_point = 1

                # if turning point is already reached, the size of the next element decreases according to the defined
                # expansion expansion_factor
            elseif turning_point == 1 && i >= (turning_point_index + 2)
                difference /= expansion_factor
                dy_1[i] .= min(difference, max_mesh_width)
            end
        end
    end

    # segment 2: pipe node 
    dy_2 = [pipe_diameter_outer / 2 pipe_diameter_outer / 2]

    # segment 3: below pipe node to lower simulation boundary
    for j in 1:Int(floor((segment_array_y[4] - segment_array_y[3]) / min_mesh_width))

        # When Sum of Elements in dy_1 is equal to the depth of this segment, the loop will be stoped.
        # Get number of nodes to cut zeros later.
        if round(sum(dy_3[1:j]); digits=2) >= round(segment_array_y[4] - segment_array_y[3]; digits=2)
            number_of_nodes_dy_3 = j - 1

            # reset variables
            turning_point = 0
            turning_point_index = 0
            break

            # first Element of dy_3 (directly next to the pipe) starts with minimal mesh width
        elseif j == 1
            dy_3[j] = min_mesh_width
            difference = copy(dy_1[1])

            # as long as the depth (sum of dy_1) in the loop iteration stays smaller than the half segment size,
            # the next element of dy increases (according to the defined expansion factor)
        elseif sum(dy_3[1:j]) + min(dy_3[j - 1] * expansion_factor, max_mesh_width) <
               (segment_array_y[4] - segment_array_y[3]) / 2
            difference *= expansion_factor
            dy_3[j] = min(difference, max_mesh_width)
            turning_point = 0

        elseif sum(dy_3[1:j]) + min(dy_3[j - 1] * expansion_factor, max_mesh_width) >=
               (segment_array_y[4] - segment_array_y[3]) / 2

            # if sum of dy_1 would be greater than the half segment size for the first time, the turning point is reached.
            if turning_point == 0
                # one special dy_1-Element is calculated in between the increasing and decreasing mesh widths
                dy_3[j] = (segment_array_y[4] - segment_array_y[3]) - sum(dy_3[1:(j - 1)]) * 2
                dy_3[j + 1] = copy(dy_3[j - 1])
                turning_point = 1
                turning_point_index = copy(j)

                # if turning point is already reached, the size of the next element decreases according to the defined 
                # expansion expansion_factor
            elseif turning_point == 1 && j > (turning_point_index + 1)
                difference /= expansion_factor
                dy_3[j] = min(difference, max_mesh_width)
            end
        end
    end

    # cut zeros from dy_1 and dy_3 and initialize empty dy-vector
    dy_1_cut = view(dy_1, 1:number_of_nodes_dy_1)
    dy_3_cut = view(dy_3, 1:number_of_nodes_dy_3)
    dy = Float64[]  # Leeres Array initialisieren

    # build final dy-vector by appending
    append!(dy, dy_1_cut)  # dy_1_cut an dy anhängen
    append!(dy, dy_2)      # dy_2 an dy anhängen
    append!(dy, dy_3_cut)  # dy_3_cut an dy anhängen

    return dy
end

function create_mesh_x(min_mesh_width, max_mesh_width, expansion_factor, pipe_diameter_outer, fluid_pipe_distance)
    # generating expanding mesh in x-direction. 
    # define segments of the computing grid in x-direction 
    sx1 = 0                             # left boundary
    sx2 = 0.5 * pipe_diameter_outer       # node right next to fluid node
    sx3 = fluid_pipe_distance           # right simulation boundary
    segment_array_x = [sx1, sx2, sx3]   # array containing segment points

    number_of_nodes_dx_2 = 0            # to cut dx_2 array later

    # segment 1
    dx_1 = [0.5 * pipe_diameter_outer]
    dx_2 = zeros(100)
    # segment 2 
    for i in 1:Int(floor((segment_array_x[3] - segment_array_x[2]) / min_mesh_width))

        # When Sum of Elements in dx_1 is equal to the size of this segment, the loop will be stoped.
        # Get number of nodes to cut zeros later.
        if round(sum(dx_2[1:i]); digits=2) >= round(segment_array_x[3] - segment_array_x[2]; digits=2)
            number_of_nodes_dx_2 = i - 1
            break

            # first Element of dx_2 (right next to the fluid node) starts with minimal mesh width
        elseif i == 1
            dx_2[i] .= min_mesh_width
            difference = dx_2

            # as long as the Sum of dx_2 stays smaller than the length of segment 2,
            # the next element of dx increases (according to the defined expansion factor)
        elseif sum(dx_2[1:i]) + min(dx_2[i - 1] * expansion_factor, max_mesh_width) <
               (segment_array_x[3] - segment_array_x[2])
            difference *= expansion_factor
            dx_2[i] .= min(difference, max_mesh_width)

        elseif sum(dx_2[1:i]) + min(dx_2[i - 1] * expansion_factor, max_mesh_width) >
               (segment_array_x[3] - segment_array_x[2])
            dx_2[i] = segment_array_x[3] - sum(dx_2[:])
        end
    end

    # Leeres Array initialisieren
    # cut zeros from dx_2 and initialize empty dx-vector
    dx_2_cut = view(dx_2, 1:number_of_nodes_dx_2)
    dx = Float64[]
    # build final dx-vector by appending
    append!(dx, dx_1)
    append!(dx, dx_2_cut)
    return dx
end

# set collector properties
total_depth_simulation_domain = 10          # total depth of simulation area in y direction
pipe_laying_depth = 2                       # laying depth below surface
pipe_diameter_outer = 0.040                 # outer pipe diameter [m]
fluid_pipe_distance = 1                     # spacing between two adjacent collector pipes [m]

# define minimal mesh width - TO DO: choose between fine or rough grid
min_mesh_width = pipe_diameter_outer / 4  # minimal mesh width in x- and y-direction [- current value could be min_mesh_width for fine grid]
max_mesh_width = pipe_diameter_outer * 128
expansion_factor = 2

# accuracy mode: high, normal.
accuracy_mode = "normal"

if accuracy_mode == "normal"
    min_mesh_width = pipe_diameter_outer / 4
    max_mesh_width = pipe_diameter_outer * 32  # TO DO
    expansion_factor = 2

elseif accuracy_mode == "high"
    min_mesh_width = pipe_diameter_outer / 8
    max_mesh_width = pipe_diameter_outer * 16  # TO DO
    expansion_factor = 1.75
end

dy = create_mesh_y(min_mesh_width,
                   max_mesh_width,
                   expansion_factor,
                   pipe_laying_depth,
                   pipe_diameter_outer,
                   total_depth_simulation_domain)
dx = create_mesh_x(min_mesh_width, max_mesh_width, expansion_factor, pipe_diameter_outer, fluid_pipe_distance)

# can be deleted later, only for debugging
println(dy[:])
# println(dx[:])
