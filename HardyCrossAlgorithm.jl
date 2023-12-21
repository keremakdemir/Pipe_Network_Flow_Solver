#= 
This algorithm analyzes the flow in networks of conduits by utilizing Hardy
Cross method. Two principles should be satisfied:

- Continuity of Flow: Total flow going into any junction equals to the total 
flow leaving the junction (conservation of mass)
- Continuity of Potential: Total change in potential along any closed path is 
zero (conservation of energy)

Considering these, either flows in individual conductors or potentials at the
junction points should be taken as unknowns. 

This algorithm takes potentials as unknows and uses method of balancing heads. 
In this sense, it is assumed that flow quantities at inlets and outlets are 
known (or given an initial value). 

The algorithm works in successive corrections (iterations). These steps are 
shown below:

1. Assume any distribution of flow, which satisfies continuity of flow.
2. Compute in each pipe the loss of head.
3. Set up in each circuit a counterbalancing flow to balance the head
in that circuit.
4. Compute the revised flows and repeat the procedure.

Loss of head in any given pipe is assumed to be h1 = r*Q^n, where
r: Loss of head in the pipe for unit quantity of flow (resistance)
Q: Quantity of flow
n: Constant exponent

A secondary equation is used in each closed circuit to set up a
counterbalancing flow. This equation is h2 = n*r*Q^(n-1)

counterbalancing flow in each closed circuit is calculated as follows:

Delta = sum of h1 (reference to direction of flow) /
        sum of h2 (without reference to direction of flow)

By utilizing Delta, revised flows in each closed circuit (loop) is calculated
and process goes on until Delta for each loop becomes zero with a tolerance of
10^-2.
=#

##############################################################################
##############################################################################
##############################################################################

#= 
Defining a network structure to designate any given network with nodes
and edges. This is an immutable structure.

Inputs:
node: Dictionary of nodes
edge: Dictionary of dictionaries of edges
=#
struct Network
    node :: Dict
    edge :: Dict  
end


#= 
create_network function initializes an empty network from a 
collection (like arrays/vectors) of node labels. This function
provides an initial 0 flow value for each node and creates 
an empty dictionary of dictionaries to store the edges

Inputs:
nodes: Vector of node labels

Returns:
A network object containing node labels and edges
=#
function create_network(nodes::Vector)
    Network(Dict(i => Node(0) for i in nodes),
    Dict(i => Dict() for i in nodes))
end


#= 
Defining Node structure that store flow values. This is an 
mutable structure so that nodal flow can be updated within
iterations. 

Inputs:
flow: Flow in or out of each node
=#
mutable struct Node
    flow
end


#= 
nodal_flow function returns the flow at any node within a specificed network. 
By convention, a flow into a node is assumed to be positive (flow entering 
the network), and a flow out of a node is assumed to be negative (flow 
leaving the network).

Inputs:
net: Network to be analyzed
i: Node label within the network

Returns:
Flow at any node within the network
=#
function nodal_flow(net::Network, i)
    net.node[i].flow
end


#= 
set_flow! function alters the inflow or outflow of a node to a scalar value q

Inputs:
net: Network to be analyzed
i: Node label within the network
q: New flow value
=#
function set_flow!(net::Network, i, q)
    net.node[i].flow = q
end


#= 
Defining Edge structure that store sign of flow and pipe information

Inputs:
sign: Sign of edges (+1 or -1 depending on assumed flow direction)
pipe: Refers to a pipe object. This is done to create forward and reverse
edges. Forward and reverse edges will share the same pipe but with a different
sign.
=#
struct Edge
    sign                  
    pipe
end


#= 
Defining Pipe structure that store flow, resistance and exponent information.
This is an mutable structure so that pipe flow can be updated within iterations. 

Inputs:
flow: Initial or calculated flow within a pipe
resistance: Resistance of the pipe
exponent: Exponent to use in Hardy Cross method.
=#
mutable struct Pipe
    flow
    resistance
    exponent
end


#= 
add_pipe! function alters the network edges by creating pipe objects
within a network between nodes i and j. Flow from i to j is assumed to
be positive. This creates a symmetric edge set (both forward and reverse edges).
The pipe object itself, which maintains pipe properties like flow and resistance, 
is shared by the forward and reverse edges between any node pair.

Inputs:
net: Network to be analyzed
i: Starting node of edge
j: Ending node of edge
p: Refers to the pipe object for forward and reverse edges
=#
function add_pipe!(net::Network, i, j, p::Pipe)
    net.edge[i][j] = Edge(+1, p)
    net.edge[j][i] = Edge(-1, p)
end


#= 
flow function returns the flow of an pipe within a network

Inputs:
net: Network to be analyzed
i: Starting node of edge
j: Ending node of edge

Returns:
Flow within a pipe
=#
function flow(net::Network, i, j)
    net.edge[i][j].pipe.flow
end


#= 
resistance function returns the resistance of an pipe within a network

Inputs:
net: Network to be analyzed
i: Starting node of edge
j: Ending node of edge

Returns:
Resistance of a pipe
=#
function resistance(net::Network, i, j)
    net.edge[i][j].pipe.resistance
end


#= 
exponent function returns the exponent of an pipe within a network

Inputs:
net: Network to be analyzed
i: Starting node of edge
j: Ending node of edge

Returns:
Hardy Cross exponent of a pipe
=#
function exponent(net::Network, i, j)
    net.edge[i][j].pipe.exponent
end


#= 
sign function returns the sign of an edge within a network

Inputs:
net: Network to be analyzed
i: Starting node of edge
j: Ending node of edge

Returns:
Sign of the edge
=#
function sign(net::Network, i, j)
    net.edge[i][j].sign
end


#= 
net_flow function calculates the net flow of a node within a network

Inputs:
net: Network to be analyzed
i: Node label within the network to calculate net flow

Returns:
Net flow of a node within specified network
=#
function net_flow(net::Network, i)
    nodal_flow(net, i) - sum([(e.sign * e.pipe.flow) for (j, e) in net.edge[i]])
end


#= 
inc_flow! function alters the flow by an amount q in edge (i, j) in a network

Inputs:
net: Network to be analyzed
i: Starting node of edge
j: Ending node of edge
q: Amount of flow to be increased/decreased by
=#
function inc_flow!(net::Network, i, j, q)
    net.edge[i][j].pipe.flow += sign(net, i, j) * q
end


#= 
inc_flow_loop! function alters flow by an amount q in a sequence of edges 
(which should form a cycle)

Inputs:
net: Network to be analyzed
edges: Vector of edges that represent a loop (like [(1, 2), (2, 4), (4, 1)])
q: Amount of flow to be increased/decreased by
=#
function inc_flow_loop!(net::Network, edges::Vector, q)
    for (i, j) in edges
        inc_flow!(net, i, j, q)
    end
end


#= 
print_network function displays the network elements in a human readable format

Inputs:
net: Network to be displayed
=#
function print_network(net::Network)
    println("Nodes:")
    for (i, n) in net.node
        println("Node $i => Nodal Flow: $(n.flow), Net Flow: $(net_flow(net, i))")
    end
    println("\nEdges:")
    for (i, d) in net.edge, j in keys(d)
        if sign(net, i, j) == +1
            println("From node $i to node $j => Flow: $(flow(net, i, j)), Resistance: $(resistance(net, i, j))")
        end
    end
end


#= 
pipe_loss_terms function calculates head loss quantities (numerator 
and denominator) for pipe (i, j). Numerator considers direction of the flow
whereas denominator does not consider direction of the flow.

Inputs:
net: Network to be analyzed
i: Starting node of pipe
j: Ending node of pipe

Returns:
A tuple containing head loss quantities (numerator and denominator) 
for pipe (i, j)
=#
function pipe_loss_terms(net::Network, i, j)
    num = flow(net, i, j)*(abs(flow(net, i, j))^(exponent(net, i, j)-1))*resistance(net, i, j)*sign(net, i, j)
    denom = exponent(net, i, j)*(abs(flow(net, i, j))^(exponent(net, i, j)-1))*resistance(net, i, j)

    (num, denom)
end

function loss_terms(net::Network, i, j)
    Q = flow(net, i, j)
    r = resistance(net, i, j)
    n = exponent(net, i, j)
    # see eq. 10.7.2 in Streeter on p. 428:
    num = 
    den = abs(Q)^(n-1)
    num, den
end

#= 
ΔQ function calculates counterbalancing flow to balance the head in 
that loop.

Inputs:
net: Network to be analyzed
edges: Vector of edges that represent a loop (like [(1, 2), (2, 4), (4, 1)])

Returns:
Counterbalancing flow that consist of sum of numerator divided by sum of denominators
of pipe loss terms
=#
function ΔQ(net::Network, edges::Vector)
    total_num = 0
    total_denom = 0
    for (i, j) in edges
            total_num += pipe_loss_terms(net, i, j)[1]
            total_denom += pipe_loss_terms(net, i, j)[2]
    end

    -total_num/total_denom
end


#= 
sum_pairs function calculates the sum of starting node (i) and ending node (j) 
of a vector of edges. For example, sum_pairs([(1, 2), (2, 4), (4, 1)]) returns [7, 7].
If first and second element of returned array is the same, then the vector of edges
constitute a loop.

Inputs:
vector_of_pairs: Vector of edges that may or may not represent a loop (like [(1, 2), (2, 4), (4, 1)])

Returns:
Sum of i and j's of the vector of edges
=#
function sum_pairs(vector_of_pairs)
    total_i = 0
    total_j = 0
    for (i, j) in vector_of_pairs
            total_i += i
            total_j += j
    end

    [total_i, total_j]
end


#= 
HC_step function does the necessary head loss calculations for a single iteration and compute the ΔQ
each loop. After that, it calculates the revised flows. 

Inputs:
net: Network to be analyzed
loops: Vector containing all loops in the network

Returns:
Delta_loops: A vector containing ΔQ for each loop
=#
function HC_step(net::Network, loops::Vector)

    Delta_loops = []
    for z in 1:length(loops)
            Delta_sp = ΔQ(net,loops[z])
            push!(Delta_loops, Delta_sp)
            inc_flow_loop!(net, loops[z], Delta_sp)
            
            #Ensuring the continuity of flow for each node (net flow should be 0)
            for (i, n) in net.node
                    @assert isapprox(net_flow(net, i), 0, atol= 10^-2) "Continuity of flow at one \
                    or more nodes are violated."
            end
    end
    Delta_loops
end


#= 
solve_HardyCross function takes a network and vector of loops and runs Hardy Cross algorithm 
until ΔQ for each loop becomes zero with a tolerance of 10^-2. It initializes ΔQ values of 
each loop as 1.

Inputs:
net: Network to be analyzed
loops: Vector containing all loops in the network
=#
function solve_HardyCross(net::Network, loops::Vector)
    
    #Ensuring that provided loops are actually valid loops
    for z in 1:length(loops)
            #Selecting a loop
            sp_loop = loops[z]
            #Checking if if loop consist of at least 3 edges
            @assert length(sp_loop) >= 3 "A loop should consist of at least 3 edges."
            #Checking if a self loop is observed or not
            @assert length([false for i in sp_loop if i[1]==i[2]]) == 0 "Self loops are prohibited (i and j \
            of edges should be different)."
            #Making sure that loops end at the node they started
            @assert sp_loop[1][1] == sp_loop[end][2] "Loops should end at the node they started."
            #Making sure that there are no self returning edges in the loops
            @assert length([false for i in sp_loop if (i[2], i[1]) in sp_loop]) == 0 "Self returning edges are \
            prohibited (e.g., (2,3), (3,2) cannot exist together in a loop)."
            #Summing i and j pairs and checking if they are equal
            loop_sum = sum_pairs(sp_loop)
            @assert loop_sum[1]==loop_sum[2] "Sum of i and j of edges are not equal to each other."
    end

    delta_vals = ones(length(loops))
    steps = 1

    while any(x->!isapprox(x, 0, atol= 10^-10), delta_vals)
            delta_vals = HC_step(net, loops)
            steps += 1
    end
    println("---------- Iteration $steps ----------\n")
    for q in 1:length(delta_vals)
        println("ΔQ$q = $(delta_vals[q])")
    end
    print_network(net)
end


#= 
find_net_loops function utilizes breadth-first search algorithm and find valid loops in the
network.

Inputs:
net: Network to be analyzed
wn (optional): Index of node to use as source vertex (default = 1)

Returns:
loops: Vector of vectors that shows every loop in the network
=#
function find_net_loops(net::Network, wn=1)
    #Creating an empty vector to store loops determined by breadth-first search algorithm
    loops = []

    #Creating an adjacency dictionary that shows which node is connected to which node
    adj_dict = Dict(i => Vector{typeof(i)}() for i in keys(net.node))
    for n in keys(net.node)
        adj_dict[n] = collect(keys(net.edge[n]))
    end

    #Running breadth-first search algorithm with a source node
    d, p, cs = bfs_cycles(adj_dict, collect(keys(net.node))[wn])

    #Creating an empty vector to store encountered edge tuples 
    enc_edge_pairs = []
    for edg in cs
        i, j = edg[1], edg[2]
        if (j, i) ∉ enc_edge_pairs

            #Appending edge tuple to enc_edge_pairs vector to so that we don't double count loops
            push!(enc_edge_pairs, edg)

            #Finding the loop, and creating a vector of tuples to match our loop definition
            my_path = bfs_path(adj_dict, i, j, cs)
            my_loop = [(my_path[no], my_path[no+1]) for no in collect(1:(length(my_path)-1))]
            #Appending the loops to loops vector 
            push!(loops, my_loop)
        end
    end

    loops
end


#= 
solve_HardyCross_auto function takes a network, determines the loops and runs Hardy Cross algorithm 
until ΔQ for each loop becomes zero with a tolerance of 10^-2. It initializes ΔQ values of 
each loop as 1.

Inputs:
net: Network to be analyzed
sv (optional): Index of node to use as source vertex during breadth-first search algorithm (default = 1)
=#
function solve_HardyCross_auto(net::Network, sv=1)

    #Finding loops with find_net_loops() function
    loops = find_net_loops(net, sv)
        
    #Ensuring that provided loops are actually valid loops
    for z in 1:length(loops)
            #Selecting a loop
            sp_loop = loops[z]
            #Checking if if loop consist of at least 3 edges
            @assert length(sp_loop) >= 3 "A loop should consist of at least 3 edges."
            #Checking if a self loop is observed or not
            @assert length([false for i in sp_loop if i[1]==i[2]]) == 0 "Self loops are prohibited (i and j \
            of edges should be different)."
            #Making sure that loops end at the node they started
            @assert sp_loop[1][1] == sp_loop[end][2] "Loops should end at the node they started."
            #Making sure that there are no self returning edges in the loops
            @assert length([false for i in sp_loop if (i[2], i[1]) in sp_loop]) == 0 "Self returning edges are \
            prohibited (e.g., (2,3), (3,2) cannot exist together in a loop)."
            #Summing i and j pairs and checking if they are equal
            loop_sum = sum_pairs(sp_loop)
            @assert loop_sum[1]==loop_sum[2] "Sum of i and j of edges are not equal to each other."
    end

    delta_vals = ones(length(loops))
    steps = 1

    while any(x->!isapprox(x, 0, atol= 10^-10), delta_vals)
            delta_vals = HC_step(net, loops)
            steps += 1
    end
    println("---------- Iteration $steps ----------\n")
    for q in 1:length(delta_vals)
        println("ΔQ$q = $(delta_vals[q])")
    end
    print_network(net)
end

