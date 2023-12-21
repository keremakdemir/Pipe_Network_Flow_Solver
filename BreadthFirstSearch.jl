using DataStructures

#= 
bfs_cycles function applies breadth-first search to find loops (i.e., cycles) in a given network.
It systematically explores the edges of an adjacency matrix to discover every vertex that is 
reachable from a distinguished source vertex s. This works for both directed and undirected
graphs. This algorithm computes the distance (smallest number of edges) from vertex s to each 
reachable vertex. To keep track of progress, this algorithm gives colors to each vertex 
as white, gray, or black. All vertices start out white and later become gray and then black.

Inputs:
adj: Adjacency dictionary that shows which vertex is connected to which vertex
s: Source vertex to initiate the algorithm
ignore (optional): This is a set of a vector of vertex points to ignore during breadth-first search

Returns:
d: Dictionary containing the distances from source vertex to specific vertex
π: Dictionary showing the immediate predecessor of a specific vertex
cycle_set: Set containing edge tuples that are not used (or explored) during breadth-first search algorithm
=#
@enum Color white gray black
function bfs_cycles(adj::Dict{T, Vector{T}}, s, ignore=Set()) where T <: Any

    cycle_set = Set{Tuple{T, T}}()
    color = Dict{T, Color}()
    d = Dict{T, Int}()
    π = Dict{T, T}()

    for u in keys(adj)
        color[u] = white
        d[u] = typemax(Int)
    end
    color[s] = gray
    d[s] = 0

    Q = Queue{T}()
    enqueue!(Q, s)

    while !isempty(Q)
        u = dequeue!(Q)
        for v in adj[u]
            if (u, v) ∉ ignore
                if color[v] == white
                    color[v] = gray
                    d[v] = d[u] + 1
                    π[v] = u
                    enqueue!(Q, v)
                elseif v != π[u]
                    push!(cycle_set, (u, v))
                end
            end
        end
        color[u] = black
    end

    d, π, cycle_set
end


#= 
get_path function utilizes breadth-first search algorithm results and finds the closest path between
two points (vertices)

Inputs:
π: Dictionary showing the immediate predecessor of a specific vertex. This comes from bfs_cycles function.
s: Starting point (vertex)
v: Ending point (vertex)

Returns:
Vector containing closest path from vertex s to vertex v
=#
function get_path(π, s, v)
    path = list(v)
    while s != v
        v = π[v]
        path = cons(v, path)
    end
    collect(path)
end


#= 
bfs_results function utilizes breadth-first search algorithm results and illustrates the distances of
vertices from a source vertex.

Inputs:
adj: Adjacency dictionary that shows which vertex is connected to which vertex
d: Dictionary containing the distances from source vertex to specific vertex
π: Dictionary showing the immediate predecessor of a specific vertex

Returns:
Prints the the distances of vertices from a source vertex.
=#
function bfs_results(adj, d, π)
    println("\nu, π, d\n-------")
    for u in keys(adj)
        println("$u, $(get(π, u, '-')), $(d[u])")
    end
end


#= 
bfs_path function utilizes breadth-first search algorithm and shows the path from a source vertex to
destination vertex. This can be used to in conjuction with cycle_set to find loops (i.e., cycles) in
a given network. The user should provide an ignore input so that a valid cycle can be obtaioned.

Inputs:
adj: Adjacency dictionary that shows which vertex is connected to which vertex
u: Source vertex
v: Ending vertex
ignore: This is a set of a vector of vertex points to ignore during breadth-first search

Returns:
Vector containing closest path from vertex u to vertex v 
=#
function bfs_path(adj, u, v, ignore)
    d, π, cycle_set = bfs_cycles(adj, u, ignore)
    push!(get_path(π, u, v), u)
end

