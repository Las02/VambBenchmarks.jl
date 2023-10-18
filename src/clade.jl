# Type parameter G here is always Genome - I only make it parametric so I can have
# mutually recursive types b/w Clade and Genome
"""
    Clade{Genome}(name::AbstractString, child::Union{Clade{Genome}, Genome})

A `Clade` represents any clade above `Genome`. Every `Genome` is expected to belong
to the same number of clades, e.g. there may be exactly 7 levels of clades above every `Genome`.
`Clade`s always have at least one child (which is either a `Genome` or a `Clade` one rank lower),
and a parent, unless it's the unique top clade from which all other clades and genomes descend from.
The rank of a `Genome` is 0, clades that contain genomes have rank 1, and clades containing rank-1
clades have rank 2 etc.
By default, zero-indexed ranks correspond to OTU, species, genus, family, order, class, phylum and domain.

# Examples
```
julia> top_clade(ref)
Genus "F", 3 genomes
├─ Species "D", 2 genomes
│  ├─ Genome(gA)
│  └─ Genome(gB)
└─ Species "E", 1 genome
   └─ Genome(gC)

julia> top_clade(ref).children
2-element Vector{Clade{Genome}}:
 Species "D", 2 genomes
 Species "E", 1 genome
```
"""
mutable struct Clade{G}
    name::String
    rank::Int
    ngenomes::Int
    parent::Union{Clade{G}, Nothing}
    children::Union{Vector{Clade{G}}, Vector{G}}

    function Clade(name::String, child::Union{Clade{G}, G}) where {G}
        (rank, ngenomes) = if child isa G
            (@isinit(child.parent)) &&
                existing_parent_error(name, child.name, child.parent.name)
            (1, 1)
        else
            parent = child.parent
            parent === nothing || existing_parent_error(name, child.name, parent.name)
            (child.rank + 1, child.ngenomes)
        end
        instance = new{G}(name, rank, ngenomes, nothing, [child])
        child.parent = instance
        return instance
    end
end

@noinline function existing_parent_error(child_name, parent_name, other_parent_name)
    error(
        "Attempted to add parent \"$parent_name\" to child \"$child_name\", which already has parent \"$other_parent_name\"",
    )
end

const RANKS =
    ["strain", "species", "genus", "family", "order", "class", "phylum", "domain", "top"]

const RANK_BY_NAME = Dict(rank => i - 1 for (i, rank) in enumerate(RANKS))

function Base.show(io::IO, x::Clade)
    suffix = x.ngenomes == 1 ? "" : "s"
    print(
        io,
        titlecase(RANKS[x.rank + 1]),
        " \"",
        x.name,
        "\", ",
        x.ngenomes,
        " genome",
        suffix,
    )
end

function Base.show(io::IO, ::MIME"text/plain", x::Clade)
    if get(io, :compact, false)
        show(io, x)
    else
        buf = IOBuffer()
        AbstractTrees.print_tree(buf, x; maxdepth=3)
        seekstart(buf)
        for (i, line) in zip(1:25, eachline(buf))
            println(io, line)
            i == 25 && print(io, '⋮')
        end
    end
end

AbstractTrees.children(x::Clade) = x.children
AbstractTrees.parent(x::Clade) = x.parent
AbstractTrees.treebreadth(x::Clade) = x.ngenomes
nchildren(x::Clade) = length(x.children)
istop(x::Clade) = isnothing(x.parent)
