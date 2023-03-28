@lazy mutable struct Reference
    const genomes::Set{Genome}
    # Sequence name => (sequence, targets)
    const targets_by_name::Dict{String, Tuple{Sequence, Vector{Target}}}
    const clades::Vector{Vector{Clade{Genome}}}
    @lazy shortest_seq_len::Int
    @lazy fraction_assembled::Float64
end

function Reference()
    Reference(
        Set{Genome}(),
        Dict{String, Tuple{Sequence, Vector{Target}}}(),
        Vector{Clade{Genome}}[],
        uninit,
        uninit
    )
end

Base.show(io::IO, ::Reference) = print(io, "Reference()")
function Base.show(io::IO, ::MIME"text/plain", x::Reference)
    if get(io, :compact, false)
        show(io, x)
    else
        print(io,
            "Reference",
            "\n  Genomes:    ", ngenomes(x),
            "\n  Sequences:  ", nseqs(x),
            "\n  Ranks:      ", nranks(x),
            "\n  Seq length: ", x.shortest_seq_len,
            "\n  Assembled:  ", round(x.fraction_assembled * 100; digits=1), " %"
        )
    end
end

top_clade(x::Reference) = only(last(x.clades))
ngenomes(x::Reference) = length(x.genomes)
nseqs(x::Reference) = length(x.targets_by_name)
function nranks(x::Reference)
    isempty(x.genomes) && return 0
    length(x.clades) + 1
end

function finish!(ref::Reference)
    @isinit(ref.fraction_assembled) && return ref
    foreach(finish!, ref.genomes)
    assembly_size = genome_size = 0
    for genome in ref.genomes
        assembly_size += genome.assembly_size
        genome_size += genome.genome_size
    end
    shortest_seq_len = minimum(i -> length(first(i)), values(ref.targets_by_name); init=typemax(Int))
    shortest_seq_len == typemax(Int) && error("Cannot initialize a Reference with no sequences")
    @init! ref.shortest_seq_len = shortest_seq_len
    @init! ref.fraction_assembled = assembly_size / genome_size
    ref
end

"""
    filter_sequences(f::Function, ref::Reference)::Reference

Create a new, independent (deep copied) `Reference`, keeping only
the sequences for which `f(sequence)` is `true`.

# Examples
```julia
julia> ref
Reference
  Genomes:   1057
  Sequences: 1247324
  Ranks:     8

julia> filter_size(s -> s.length ≥ 5000, ref)
Reference
  Genomes:   1057
  Sequences: 30501
  Ranks:     8
```
"""
function filter_sequences(@nospecialize(f::Function), ref::Reference)
    ref = uninit!(deepcopy(ref))
    filter!(ref.targets_by_name) do (_, v)
        f(first(v))::Bool
    end
    for genome in ref.genomes, source in genome.sources
        filter!(i -> f(first(i))::Bool, source.sequences)
    end
    finish!(ref)
end

function filter_genomes(@nospecialize(f::Function), ref::Reference)
    ref = uninit!(deepcopy(ref))
    genomes_to_remove = Genome[]
    sources_to_remove = Set{Source}()
    filter!(ref.genomes) do g
        keep = f(g)::Bool
        keep || push!(genomes_to_remove, g)
        keep
    end
    for genome in genomes_to_remove
        union!(sources_to_remove, genome.sources)
    end
    for (_, targets) in values(ref.targets_by_name)
        filter!(targets) do (source, _)
            !in(source, sources_to_remove)
        end
    end
    for genome in genomes_to_remove
        recursively_delete_child!(genome)
    end
    for i in length(ref.clades)-1:-1:1
        empty!(ref.clades[i])
        for parent in ref.clades[i+1]
            union!(ref.clades[i], parent.children)
        end
    end
    finish!(ref)
end

function uninit!(ref::Reference)
    @uninit! ref.fraction_assembled
    @uninit! ref.shortest_seq_len
    for genome in ref.genomes
        @uninit! genome.genome_size
        @uninit! genome.assembly_size
        for source in genome.sources
            @uninit! source.assembly_size
        end
    end
    ref
end

function add_genome!(ref::Reference, genome::Genome)
    in(genome, ref.genomes) && error(lazy"Genome $(genome.name) already in reference")
    push!(ref.genomes, genome)
    ref
end

function add_sequence!(ref::Reference, seq::Sequence, targets::Vector{Target})
    for (source, span) in targets
        add_sequence!(source, seq, span)
    end
    if last(get!(ref.targets_by_name, seq.name, (seq, targets))) !== targets
        error(lazy"Duplicate sequence in reference: $(seq.name)")
    end
    ref
end

function parse_bins(
    io::IO,
    ref::Reference,
    binsplit_sep::Union{Nothing, AbstractString, Char}=nothing
)
    itr = tab_pairs(eachline(io))
    if binsplit_sep !== nothing
        itr = binsplit_tab_pairs(itr, binsplit_sep)
    end
    seqs_by_binname = Dict{SubString{String}, Vector{Sequence}}()
    for (binname, seqname) in itr
        (seq, _) = ref.targets_by_name[seqname]
        push!(get!(valtype(seqs_by_binname), seqs_by_binname, binname), seq)
    end
    [Bin(binname, seqs, ref.targets_by_name) for (binname, seqs) in seqs_by_binname]
end

const JSON_VERSION = 1
struct ReferenceJSON
    version::Int
    # [(name, flags, [(sourcename, length)])]
    genomes::Vector{Tuple{String, Int, Vector{Tuple{String, Int}}}}
    # [Sequence => sequence_length, [(subject, from, to)]]
    sequences::Vector{Tuple{String, Int, Vector{Tuple{String, Int, Int}}}}
    # [[child => parent] ...]
    taxmaps::Vector{Vector{Tuple{String, Union{String, Nothing}}}}
end
StructTypes.StructType(::Type{ReferenceJSON}) = StructTypes.Struct()

function Reference(io::IO; min_seq_length::Integer=1)
    Reference(JSON3.read(io, ReferenceJSON), Int(min_seq_length))
end

function Reference(json_struct::ReferenceJSON, min_seq_length::Int)
    if json_struct.version != JSON_VERSION
        @warn (
            "Deserializing reference JSON of version $(json_struct.version), " *
            "but the supported version of the currently loaded version of VambBenchmarks is $(JSON_VERSION)."
        )
    end
    ref = Reference()

    # Parse genomes
    for (genomename, flags, sourcesdict) in json_struct.genomes
        genome = Genome(genomename, FlagSet(UInt(flags)))
        add_genome!(ref, genome)
        for (source_name, source_length) in sourcesdict
            add_source!(genome, source_name, source_length)
        end
    end

    # Check for unique sources
    source_by_name = Dict{String, Source{Genome}}()
    for genome in ref.genomes, source in genome.sources
        if haskey(source_by_name, source.name)
            error(
                lazy"Duplicate source: \"$(source.name)\" belongs to both genome ",
                lazy"\"$(source_by_name[source.name].genome.name)\" and \"$(genome.name)\"."
            )
        end
        source_by_name[source.name] = source
    end

    # Parse sequences
    for (seq_name, seq_length, targs) in json_struct.sequences
        seq_length ≥ min_seq_length || continue
        targets = map(targs) do (source_name, from, to)
            if to < from
                error(lazy"Sequence \"$(seq_name)\" spans $(from)-$(to), must span at least 1 base.")
            end
            source = get(source_by_name, source_name, nothing)
            if source === nothing
                error(lazy"Sequence \"$(seq_name)\" maps to source \"$(source_name)\", but no such source in reference")
            end
            (source, Int(from):Int(to))
        end
        seq = Sequence(seq_name, seq_length)
        add_sequence!(ref, seq, targets)
    end

    # Add taxonomy
    copy!(ref.clades, parse_taxonomy(ref.genomes, json_struct.taxmaps))

    # Finalize the reference
    finish!(ref)
end

function save(io::IO, ref::Reference)
    json_dict = Dict{Symbol, Any}()
    json_dict[:version] = JSON_VERSION
    # Genomes
    json_dict[:genomes] = [
        (genome.name, Int(genome.flags.x), [(s.name, s.length) for s in genome.sources])
    for genome in ref.genomes]

    # Sequences
    json_dict[:sequences] = [
        (seq.name, seq.length, [(source.name, first(span), last(span)) for (source, span) in targets])
        for (_, (seq, targets)) in ref.targets_by_name
    ]

    # Taxmaps
    taxmaps = [Tuple{String, Union{String, Nothing}}[(genome.name, genome.parent.name) for genome in ref.genomes]]
    json_dict[:taxmaps] = taxmaps
    for clades in ref.clades
        v = eltype(taxmaps)()
        push!(taxmaps, v)
        for clade in clades
            parent = clade.parent
            value = parent === nothing ? nothing : parent.name
            push!(v, (clade.name, value))
        end
    end

    JSON3.write(io, json_dict)
end

function parse_taxonomy(
    genomes::Set{Genome},
    dict::Vector{Vector{Tuple{String, Union{String, Nothing}}}}
)::Vector{Vector{Clade{Genome}}}
    child_by_name = Dict{String, Node}(g.name => g for g in genomes)
    parent_by_name = empty(child_by_name)
    result = Vector{Vector{Clade{Genome}}}()
    for (rank, taxmap) in enumerate(dict)
        cladeset = Clade{Genome}[]
        for (child_name, maybe_parent_name) in taxmap
            # Get child
            child = get(child_by_name, child_name, nothing)
            child === nothing && error(
                "At rank $rank, found child name \"$child_name\", but this does not exist " *
                "on the previous rank."
            )
            # Create parent if it does not already exist
            parent_name = maybe_parent_name === nothing ? child_name : maybe_parent_name
            parent = get(parent_by_name, parent_name, nothing)::Union{Clade{Genome}, Nothing}
            if parent === nothing
                parent = Clade(parent_name, child)
                parent_by_name[parent_name] = parent
                push!(cladeset, parent)
            else
                add_child!(parent, child)
                child.parent = parent::Clade{Genome}
            end
        end
        # Enforce all childrens have parents
        for child in values(child_by_name)
            isdefined(child, :parent) || error(
                "At rank $rank, child $(child.name) has no parent"
            )
        end
        parent_by_name, child_by_name = child_by_name, parent_by_name
        empty!(parent_by_name)
        push!(result, cladeset)
    end
    # If there isn't exactly one remaining clade at the top, make one.
    if length(child_by_name) != 1 || only(values(child_by_name)) isa Genome
        top = Clade("top", first(values(child_by_name)))
        for child in values(child_by_name)
            child === first(top.children) || add_child!(top, child)
        end
        push!(result, [top])
    else
        top = only(values(child_by_name))::Clade{Genome}
        # If multiple redundant top clades, go down to the useful level
        while nchildren(top) == 1 && only(top.children).name == top.name
            top = only(top.children)::Clade{Genome}
            top.parent = nothing
            pop!(result)
        end
    end
    @assert length(last(result)) == 1
    foreach(i -> sort!(by=j -> j.name, i), result)
    return result
end
