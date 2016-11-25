module GenbankParse

using Bio.Seq,
      Bio.Intervals,
      DataFrames,
      JLD,
      Logging

@Logging.configure(level=DEBUG)

type GenbankFeature
    start::Int
    stop::Int
    strand::Int
    featuretype::AbstractString
    tags::Dict

    function GenbankFeature(start::Int, stop::Int, strand::Int, featuretype::AbstractString, tags::Dict)
        start < stop ? true : error("start must be before stop")
        strand == 1 || strand == -1 ? true : error("strand must be 1 or -1")
        new(start, stop, strand, featuretype, tags)
    end
end


type GenbankLocus
    name::AbstractString
    len::Int
    moltype::AbstractString
    topo::Nullable{AbstractString}
    gbkdivision::Nullable{AbstractString}
    date::Nullable{Date}

    function GenbankLocus(name::AbstractString, len::Int, moltype::AbstractString, topo=nothing, gbkdivision=nothing, date=nothing)
        moltype == "DNA" || moltype == "RNA" ? true : @warn("Invalid molecule type: $moltype")
        topo == nothing || topo == "linear" || topo == "circular" ? true : @warn("Invalid topology: $topo")
        gbkdivision = Set(["PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "VRL", "PHG", "SYN", "UNA", "EST", "PAT", "STS", "GSS", "HTG", "HTC", "ENV"])
        gbkdivision == nothing || gbkdivision in gbkdivision ? true : @warn("Invalid GenBank Division: $gbkdivision")
        new(name, len, moltype, Nullable{AbstractString}(topo), Nullable{AbstractString}(gbkdivision), Nullable{Date}(date))
    end
end


type GenbankMetadata
    locus::GenbankLocus
    definition::AbstractString
    accession::AbstractString
    organism::AbstractString
    features::Vector{GenbankFeature}
end


function parsefeatures(s::IOStream, line::AbstractString)
    startswith(line, "FEATURES") ? true : @err("Not the beginning of FEATURES")
    function _iter()
        @debug("Parse Features iteration started")
        line = readline(s)

        while ismatch(r"^\s{5}\w+", line)
            @debug(" Feature typeLine: $line")
            currentfeature = Vector{AbstractString}([line])
            line = readline(s)
            while ismatch(r"^\s{21}", line)
                push!(currentfeature, line)
                line = readline(s)
            end

            feature = capturefeature(currentfeature)
            @debug("Producing Feature: $feature")
            produce((feature, line))

        end

        @debug("DEBUG: out of features, line = $line")
        if ismatch(r"^\w", line)
            @info("End of Features")
        else
            @err("Something went wrong")
        end
    end
    Task(_iter)
end


function capturefeature{T <: AbstractString}(featurelines::Vector{T})
    featstart = match(r"^\s{5}(\w+)\s+(complement\()?<?(\d+)\.{2}(\d+)>?\)?", featurelines[1])
    featstart.captures[1] != nothing ? tp = featstart.captures[1] : @err("Not a valid feature")
    featstart.captures[2] == nothing ? strand = 1 : strand = -1
    start = parse(Int64, featstart.captures[3])
    stop = parse(Int64, featstart.captures[4])
    @info("Feature: $tp\nLocation: $start-$stop on strand $strand")
    feature = GenbankFeature(start, stop, strand, tp, Dict())

    tag = nothing

    for line in featurelines[2:length(featurelines)]
        parser = match(r"\s{21}(\/(\w+)=)?(.+)", line)

        if parser != nothing && parser.captures[1] != nothing
            tag = parser.captures[2]

            tag in keys(feature.tags) ? @warn("WARNING: Overwriting $tag tag in feature at $start") : true
            feature.tags[tag] = parser.captures[3]
        elseif parser != nothing
            feature.tags[tag] = "$(feature.tags[tag]) $(parser.captures[2])"
        else
            @warn("empty tag line")
        end
    end

    for key in keys(feature.tags)
        @debug("key: $key")
        if key == "translation"
            feature.tags[key] = AminoAcidSequence(replace(feature.tags[key], r"[\"\s]", ""))
        else
            feature.tags[key] = replace(feature.tags[key], "\"", "")
        end
    end
    return feature
end


function findorigin(s::IOStream, startindex::Int)
  seek(s, startindex)
  while !eof(s)
      p = position(s)
      l = readline(s)
      if ismatch(r"^ORIGIN", l)
          return p
      end
  end
  @warn("No origin found after position $startindex")
end


function findcontigend(s::IOStream, startindex::Int)
  seek(s, startindex)
  while !eof(s)
      l = readline(s)
      if l == "//\n"
          return position(s)
      end
  end
  @warn("No contigs found after position $startindex")
end


function getseq(s::IOStream, startindex::Int)
  o = findorigin(s, startindex)
  seqstart = position(s)
  seqend = findcontigend(s, o)
  return parseseq(s, seqstart, seqend)
end


function parseseq(s::IOStream, seqstart::Int, seqend::Int)
  seek(s, seqstart)
  seq = String(read(s, (seqend - seqstart - 3)))
  return DNASequence(replace(seq, r"[ \n\d]+", ""))
end


function parsegbk(s::IOStream)
    seek(s, 0)
end


end # module GenbankParse
