module GenbankParse

using Bio.Seq,
      Bio.Intervals,
      DataFrames,
      JLD

function findfeatures(s::IOStream, startindex::Int)
  seek(s, startindex)
  while !eof(s)
      p = position(s)
      l = readline(s)
      if ismatch(r"^FEATURES", l)
          return p
      end
  end
  error("No features list found after position $startindex")
end


function findlocus(s::IOStream, startindex::Int)
  seek(s, startindex)
  while !eof(s)
      p = position(s)
      l = readline(s)
      if ismatch(r"^LOCUS", l)
          return p
      end
  end
  error("No locus found after position $startindex")
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
  error("No origin found after position $startindex")
end


function findcontigend(s::IOStream, startindex::Int)
  seek(s, startindex)
  while !eof(s)
      l = readline(s)
      if l == "//\n"
          return position(s)
      end
  end
  error("No contigs found after position $startindex")
end


function getinfo(s::IOStream, startindex::Int)
  md = Dict()
  featuresstart = findfeatures(s, startindex)
  seek(s, startindex)
  contiginfo = String(read(s, featuresstart - startindex))

  locusline = match(
      r"LOCUS\s{7}(\S+)\s+\d+ bp\s+(\w+)\s+(\w+)\s+([A-Z]{3})(?:\s+(\d\d\-[A-Z]{3}\-\d{4}))?",
      contiginfo
      )
  definition = replace(match(r"(?<=\n)DEFINITION  ([\s\S]+?)(?=\n[A-Z])", contiginfo).captures[1], r"\n\s+", " ")
  accession = match(r"(?<=\n)ACCESSION   ([\S]+)", contiginfo).captures[1]

  source = match(r"(?<=\n)SOURCE      ([\S\s]+?)\n", contiginfo)
      if source != nothing
          source = source.captures[1]
      end
  organism = match(r"(?<=\n)  ORGANISM  ([\S\s]+?)\n", contiginfo)
      if organism != nothing
          organism = organism.captures[1]
      end

  return locusline.captures[1], Dict([
                  ("moltype", locusline.captures[2]),
                  ("topology", locusline.captures[3]),
                  ("gbk division", locusline.captures[4]),
                  ("definition", definition),
                  ("souce", source),
                  ("organism", organism)
                  ])
end


function getfeatures(s::IOStream, startindex::Int)
    findfeatures(s, startindex)
    function _iter()
        featstart = r"^\s{5}(\w+)\s+([complent(\d.<>]+)"
        tagstart = r"^\s{21}\/(\w+)=(.+)"
        tagcont = r"^\s{21}(.+)"

        pos = position(s)
        ln = readline(s)

        @label newfeature
        feature = Dict()
        feat = match(featstart, ln)
        if feat == nothing
            println("exit")
            @goto ex
        else
            feature["type"] = feat.captures[1]
            loc = match(r"(\d+)..(\d+)", feat.captures[2])
            feature["location"] = Dict{AbstractString, Int}([
                        ("start", parse(Int, loc.captures[1])),
                        ("end", parse(Int, loc.captures[2]))
                        ])
            if startswith(feat.captures[2], "complement")
                feature["location"]["strand"] = -1
            else
                feature["location"]["strand"] = 1
            end
        end


        tagtype = nothing
        @label tags

        ln = readline(s)

        if !ismatch(r"^\s{21}", ln)
            @goto addfeature
        end

        nexttag = match(tagstart, ln)
        sametag = match(tagcont, ln)

        if nexttag != nothing
            tagtype = nexttag.captures[1]
            feature[tagtype] = nexttag.captures[2]
        elseif sametag != nothing
            feature[tagtype] = "$(feature[tagtype]) $(sametag.captures[1])"
        end

        if tagtype == "translation"
            feature[tagtype] = AminoAcidSequence(replace(feature[tagtype], r"[\"\s]", ""))
        else
            feature[tagtype] = replace(feature[tagtype], "\"", "")
        end
        @goto tags

        @label addfeature
        produce(feature)

        @goto newfeature

        @label ex
    end
    Task(_iter)
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

end # module GenbankParse
