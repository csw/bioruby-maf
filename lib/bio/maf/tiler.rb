module Bio::MAF

  # Tiles a given genomic interval.
  # Inspired by: lib/bx/align/tools/tile.py in bx-python

  class Tiler

    attr_accessor :index
    attr_accessor :parser
    attr_accessor :reference
    # GenomicInterval
    attr_accessor :interval
    attr_accessor :species
    attr_accessor :species_map

    def tile
      parser.sequence_filter[:only_species] = @species
      # TODO: remove gaps
      blocks = index.find([interval], parser).sort_by { |b| b.vars[:score] }
      mask = Array.new(interval.length, :ref)
      i_start = interval.zero_start
      i_end = interval.zero_end
      blocks.each do |block|
        ref = block.ref_seq
        slice_start = [i_start, ref.start].max
        slice_end = [i_end, ref.end].min
        mask.fill(block,
                  (slice_start - i_start)...(slice_end - i_start))
      end
      text = []
      species.each { |s| text << '' }
      nonref_text = text[1...text.size]
      runs(mask) do |range, block|
        g_range = (range.begin + i_start)...(range.end + i_start)
        if block == :ref
          # not covered by an alignment block
          # use the reference sequence if given, otherwise 'N'
          range_size = range.end - range.begin
          text[0] << if reference
                       reference.slice(g_range)
                     else
                       'N' * range_size
                     end
          stars = '*' * range_size
          nonref_text.each { |t| t << stars }
        else
          # covered by an alignment block
          # TODO: slice block, accounting for gaps
          # for empty sequences: fill in
          t_range = block.ref_seq.text_range(g_range)
          species.each_with_index do |species, i|
            sp_text = text[i]
            seq = block.sequences.find { |s| s.source == species || s.species == species }
            if seq
              # got alignment text
              sp_text << seq.text.slice(t_range)
            else
              # no alignment for this one here, use '*'
              sp_text << '*' * (t_range.end - t_range.begin)
            end
          end
        end
      end
      text
    end

    def write_fasta(f)
      species.zip(tile()) do |seq, text|
        f.puts "> #{seq}"
        f.puts text
      end
    end

    def runs(mask)
      cur = nil
      cur_start = nil
      mask.each_with_index do |obj, i|
        if ! cur.equal?(obj)
          yield(cur_start...i, cur) if cur
          cur = obj
          cur_start = i
        end
      end
      yield(cur_start...mask.size, cur)
    end

  end
end
