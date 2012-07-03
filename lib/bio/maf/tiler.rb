module Bio::MAF

  # Tiles a given genomic interval.
  # Inspired by: lib/bx/align/tools/tile.py in bx-python

  class Tiler

    attr_accessor :index
    attr_accessor :parser
    attr_accessor :reference
    # GenomicInterval
    attr_accessor :interval

    def species=(v)
      if v.is_a? Hash
        @species_map = v
        @species = v.keys
      else
        @species = v
      end
    end

    def tile
      parser.sequence_filter[:only_species] = @species
      # TODO: remove gaps
      blocks = index.find(parser, [interval]).sort_by { |b| b.vars[:score] }
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
      runs(mask) do |range, entry|
        if entry == :ref
          # not covered by an alignment block
          # use the reference sequence if given
        else
          # covered by an alignment block
          # TODO: slice block, accounting for gaps
          # empty: fill in
        end
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
