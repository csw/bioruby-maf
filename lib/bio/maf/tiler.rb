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

    def initialize
      @species_map = {}
    end

    def ref_data(range)
      if reference
        if reference.respond_to? :read_interval
          reference.read_interval(range.begin, range.end)
        elsif reference.is_a? String
          reference.slice(range)
        else
          raise "Unhandled reference data source: #{reference}"
        end
      else
        nil
      end
    end

    def tile
      parser.sequence_filter[:only_species] = @species
      # TODO: remove gaps
      blocks = index.find([interval], parser).sort_by { |b| b.vars[:score] }
      mask = Array.new(interval.length, :ref)
      i_start = interval.zero_start
      i_end = interval.zero_end
      if reference
        ref_region = ref_data(i_start...i_end)
      end
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
          text[0] << if ref_region
                       ref_region.slice(range)
                     else
                       'N' * range_size
                     end
          stars = '*' * range_size
          nonref_text.each { |t| t << stars }
        else
          # covered by an alignment block
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
      species.zip(tile()) do |species, text|
        sp_out = species_map[species] || species
        f.puts "> #{sp_out}"
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

  class FASTARangeReader
    attr_reader :f

    def initialize(fspec)
      if fspec.respond_to? :seek
        @f = fspec
      else
        @f = File.open(fspec)
      end
    end

    def read_interval(z_start, z_end)
      data = ''
      f.seek(0)
      first = f.readline
      raise "expected FASTA comment" unless first =~ /^>/
      region_size = z_end - z_start
      in_region = false
      pos = 0
      f.each_line do |line_raw|
        raise "unexpected line start" unless line_raw =~ /^\w/
        line = line_raw.strip
        end_pos = pos + line.size
        if (! in_region) && pos <= z_start && z_start < end_pos
          data << line.slice((z_start - pos)...(line.size))
          in_region = true
        elsif in_region
          need = region_size - data.size
          if need > line.size
            data << line
          else
            # last line
            data << line.slice(0, need)
            break
          end
        end
        pos = end_pos
      end
      return data
    end
  end
end
