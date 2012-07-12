require 'pathname'
require 'zlib'

module Bio::MAF

  # Tiles a given genomic interval.
  # Inspired by: lib/bx/align/tools/tile.py in bx-python

  class Tiler

    attr_accessor :index
    attr_accessor :parser
    attr_reader :reference
    # GenomicInterval
    attr_accessor :interval
    attr_accessor :species
    attr_accessor :species_map

    def initialize
      @species_map = {}
    end

    # Set the reference sequence.
    #
    # @param source [FASTARangeReader, String, Pathname]
    def reference=(source)
      ref = case
            when source.is_a?(FASTARangeReader)
              source
            when source.respond_to?(:seek)
              # open file
              FASTARangeReader.new(source)
            when source.respond_to?(:start_with?) && source.start_with?('>')
              # FASTA string
              FASTARangeReader.new(StringIO.new(source))
            else
              FASTARangeReader.new(source.to_s)
            end
      @reference = ref
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

    def species_to_use
      species || species_map.keys
    end

    def tile
      parser.sequence_filter[:only_species] = species_to_use
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
      species_to_use.each { |s| text << '' }
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
          species_to_use.each_with_index do |species, i|
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
      species_to_use.zip(tile()) do |species, text|
        sp_out = species_map[species] || species
        f.puts ">#{sp_out}"
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
    attr_reader :f, :pos

    def initialize(fspec)
      if fspec.respond_to? :seek
        @f = fspec
      else
        reader_class = if fspec =~ /.gz$/
                         Zlib::GzipReader
                       else
                         File
                       end
        @f = reader_class.open(fspec)
      end
      position_at_start
    end

    GT = '>'.getbyte(0)

    def position_at_start
      first = f.readline
      raise "expected FASTA comment" unless first =~ /^>/
      @pos = 0
    end

    def read_interval(z_start, z_end)
      if z_start < pos
        position_at_start
      end
      data = ''
      region_size = z_end - z_start
      in_region = false
      f.each_line do |line_raw|
        if line_raw.getbyte(0) == GT
          raise "unexpected description line: #{line_raw.inspect}"
        end
        line = line_raw.strip
        end_pos = pos + line.size
        if (! in_region) && pos <= z_start && z_start < end_pos
          offset = z_start - pos
          end_offset = [(offset + region_size), line.size].min
          data << line.slice(offset...end_offset)
          in_region = true
        elsif in_region
          need = region_size - data.size
          raise "should not happen: region #{region_size}, data #{data.size}, need #{need}" if need < 0
          if need > line.size
            data << line
          else
            # last line
            data << line.slice(0, need)
            break
          end
        end
        @pos = end_pos
      end
      return data
    end
  end
end
