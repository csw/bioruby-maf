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

    # The species of interest to extract from the MAF file. Will be
    # set as a {Parser#sequence_filter} for parsing. Defaults to the
    # keys of {#species_map}.
    #
    # @return [Array<String>]
    attr_accessor :species

    # A hash mapping species to their desired output names.
    #
    # @return [Hash]
    attr_accessor :species_map

    # The character used to fill regions where no sequence data is available for a particular species. Defaults to `*`.
    # @return [String]
    attr_reader   :fill_char

    attr_accessor :remove_absent_species

    def initialize
      @species_map = {}
      self.fill_char = '*'
      self.remove_absent_species = true
    end

    # Set the character to be used for filling regions with no
    # sequence data from the MAF file or a reference sequence.
    # @param c [String] a one-character String to fill with
    def fill_char=(c)
      unless c.is_a?(String) && c.length == 1
        raise ArgumentError, "not a single character: #{c.inspect}" 
      end
      @fill_char = c
    end

    # Set the reference sequence. This can be a {Pathname} or a
    # {String} giving the path to an optionally-gzipped FASTA file, an
    # open {IO} stream to a FASTA file, a String containing FASTA
    # data, or a {FASTARangeReader} instance.
    #
    # @param source [FASTARangeReader, String, Pathname, #readline]
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

    def species_for_output
      species_to_use.collect { |s| species_map[s] || s }
    end

    # Return an array of tiled sequence data, in the order given by
    # {#species_to_use}.
    # @return [Array<String>]
    def tile
      parser.sequence_filter[:only_species] = species_to_use
      parser.opts[:remove_gaps] = true
      LOG.debug { "finding blocks covering interval #{interval}." }
      blocks = index.find([interval], parser).sort_by { |b| b.vars[:score] }
      mask = Array.new(interval.length, :ref)
      i_start = interval.zero_start
      i_end = interval.zero_end
      if reference
        LOG.debug { "using a #{reference.class} reference." }
        ref_region = ref_data(i_start...i_end)
      end
      LOG.debug "tiling #{blocks.count} blocks."
      blocks.each do |block|
        ref = block.ref_seq
        LOG.debug { "tiling with block #{ref.start}-#{ref.end}" }
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
          fill_text = fill_char * range_size
          nonref_text.each { |t| t << fill_text }
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
              # no alignment for this one here, use the fill char
              sp_text << fill_char * (t_range.end - t_range.begin)
            end
          end
        end
      end
      if remove_absent_species
        non_fill = non_fill_re
        LOG.debug { "searching for non-fill characters with #{non_fill}" }
        text.each_with_index do |seq, i|
          unless non_fill.match(seq)
            text[i] = nil
          end
        end
      end
      text
    end

    def non_fill_re
      fill_esc = Regexp.escape(fill_char)
      Regexp.compile("[^#{fill_esc}]")
    end

    def output_text
      species_for_output.zip(tile()).reject { |s, t| t.nil? }
    end

    # Tile sequences to build a new {Bio::BioAlignment::Alignment
    # Alignment} object. This will have one
    # {Bio::BioAlignment::Sequence Sequence} per entry in {#species}
    # or {#species_map}, in the same order. Each sequence will have an
    # {Bio::BioAlignment::Sequence#id id} given by {#species_map} or,
    # if none is present, the identifier from {#species}.
    #
    # @return [Bio::BioAlignment::Alignment]
    # @api public
    def build_bio_alignment
      out = output_text.to_a
      Bio::BioAlignment::Alignment.new(out.collect { |e| e[1] },
                                       out.collect { |e| e[0] })
    end

    # Write a FASTA representation of the tiled sequences to the given
    # output stream.
    #
    # @param [#puts] f the output stream to write the FASTA data to.
    # @api public
    def write_fasta(f)
      output_text.each do |sp_out, text|
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
