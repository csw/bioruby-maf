require 'bio-alignment'

module Bio
  class GenomicInterval
    def intersection(other)
      raise ArgumentError unless self.chrom == other.chrom
      GenomicInterval.new(self.chrom,
                          [self.chr_start, other.chr_start].max,
                          [self.chr_end, other.chr_end].min)
    end
  end

  module MAF

    # A MAF header, containing the variable-value pairs from the first
    # line of the file as well as the alignment parameters.
    # @api public
    class Header
      # Variable-value pairs from the ##maf line
      # @return [Hash]
      attr_accessor :vars
      # Alignment parameters from the MAF header.
      # @return [Hash]
      attr_accessor :alignment_params

      # Create a default header with version=1.
      # @return [Header]
      def Header.default
        Header.new({:version => 1}, nil)
      end

      def initialize(vars, params)
        @vars = vars
        @alignment_params = params
      end

      # The required version parameter.
      # @return [String]
      def version
        vars[:version]
      end

      # The optional scoring parameter, if present.
      # @return [String]
      def scoring
        vars[:scoring]
      end

    end

    # A MAF alignment block.
    # @api public
    class Block
      # Parameters from the 'a' line starting the alignment block.
      attr_reader :vars
      # Sequences, one per 's' or 'e' line.
      # @return [Array<Sequence>]
      attr_reader :sequences
      # Offset of the alignment block within the MAF file, in bytes.
      # @return [Integer]
      attr_reader :offset
      # Size of the alignment block within the MAF file, in bytes.
      # @return [Integer]
      attr_reader :size

      def initialize(vars, sequences, offset, size, filtered)
        #raise ArgumentError, "no sequences given for block at offset #{offset}!" unless sequences && sequences.first
        @vars = vars
        @sequences = sequences
        @offset = offset
        @size = size
        @filtered = filtered
      end

      def ref_seq
        sequences[0]
      end

      def raw_seq(i)
        sequences.fetch(i)
      end

      def each_raw_seq
        sequences.each { |s| yield s }
      end

      # Text size of the alignment block. This is the number of text
      # characters in each line of sequence data, including dashes and
      # other gaps in the sequence.
      def text_size
        sequences.first.text.size
      end

      # Whether this block has been modified by a parser filter.
      # @return [Boolean]
      def filtered?
        @filtered
      end

      def to_bio_alignment
        ba_seq = sequences.collect { |s| s.to_bio_alignment }
        Bio::BioAlignment::Alignment.new(ba_seq)
      end

      GAP = /-+/

      # Remove gaps present in all sequences. These would generally
      # occur when some sequences have been filtered out.
      # @see #remove_gaps!
      # @see Parser#sequence_filter
      def find_gaps
        ref_s = StringScanner.new(sequences.first.text)
        others = sequences.slice(1, sequences.size - 1).reject { |s| s.empty? }.collect { |s| StringScanner.new(s.text) }
        gaps = []
        while ref_s.scan_until(GAP)
          offset = ref_s.pos - ref_s.matched_size
          others.each { |s| s.pos = offset }
          unless others.find { |s| ! s.scan(GAP) }
            # all matched
            gap_size = [ref_s.matched_size,
                        others.map {|s| s.matched_size}.min].min
            gaps << [offset, gap_size]
          end
        end
        gaps
      end

      # Remove gaps present in all sequences. These would generally
      # occur when some sequences have been filtered out.
      # @see #find_gaps
      # @see Parser#sequence_filter
      def remove_gaps!
        gaps = find_gaps()
        gaps.reverse_each do |offset, len|
          sequences.each do |seq|
            seq.delete_text(offset, len)
          end
        end
        gaps.size
      end

      # Returns a new Block covering only the region where it overlaps
      # the given interval.
      # @param [Bio::GenomicInterval] interval to slice the block with
      # @return [Block] block covering intersection with interval
      def slice(interval)
        case interval.compare(ref_seq.interval)
        when :equal
          return self
        when :contains, :contained_by, :left_overlapped, :right_overlapped
          _slice(interval.intersection(ref_seq.interval))
        when :left_adjacent, :right_adjacent, :left_off, :right_off
          raise "Cannot slice a block with a non-overlapping interval! Block #{ref_seq.interval}, interval #{interval}"
        when :different_chrom
          raise "Cannot slice a block with reference sequence #{ref_seq.source} using an interval on #{interval.chrom}!"
        else
          raise "Unhandled comparison result: #{interval.compare(ref_seq.interval)}"
        end
      end

      def _slice(interval)
        #offset = interval.zero_start - ref_seq.start
        #i_len  = interval.length
        range = _slice_text_range(interval)
        s2 = sequences.collect { |s| s.slice(range) }
        v2 = vars.dup
        #v2[:score] = '0.0'
        # TODO: should the filtered param be #modified? instead?
        Block.new(v2, s2, offset, size, @filtered)
      end

      def _slice_text_range(interval)
        i_start  = interval.zero_start
        i_end    = interval.zero_end
        g_pos    = ref_seq.start
        t_start  = nil
        t_end    = nil
        ref_seq.text.each_char.each_with_index do |c, t_pos|
          if c != '-'
            # non-gap
            if g_pos == i_start
              t_start = t_pos
            end
            g_pos += 1
            if t_start && g_pos == i_end
              t_end = t_pos + 1
              break
            end
          end
        end
        unless t_start && t_end
          raise "did not find start and end for #{interval} in #{ref_seq.inspect}!"
        end
        return t_start...t_end
      end

      def joinable_with?(other)
        if sequences.size == other.sequences.size
          r1 = ref_seq
          r2 = other.ref_seq
          return false if r1.source != r2.source
          return false if r1.end != r2.start
          rest = sequences.each_with_index
          rest.next
          mismatch = rest.find do |s1, i|
            s2 = other.seq_from(s1.source, i)
            (! s2) || ! s1.joinable_with?(s2)
          end
          return (! mismatch)
        else
          return false
        end
      end

      def join(other)
        nseq = sequences.each_with_index.collect do |s1, i|
          s2 = other.seq_from(s1.source, i)
          s1.join(s2)
        end
        v2 = vars.dup
        v2[:score] = '0.0'
        Block.new(v2, nseq, offset, nil, @filtered)
      end

      def seq_from(src, pos_guess)
        sg = sequences[pos_guess]
        if sg.source == src
          sg
        else
          sequences.find { |s| s.source == src }
        end
      end

    end

    # A sequence within an alignment block.
    # @api public
    class Sequence
      # @return [String] Source sequence name.
      attr_reader :source
      # @return [Integer] Zero-based start position.
      attr_reader :start
      # @return [Integer] Size of aligning region in source sequence.
      attr_reader :size
      # :+ or :-, indicating which strand the alignment is to.
      # @return [Symbol]
      attr_reader :strand
      # Size of the entire source sequence, not just the aligning
      # region.
      # @return [Integer]
      attr_reader :src_size
      # Sequence data for the alignment, including insertions.
      # @return [String]
      attr_reader :text
      # Array of raw synteny information from 'i' line.
      # @return [Array<String>]
      attr_accessor :i_data
      # Quality string from 'q' line.
      # @return [String]
      attr_accessor :quality
      alias_method :source_size, :src_size

      def initialize(*args)
        @source, @start, @size, @strand, @src_size, @text = args
      end

      def end
        start + size
      end

      def interval
        GenomicInterval.zero_based(self.source, self.start, self.end)
      end

      def slice(range)
        before = text.slice(0...(range.begin))
        non_gap_prev = before.delete("-").size
        new_text = text.slice(range)
        unless new_text
          raise "could not extract slice #{range} from #{self.inspect}!"
        end
        non_gap_text = new_text.delete("-").size
        s2 = Sequence.new(source,
                          start + non_gap_prev,
                          non_gap_text,
                          strand,
                          src_size,
                          new_text)
        s2.quality = quality.slice(range) if quality
        # TODO: what to do with synteny data?
        s2
      end

      # Whether this sequence is empty. Only true for {EmptySequence}
      # instances from 'e' lines.
      def empty?
        false
      end

      def gapped?
        size != text.size
      end

      I_STATUS = {
        'C' => :contiguous,
        'I' => :intervening,
        'N' => :first,
        'n' => :first_bridged,
        'M' => :missing_data,
        'T' => :tandem
      }

      def decode_status_char(c)
        I_STATUS[c] || raise("Unsupported status character #{c}!")
      end

      def left_status_char
        i_data && i_data[0]
      end

      def left_status
        i_data && decode_status_char(left_status_char())
      end

      def left_count
        i_data && i_data[1].to_i
      end

      def right_status_char
        i_data && i_data[2]
      end

      def right_status
        i_data && decode_status_char(right_status_char())
      end

      def right_count
        i_data && i_data[3].to_i
      end

      def species
        parts = source.split('.', 2)
        parts.size == 2 ? parts[0] : nil
      end

      def delete_text(offset, len)
        unless empty?
          text.slice!(offset, len)
          if quality
            quality.slice!(offset, len)
          end
        end
      end

      def to_bio_alignment
        Bio::BioAlignment::Sequence.new(source, text)
      end

      def write_fasta(writer)
        writer.write("#{source}:#{start}-#{start + size}",
                     text)
      end

      def joinable_with?(o)
        (self.end == o.start) \
        && (self.strand == o.strand) \
        && (self.empty? == o.empty?)
      end

      def join(o)
        s2 = Sequence.new(source,
                          start,
                          size + o.size,
                          strand,
                          src_size,
                          text + o.text)
        if quality && o.quality
          s2.quality = quality + o.quality
        end
        s2
      end

      # Maps the given zero-based genomic range onto a range of string
      # offsets, suitable for extracting the text for the given range
      # from #text.
      #
      # @see String#slice
      def text_range(range)
        r_end = range.exclude_end? ? range.end : range.end + 1
        r_size = r_end - range.begin
        if range.begin == start && r_size == size
          # special case, entire text
          0...text.size
        else
          if range.begin < start || r_end > self.end
            raise "Range #{range} outside sequence bounds; start #{start}, size #{size}"
          end
          if ! gapped?
            # no gaps, can map indexes directly
            (range.begin - start)...(r_end - start)
          else
            # gaps present
            g_start = start     # genomic position of the start
            t_start = 0         # text position of the start
            m_begin = nil       # beginning of match
            match = nil
            text.scan(/(\w+|-+)/) do |parts|
              part = parts[0]
              if part[0] != '-'
                # sequence text
                g_end = g_start + part.size
                if g_start <= range.begin && range.begin < g_end
                  offset_in_part = range.begin - g_start
                  m_begin = offset_in_part + t_start
                end
                if g_start <= r_end && r_end <= g_end
                  raise "reached end before start!" unless m_begin
                  offset_in_part = r_end - g_start
                  m_end = offset_in_part + t_start
                  match = m_begin...m_end
                  break
                end
                g_start = g_end
              else
                # gap
              end
              t_start += part.size
            end
            raise "no match found!" unless match
            return match
          end
        end
      end
    end

    # An empty sequence record from an 'e' line.
    #
    # This indicates that "there isn't aligning DNA for a species but
    # that the current block is bridged by a chain that connects
    # blocks before and after this block" (MAF spec).
    # @api public
    class EmptySequence < Sequence
      attr_reader :status

      def initialize(*args)
        super(*args[0..4])
        @status = args[5]
      end

      def text
        ''
      end

      def slice(offset, len)
        self
      end

      def join(o)
        EmptySequence.new(source,
                          start,
                          size + o.size,
                          strand,
                          src_size,
                          @status)
      end

      def empty?
        true
      end

      def write_fasta(writer)
        raise "empty sequence output not implemented!"
      end
    end

  end
  
end
