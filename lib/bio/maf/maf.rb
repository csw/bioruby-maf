module Bio
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
        @vars = vars
        @sequences = sequences
        @offset = offset
        @size = size
        @filtered = filtered
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

      # Whether this sequence is empty. Only true for {EmptySequence}
      # instances from 'e' lines.
      def empty?
        false
      end

      def delete_text(offset, len)
        unless empty?
          text.slice!(offset, len)
          if quality
            quality.slice!(offset, len)
          end
        end
      end

      def write_fasta(writer)
        writer.write("#{source}:#{start}-#{start + size}",
                     text)
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

      def empty?
        true
      end

      def write_fasta(writer)
        raise "empty sequence output not implemented!"
      end
    end

  end
  
end
