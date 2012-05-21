module Bio
  module MAF

    class ParseError < Exception; end

    class Header
      attr_accessor :vars, :alignment_params

      def initialize(vars, params)
        @vars = vars
        @alignment_params = params
      end

      def version
        vars[:version]
      end

      def scoring
        vars[:scoring]
      end

    end

    class Block
      attr_reader :vars, :sequences

      def initialize(*args)
        @vars, @sequences = args
      end

      def size
        sequences.size
      end

      def raw_seq(i)
        sequences[i]
      end

      def each_raw_seq
        sequences.each { |s| yield s }
      end
    end

    class Sequence
      attr_reader :source, :start, :size, :strand, :src_size, :text

      def initialize(*args)
        @source, @start, @size, @strand, @src_size, @text = args
      end

      def write_fasta(writer)
        writer.write("#{source}:#{start}-#{start + size}",
                     text)
      end
    end

    class Parser

      attr_accessor :header, :file_spec, :f, :r

      def initialize(file_spec)
        @file_spec = file_spec
        @f = File.open(file_spec)
        @r = LineReader.new(@f)
        _parse_header
      end

      def parse_maf_vars(line)
        vars = {}
        parts = line.strip.split
        parts[1..parts.size].each do |var_spec|
          if var_spec =~ /(\w+)=(.*)/
            vars[$1.to_sym] = $2
          else
            raise ParseError, "Malformed MAF variable: #{var_spec}"
          end
        end
        return vars
      end

      def _parse_header
        l1 = r.next_line
        raise ParseError, "Not a MAF file: #{file_spec}" unless l1 =~ /^##maf\s/
        vars = parse_maf_vars(l1)
        raise ParseError, "Missing MAF version!" unless vars.has_key? :version
        params = nil
        while true
          line = r.next_line.strip
          case line
          when /^#\s*(\S.*)/
            if params
              params << ' '
            else
              params = ''
            end
            params << $1
          when /^#/, /^$/
            break
          end
          if f =~ /#\s*.*/
            r.rewind
            break
          end
        end
        @header = Header.new(vars, params)
      end

      def skip_empty_lines
        while true do
          l = r.next_line
          if l =~ /^[\w]/
            r.rewind
            break
          end
        end
      end

      def parse_block
        skip_empty_lines
        a_line = r.next_line
        raise ParseError, "expected 'a'' line" unless a_line =~ /^a/
        vars = parse_maf_vars(a_line)
        seqs = []
        while true do
          begin
            line = r.next_line
          rescue EOFError
            break
          end
          case line
          when /^s/
            seqs << parse_seq(line)
          when /^[ieq]/
            # unhandled for now: synteny, empty regions, quality
          when /^a/
            # start of next alignment block
            r.rewind
            break
          end
        end
        return Block.new(vars, seqs)
      end

      def parse_seq(line)
        _, src, start, size, strand, src_size, text = line.split
        strand_sym =
          case strand
          when '+'
            :+
          when '-'
            :-
          else
            raise ParseError, "invalid strand #{strand}"
          end
        return Sequence.new(src,
                            start.to_i,
                            size.to_i,
                            strand_sym,
                            src_size.to_i,
                            text)
      end

      def each_block
        until f.eof?
          yield parse_block()
        end
      end

    end

    class LineReader
      attr_reader :f

      def initialize(f)
        @f = f
      end

      def next_line
        if @again
          @again = false
          return @last
        else
          return @last = @f.readline
        end
      end

      def rewind
        @again = true
      end
    end
    
  end
  
end
