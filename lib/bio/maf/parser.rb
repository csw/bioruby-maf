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
      attr_reader :score, :sequences
    end

    class Sequence
      attr_accessor 
    end

    class Parser

      attr_accessor :header, :file_spec, :f, :r

      def initialize(file_spec)
        @file_spec = file_spec
        @f = File.open(file_spec)
        @r = LineReader.new(@f)
        _parse_header
      end

      def _parse_header
        l1 = r.next_line
        raise ParseError, "Not a MAF file: #{file_spec}" unless l1 =~ /^##maf\s/
        vars = {}
        parts = l1.strip.split
        parts[1..parts.size].each do |var_spec|
          if var_spec =~ /(\w+)=(.*)/
            vars[$1.to_sym] = $2
          else
            raise ParseError, "Malformed MAF variable: #{var_spec}"
          end
        end
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
