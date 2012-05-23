require 'strscan'

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
        sequences.fetch(i)
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
          when /^a/
            r.rewind
            break
          end
          #if f =~ /#\s*.*/
          #  r.rewind
          #  break
          #end
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
        raise ParseError, "expected 'a' line" unless a_line =~ /^a/
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

    class ChunkParser

      ## Parses alignment blocks by reading a chunk of the file at a time.

     attr_reader :header, :file_spec, :f, :s, :at_end, :last_block_pos

      CHUNK_SIZE = 8 * 1024 * 1024

      def initialize(file_spec)
        @file_spec = file_spec
        @f = File.open(file_spec)
        @s = StringScanner.new(read_chunk())
        set_last_block_pos!
        @at_end = false
        _parse_header()
      end

      BLOCK_START = /^(?=a)/
      BLOCK_START_OR_EOS = /(?:^(?=a))|\z/
      EOL_OR_EOF = /\n|\z/

      def read_chunk
        f.read(CHUNK_SIZE)
      end

      def set_last_block_pos!
        @last_block_pos = s.string.rindex(BLOCK_START)
      end

      def _parse_header
        parse_error("not a MAF file") unless s.scan(/##maf\s*/)
        vars = parse_maf_vars()
        align_params = nil
        while s.scan(/^#\s*(.+?)\n/)
          if align_params == nil
            align_params = s[1]
          else
            align_params << ' ' << s[1]
          end
        end
        @header = Header.new(vars, align_params)
        s.skip_until BLOCK_START || parse_error("Cannot find block start!")
      end

      ## Invariants:
      ##
      ## cur_chunk: current chunk
      ## s: StringScanner, positioned at start of the block to parse
      ##   (start_pos = s.pos)
      ##
      ## On finding the start of a block:
      ## Look for the start of the next block.
      ##   If not found, see whether we are at EOF.
      ##     If at EOF: last block.
      ##     If not:
      ##       Read the next chunk
      ##       Find the start of the next block in that chunk
      ##       Concatenate the two block fragments
      ##       Parse the resulting block
      ##       Promote the next scanner, positioned

      def peek_for_next_block
        start_pos = s.pos
        s.pos += 1
        next_block_offset = s.search_full(BLOCK_START, false, false)
        next_block_pos = s.pos + next_block_offset if next_block_offset
        s.pos = start_pos
        return next_block_pos
      end

      def parse_block
        return nil if at_end
        block = nil
        #next_block_pos = peek_for_next_block()
        s.skip_until BLOCK_START
        if s.pos != last_block_pos
          # in non-trailing block
          block = parse_block_data
          #s.pos = next_block_pos
        else
          # in trailing block fragment
          if f.eof?
            # last block, parse it as is
            block = parse_block_data
            @at_end = true
          else
            next_chunk = read_chunk
            next_scanner = StringScanner.new(next_chunk)
            if s.string[s.string.size - 1] == "\n"
              pat = BLOCK_START_OR_EOS
            else
              pat = /(?:\n(?=a))|\z/
            end
            leading_frag = next_scanner.scan_until(pat)
            unless leading_frag
              parse_error("no leading fragment match!")
            end
            joined_block = s.rest + leading_frag
            @s = StringScanner.new(joined_block)
            block = parse_block_data
            @s = next_scanner
            if s.eos?
              @at_end = true
            else
              set_last_block_pos!
            end
          end
        end
        return block
      end

      def rest_of_line
        s.scan_until(EOL_OR_EOF) || parse_error("Cannot scan to newline")
      end

      def parse_error(msg)
        s_start = [s.pos - 10, 0].max
        s_end = [s.pos + 10, s.string.length].min
        if s_start > 0
          left = s.string[s_start..(s.pos - 1)]
        else
          left = ''
        end
        right = s.string[s.pos..s_end]
        extra = "pos #{s.pos}, last #{last_block_pos}"

        raise ParseError, "#{msg} at: '#{left}>><<#{right}' (#{extra})"
      end

      def parse_block_data
        s.scan(/^a\s*/) || parse_error("bad a line")
        block_vars = parse_maf_vars()
        seqs = []
        while s.scan(/^([sieqa])\s+/)
          case s[1]
          when 's'
            seqs << parse_seq
          when 'i', 'e', 'q'
            # ignore
            s.skip_until(EOL_OR_EOF)
          when 'a'
            parse_error "unexpectedly reached next block"
          end
        end
        return Block.new(block_vars, seqs)
      end

      def parse_maf_vars
        vars = {}
        while s.scan(/(\w+)=(\S*)\s+/) do
          vars[s[1].to_sym] = s[2]
        end
        vars
      end

      STRAND_SYM = {
        '+' => :+,
        '-' => :-
      }

      def parse_seq
        src, start, size, strand, src_size, text = rest_of_line.split
        return Sequence.new(src,
                            start.to_i,
                            size.to_i,
                            STRAND_SYM.fetch(strand),
                            src_size.to_i,
                            text)
      end

      def each_block
        until at_end
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
