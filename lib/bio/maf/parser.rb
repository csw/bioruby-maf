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

      ## On finding the start of a block:
      ## See whether we are at the last block in the chunk.
      ##   If at the last block:
      ##     If at EOF: last block.
      ##     If not:
      ##       Read the next chunk
      ##       Find the start of the next block in that chunk
      ##       Concatenate the two block fragments
      ##       Parse the resulting block
      ##       Promote the next scanner, positioned

      def parse_block
        return nil if at_end
        block = nil
        s.skip_until BLOCK_START
        if s.pos != last_block_pos
          # in non-trailing block
          block = parse_block_data
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
        payload = s.scan_until(/^(?=a)/)
        unless payload
          payload = s.rest
          s.pos = s.string.size # jump to EOS
        end
        payload.split("\n").each do |line|
          case line[0]
          when 's'
            _, src, start, size, strand, src_size, text = line.split
            seqs << Sequence.new(src,
                                 start.to_i,
                                 size.to_i,
                                 STRAND_SYM.fetch(strand),
                                 src_size.to_i,
                                 text)
          when 'i', 'e', 'q', '#', nil
            next
          else
            parse_error "unexpected line: '#{line}'"
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
