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
      attr_reader :vars, :sequences, :offset, :size

      def initialize(*args)
        @vars, @sequences, @offset, @size = args
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
      alias_method :source_size, :src_size

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
      attr_reader :chunk_start

      CHUNK_SIZE = 8 * 1024 * 1024

      def initialize(file_spec)
        @file_spec = file_spec
        @f = File.open(file_spec)
        @s = StringScanner.new(read_chunk())
        @chunk_start = 0
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

      def merge_fetch_list(orig_fl)
        fl = orig_fl.dup
        r = []
        until fl.empty? do
          cur = fl.shift
          if r.last && (r.last[0] + r.last[1]) == cur[0]
            # contiguous with the previous one
            # add to length and increment count
            r.last[1] += cur[1]
            r.last[2] += 1
          else
            r << cur
            r.last << 1
          end
        end
        return r
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
            # Trailing block fragment, but there is another chunk.
            # Read the next chunk, find the start of the first
            # alignment block, and concatenate this trailing block
            # fragment with the leading fragment before the start of
            # that next block. Parse the resulting joined block, then
            # position the scanner to parse the next block.
            next_chunk_start = chunk_start + s.string.bytesize
            next_chunk = read_chunk
            # Find the next alignment block
            next_scanner = StringScanner.new(next_chunk)
            # If this trailing fragment ends with a newline, then an
            # 'a' at the beginning of the leading fragment is the
            # start of the next alignment block.
            if s.string[s.string.size - 1] == "\n"
              pat = BLOCK_START_OR_EOS
            else
              pat = /(?:\n(?=a))|\z/
            end
            leading_frag = next_scanner.scan_until(pat)
            unless leading_frag
              parse_error("no leading fragment match!")
            end
            # Join the fragments and parse them
            joined_block = s.rest + leading_frag
            @chunk_start = chunk_start + s.pos
            @s = StringScanner.new(joined_block)
            block = parse_block_data
            # Set up to parse the next block
            @s = next_scanner
            @chunk_start = next_chunk_start
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
        block_start_pos = s.pos
        block_offset = chunk_start + block_start_pos
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
        return Block.new(block_vars,
                         seqs,
                         block_offset,
                         s.pos - block_start_pos)
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

  end
  
end
