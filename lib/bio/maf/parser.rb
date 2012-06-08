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

      def text_size
        sequences.first.text.size
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

    class ChunkReader
      attr_accessor :chunk_size, :chunk_shift, :pos
      attr_reader :f
      def initialize(f, chunk_size)
        @f = f
        self.chunk_size = chunk_size
        @pos = 0
      end

      def chunk_size=(size)
        check_chunk_size(size)
        @chunk_size = size
        # power of 2 so don't worry about rounding
        @chunk_shift = Math.log2(size).to_i
      end

      def check_chunk_size(size)
        if size < 1
          raise "Invalid chunk size: #{size}"
        end
        ## test whether it is a power of 2
        ## cf. http://bit.ly/JExNc4
        if size & (size - 1) != 0
          raise "Invalid chunk size (not a power of 2): #{size}}"
        end
      end

      # Reads the next chunk of the file.
      def read_chunk
        chunk = f.read(@chunk_size)
        @pos += chunk.bytesize
        return chunk
      end

      # Reads a chunk of the file, starting at the specified
      # offset. If the offset is not the start of a chunk, the read
      # will be shortened so that it ends on a chunk boundary.
      def read_chunk_at(offset)
        f.seek(offset)
        next_chunk_start = ((offset >> chunk_shift) + 1) << chunk_shift
        chunk = f.read(next_chunk_start - offset)
        @pos = offset + chunk.bytesize
        return chunk
      end
    end

    class Parser

      ## Parses alignment blocks by reading a chunk of the file at a time.

      attr_reader :header, :file_spec, :f, :s, :cr, :at_end
      attr_reader :chunk_start, :chunk_size, :last_block_pos
      attr_accessor :sequence_filter

      SEQ_CHUNK_SIZE = 8 * 1024 * 1024

      def initialize(file_spec, opts={})
        @chunk_size = opts[:chunk_size] || SEQ_CHUNK_SIZE
        @chunk_start = 0
        @file_spec = file_spec
        @f = File.open(file_spec)
        @cr = ChunkReader.new(@f, chunk_size)
        @s = StringScanner.new(read_chunk())
        set_last_block_pos!
        @at_end = false
        _parse_header()
      end

      BLOCK_START = /^(?=a)/
      BLOCK_START_OR_EOS = /(?:^(?=a))|\z/
      EOL_OR_EOF = /\n|\z/

      def read_chunk
        cr.read_chunk
      end

      def fetch_blocks(fetch_list, filters=nil)
        ## fetch_list: array of [offset, length, block_count] tuples
        ## returns array of Blocks
        return fetch_blocks_merged(merge_fetch_list(fetch_list))
      end

      def fetch_blocks_merged(fetch_list)
        r = []
        old_chunk_size = @chunk_size
        @chunk_size = 4096
        @at_end = false
        begin
          fetch_list.each do |offset, len, block_offsets|
            if (chunk_start <= offset) \
              && (offset < (chunk_start + s.string.size))
              ## the selected offset is in the current chunk
              s.pos = offset - chunk_start
            else
              chunk = cr.read_chunk_at(offset)
              @chunk_start = offset
              @s = StringScanner.new(chunk)
            end
            block_offsets.each do |expected_offset|
              block = parse_block
              unless block
                raise "Parse error: expected a block at offset #{expected_offset} but could not parse one!"
              end
              unless block.offset == expected_offset
                raise "Parse error: got block with offset #{block.offset}, expected #{expected_offset}!"
              end
              r << block
            end
          end
          return r
        ensure
          @chunk_size = old_chunk_size
        end
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
            r.last[2] << cur[0]
          else
            cur << [cur[0]]
            r << cur
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
            next_chunk_start = cr.pos
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
            if sequence_filter
              if sequence_filter[:only_species]
                src_sp = src.split('.', 2)[0]
                m = sequence_filter[:only_species].find { |sp| src_sp == sp }
                next unless m
              end
            end
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
