require 'strscan'
require 'java' if RUBY_PLATFORM == 'java'

# @api public
module Bio
  # @api public
  module MAF

    # @api public
    class ParseError < Exception; end

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

      def initialize(*args)
        @vars, @sequences, @offset, @size = args
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

      # Whether this sequence is empty. Only true for {EmptySequence}
      # instances from 'e' lines.
      def empty?
        false
      end

      def gapped?
        size != text.size
      end

      def species
        parts = source.split('.', 2)
        parts.size == 2 ? parts[0] : nil
      end

      def write_fasta(writer)
        writer.write("#{source}:#{start}-#{start + size}",
                     text)
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

      def empty?
        true
      end

      def write_fasta(writer)
        raise "empty sequence output not implemented!"
      end
    end

    # Reads MAF files in chunks.
    # @api private
    class ChunkReader
      # Size, in bytes, of the chunks to read. Must be a power of 2.
      # @return [Integer]
      attr_accessor :chunk_size
      # Current position in the file.
      # @return [Integer]
      attr_accessor :pos
      # {File} from which chunks are read.
      # @return [File]
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
        # @chunk_shift = Math.log2(size).to_i
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
      # @return [String] Next {#chunk_size} bytes of MAF data.
      def read_chunk
        chunk = f.read(@chunk_size)
        @pos += chunk.bytesize if chunk
        return chunk
      end

      # Reads a chunk of the file.
      #
      # Currently always reads size_hint bytes but this may change
      # with BGZF support.
      #
      # @param [Integer] offset file offset to read from.
      # @param [Integer] size_hint desired size of chunk.
      # @return [String] Chunk of MAF data.
      def read_chunk_at(offset, size_hint=@chunk_size)
        f.seek(offset)
        chunk = f.read(size_hint)
        @pos = offset + chunk.bytesize
        return chunk
      end
    end

    # Variant ChunkReader using a read-ahead thread with internal
    # queue for sequential parsing. Not useful for random-access
    # parsing.
    #
    # Only beneficial on JRuby.
    class ThreadedChunkReader < ChunkReader

      def initialize(f, chunk_size, buffer_size=64)
        super(f, chunk_size)
        @buffer = SizedQueue.new(buffer_size)
        @eof_reached = false
        start_read_ahead
      end

      # Spawn a read-ahead thread. Called from {#initialize}.
      def start_read_ahead
        @read_thread = Thread.new { read_ahead }
      end

      # Read ahead into queue.
      def read_ahead
        # n = 0
        begin
          f_pos = 0
          until f.eof?
            chunk = f.read(@chunk_size)
            @buffer << [f_pos, chunk]
            f_pos += chunk.bytesize
            # n += 1
            # if (n % 100) == 0
            #   $stderr.puts "buffer size: #{@buffer.size}"
            # end
          end
          @eof_reached = true
        rescue Exception
          @read_ahead_ex = $!
          $stderr.puts "read_ahead aborting: #{$!}"
        end
      end

      # (see ChunkReader#read_chunk)
      def read_chunk
        raise "readahead failed: #{@read_ahead_ex}" if @read_ahead_ex
        if @eof_reached && @buffer.empty?
          return nil
        else
          c_pos, chunk = @buffer.shift()
          @pos = c_pos
          return chunk
        end
      end

    end

    # MAF parsing code useful for sequential and random-access parsing.
    module MAFParsing
      
      BLOCK_START = /^(?=a)/
      BLOCK_START_OR_EOS = /(?:^(?=a))|\z/
      EOL_OR_EOF = /\n|\z/

      def set_last_block_pos!
        @last_block_pos = s.string.rindex(BLOCK_START)
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

      # Parse the block at the current position, joining fragments
      # across chunk boundaries if necessary.
      #
      # @return [Block] alignment block
      # @api public
      def parse_block
        return nil if at_end
        if s.pos != last_block_pos
          # in non-trailing block
          parse_block_data
        else
          # in trailing block fragment
          parse_trailing_fragment
        end
      end

      ## Read chunks and accumulate a leading fragment until we
      ## encounter a block start or EOF.
      def gather_leading_fragment
        leading_frag = ''
        while true
          next_chunk_start = cr.pos
          next_chunk = cr.read_chunk
          if next_chunk
            next_scanner = StringScanner.new(next_chunk)
            # If this trailing fragment ends with a newline, then an
            # 'a' at the beginning of the leading fragment is the
            # start of the next alignment block.
            if trailing_nl?(leading_frag) || trailing_nl?(s.string)
              pat = BLOCK_START
            else
              pat = /(?:\n(?=a))/
            end
            frag = next_scanner.scan_until(pat)
            if frag
              # got block start
              leading_frag << frag
              break
            else
              # no block start in this
              leading_frag << next_chunk
            end
          else
            # EOF
            @at_end = true
            break
          end
        end
        return leading_frag, next_scanner, next_chunk_start
      end

      # Join the trailing fragment of the current chunk with the
      # leading fragment of the next chunk and parse the resulting
      # block.
      #
      # @return [Block] the alignment block.

      def parse_trailing_fragment
        leading_frag, next_scanner, next_chunk_start = gather_leading_fragment
        # join fragments and parse
        trailing_frag = s.rest
        joined_block = trailing_frag + leading_frag
        @chunk_start = chunk_start + s.pos
        @s = StringScanner.new(joined_block)
        begin
          block = parse_block_data
        rescue ParseError => pe
          parse_error "Could not parse joined fragments: #{pe}\nTRAILING: #{trailing_frag}\nLEADING: #{leading_frag}"
        end
        # Set up to parse the next block
        @s = next_scanner
        @chunk_start = next_chunk_start
        unless @at_end
          set_last_block_pos!
        end
        return block
      end

      # Raise a {ParseError}, indicating position within the MAF file
      # and the chunk as well as the text surrounding the current
      # scanner position.
      #
      # @param [String] msg the error message
      def parse_error(msg)
        s_start = [s.pos - 10, 0].max
        s_end = [s.pos + 10, s.string.length].min
        if s_start > 0
          left = s.string[s_start..(s.pos - 1)]
        else
          left = ''
        end
        right = s.string[s.pos..s_end]
        extra = "pos #{s.pos} [#{chunk_start + s.pos}], last #{last_block_pos}"

        raise ParseError, "#{msg} at: '#{left}>><<#{right}' (#{extra})"
      end

      S = 's'.getbyte(0)
      I = 'i'.getbyte(0)
      E = 'e'.getbyte(0)
      Q = 'q'.getbyte(0)
      COMMENT = '#'.getbyte(0)

      # Parse a {Block} from the current position. Requires that {#s}
      # and {#chunk_start} be set correctly.
      #
      # @return [Block] the alignment block.
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
        lines = payload.split("\n")
        until lines.empty?
          line = lines.shift
          first = line.getbyte(0)
          if first == S
            seq = parse_seq_line(line, sequence_filter)
            seqs << seq if seq
          elsif first == E && parse_empty
            e_seq = parse_empty_line(line, sequence_filter)
            seqs << e_seq if e_seq
          elsif first == I && parse_extended
            parts = line.split
            parse_error("wrong i source #{parts[1]}!") unless seqs.last.source == parts[1]
            seqs.last.i_data = parts.slice(2..6)
          elsif first == Q && parse_extended
            _, src, quality = line.split
            parse_error("wrong q source #{src}!") unless seqs.last.source == src
            seqs.last.quality = quality
          elsif [I, E, Q, COMMENT, nil].include? first
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

      # Parse an 's' line.
      # @return [Sequence]
      def parse_seq_line(line, filter)
        _, src, start, size, strand, src_size, text = line.split
        return nil if filter && ! seq_filter_ok?(src, filter)
        begin
          Sequence.new(src,
                       start.to_i,
                       size.to_i,
                       STRAND_SYM.fetch(strand),
                       src_size.to_i,
                       text)
        rescue KeyError
          parse_error "invalid sequence line: #{line}"
        end
      end

      # Parse an 'e' line.
      # @return [EmptySequence]
      def parse_empty_line(line, filter)
        _, src, start, size, strand, src_size, status = line.split
        return nil if filter && ! seq_filter_ok?(src, filter)
        begin
          EmptySequence.new(src,
                            start.to_i,
                            size.to_i,
                            STRAND_SYM.fetch(strand),
                            src_size.to_i,
                            status)
        rescue KeyError
          parse_error "invalid empty sequence line: #{line}"
        end
      end

      # Indicates whether the given sequence source should be parsed,
      # given the current sequence filters.
      def seq_filter_ok?(src, filter)
        if filter[:only_species]
          src_sp = src.split('.', 2)[0]
          m = filter[:only_species].find { |sp| src_sp == sp }
          return m
        else
          return true
        end
      end

      # Parse key-value pairs from the MAF header or an 'a' line.
      # @return [Hash]
      def parse_maf_vars
        vars = {}
        while s.scan(/(\w+)=(\S*)\s+/) do
          vars[s[1].to_sym] = s[2]
        end
        vars
      end

      # Does `string` have a trailing newline?
      def trailing_nl?(string)
        if string.empty?
          false
        else
          s.string[s.string.size - 1] == "\n"
        end
      end

      STRAND_SYM = {
        '+' => :+,
        '-' => :-
      }
    end

    # A MAF parsing context, used for random-access parsing.
    class ParseContext
      include MAFParsing
      attr_accessor :f, :s, :cr, :parser
      attr_accessor :chunk_start, :last_block_pos, :at_end

      def initialize(fd, chunk_size, parser, opts)
        @f = fd
        @parser = parser
        reader = opts[:chunk_reader] || ChunkReader
        @cr = reader.new(@f, chunk_size)
        @last_block_pos = -1
      end

      def sequence_filter
        parser.sequence_filter
      end

      def parse_empty
        parser.parse_empty
      end

      def parse_extended
        parser.parse_extended
      end

      def set_last_block_pos!
        @last_block_pos = s.string.rindex(BLOCK_START)
      end

      # Fetch and parse blocks at given `offset` and `len`
      # @param [Integer] offset Offset to start parsing at.
      # @param [Integer] len Number of bytes to read.
      # @param [Array] block_offsets Offsets of blocks to parse.
      # @return [Array<Block>]
      def fetch_blocks(offset, len, block_offsets)
        start_chunk_read_if_needed(offset, len)
        # read chunks until we have the entire merged set of
        # blocks ready to parse
        # to avoid fragment joining
        append_chunks_to(len)
        # parse the blocks
        Enumerator.new do |y|
          block_offsets.each do |expected_offset|
            block = parse_block
            ctx.parse_error("expected a block at offset #{expected_offset} but could not parse one!") unless block
            ctx.parse_error("got block with offset #{block.offset}, expected #{expected_offset}!") unless block.offset == expected_offset
            y << block
          end
        end
      end

      def start_chunk_read_if_needed(offset, len)
        if chunk_start \
          && (chunk_start <= offset) \
          && (offset < (chunk_start + s.string.size))
          ## the selected offset is in the current chunk
          s.pos = offset - chunk_start
        else
          chunk = cr.read_chunk_at(offset, len)
          @chunk_start = offset
          @s = StringScanner.new(chunk)
        end
      end

      def append_chunks_to(len)
        # XXX: need to rethink this for BGZF; prefetching ChunkReader
        while s.string.size < len
          s.string << cr.read_chunk()
        end
      end

    end

    # MAF parser, used for sequential and random-access parsing.
    #
    # Options:
    #
    #  * `:parse_extended`: whether to parse 'i' and 'q' lines
    #  * `:parse_empty`: whether to parse 'e' lines
    #  * `:chunk_size`: read MAF file in chunks of this many bytes
    #  * `:random_chunk_size`: as above, but for random access ({#fetch_blocks})
    #  * `:merge_max`: merge up to this many bytes of blocks for
    #    random access
    #  * `:chunk_reader`: use the specified class to read
    #    chunks. (Only useful with {ThreadedChunkReader}).
    #  * `:threads`: number of threads to use for parallel
    #    parsing. Only useful under JRuby.
    # @api public
    
    class Parser
      include MAFParsing

      # @return [Header] header of the MAF file being parsed.
      attr_reader :header
      # @return [String] path of MAF file being parsed.
      attr_reader :file_spec
      # @return [File] file handle for MAF file.
      attr_reader :f
      # @return [StringScanner] scanner for parsing.
      attr_reader :s
      # @return [ChunkReader] ChunkReader.
      attr_reader :cr
      # @return [Boolean] whether EOF has been reached.
      attr_reader :at_end
      # @return [Hash] parser options.
      attr_reader :opts
      # @return [Integer] starting offset of the current chunk.
      attr_reader :chunk_start
      # @return [Integer] offset of the last block start in this chunk.
      attr_reader :last_block_pos

      # @api private
      attr_accessor :parse_extended
      attr_accessor :parse_empty

      SEQ_CHUNK_SIZE = 131072
      RANDOM_CHUNK_SIZE = 4096
      MERGE_MAX = SEQ_CHUNK_SIZE

      # Create a new parser instance.
      #
      # @param [String] file_spec path of file to parse.
      # @param [Hash] opts parser options.
      # @api public
      def initialize(file_spec, opts={})
        @opts = opts
        if RUBY_PLATFORM == 'java'
          opts[:threads] ||= java.lang.Runtime.runtime.availableProcessors
        end
        chunk_size = opts[:chunk_size] || SEQ_CHUNK_SIZE
        @random_access_chunk_size = opts[:random_chunk_size] || RANDOM_CHUNK_SIZE
        @merge_max = opts[:merge_max] || MERGE_MAX
        @parse_extended = opts[:parse_extended] || false
        @parse_empty = opts[:parse_empty] || false
        @chunk_start = 0
        @file_spec = file_spec
        @f = File.open(file_spec)
        reader = opts[:chunk_reader] || ChunkReader
        @cr = reader.new(@f, chunk_size)
        @s = StringScanner.new(cr.read_chunk())
        set_last_block_pos!
        @at_end = false
        _parse_header()
      end

      # Create a {ParseContext} for random access, using the given
      # chunk size.
      #
      # @return [ParseContext]
      # @api private
      def context(chunk_size)
        # IO#dup calls dup(2) internally, but seems broken on JRuby...
        fd = File.open(file_spec)
        ParseContext.new(fd, chunk_size, self, @opts)
      end

      # Execute the given block with a {ParseContext} using the given
      # `chunk_size` as an argument.
      #
      # @see #context
      # @api private
      def with_context(chunk_size)
        ctx = context(chunk_size)
        begin
          yield ctx
        ensure
          ctx.f.close
        end
      end

      # Sequence filter to apply.
      # @api public
      # @return [Hash]
      def sequence_filter
        @sequence_filter ||= {}
      end

      # Set the sequence filter.
      # @api public
      # @param [Hash] filter the new filter
      def sequence_filter=(filter)
        @sequence_filter = filter
      end

      # Fetch and parse blocks given by `fetch_list`.
      #
      # `fetch_list` should be an array of `[offset, length]` tuples.
      #
      # @param [Array] fetch_list the fetch list
      # @return [Array<Block>] the requested alignment blocks
      def fetch_blocks(fetch_list)
        merged = merge_fetch_list(fetch_list)
        if RUBY_PLATFORM == 'java' && @opts.fetch(:threads, 1) > 1
          fetch_blocks_merged_parallel(merged)
        else
          fetch_blocks_merged(merged)
        end
      end

      # Fetch and parse the blocks given by the merged fetch list.
      #
      # @param [Array] fetch_list merged fetch list from {#merge_fetch_list}.
      # @return [Array<Block>] the requested alignment blocks
      def fetch_blocks_merged(fetch_list)
        Enumerator.new do |y|
          start = Time.now
          total_size = fetch_list.collect { |e| e[1] }.reduce(:+)
          with_context(@random_access_chunk_size) do |ctx|
            fetch_list.each do |e|
              ctx.fetch_blocks(*e).each do |block|
                y << block
                #total_size += block.size
              end
            end
          end
          elapsed = Time.now - start
          rate = (total_size / 1048576.0) / elapsed
          $stderr.printf("Fetched blocks in %.3fs, %.1f MB/s.\n",
                         elapsed, rate)
        end
      end

      # Fetch and parse the blocks given by the merged fetch list, in
      # parallel. Uses the number of threads specified by the
      # `:threads` parser option.
      #
      # @param [Array] fetch_list merged fetch list from {#merge_fetch_list}.
      # @return [Array<Block>] the requested alignment blocks
      def fetch_blocks_merged_parallel(fetch_list)
        Enumerator.new do |y|
          total_size = fetch_list.collect { |e| e[1] }.reduce(:+)
          start = Time.now
          n_threads = @opts.fetch(:threads, 1)
          # TODO: break entries up into longer runs for more
          # sequential I/O
          jobs = java.util.concurrent.ConcurrentLinkedQueue.new(fetch_list)
          completed = java.util.concurrent.LinkedBlockingQueue.new(128)
          threads = []
          n_threads.times { threads << make_worker(jobs, completed) }

          n_completed = 0
          while (n_completed < fetch_list.size)
            c = completed.poll(5, java.util.concurrent.TimeUnit::SECONDS)
            if c.nil?
              if threads.find { |t| t.alive? }
                next
              else
                raise "No threads alive, completed #{n_completed}/#{fetch_list.size} jobs!"
              end
            end
            raise "worker failed: #{c}" if c.is_a? Exception
            c.each do |block|
              y << block
            end
            n_completed += 1
          end
          threads.each { |t| t.join }
          elapsed = Time.now - start
          $stderr.printf("Fetched blocks from %d threads in %.1fs.\n",
                         n_threads,
                         elapsed)
          mb = total_size / 1048576.0
          $stderr.printf("%.3f MB processed (%.1f MB/s).\n",
                         mb,
                         mb / elapsed)
        end
      end

      # Create a worker thread for parallel parsing.
      #
      # @see #fetch_blocks_merged_parallel
      def make_worker(jobs, completed)
        Thread.new do
          with_context(@random_access_chunk_size) do |ctx|
            while true
              req = jobs.poll
              break unless req
              begin
                n_blocks = req[2].size
                blocks = ctx.fetch_blocks(*req).to_a
                if blocks.size != n_blocks
                  raise "expected #{n_blocks}, got #{blocks.size}: #{e.inspect}"
                end
                completed.put(blocks)
              rescue Exception => e
                completed.put(e)
                $stderr.puts "Worker failing: #{e.class}: #{e}"
                $stderr.puts e.backtrace.join("\n")
                raise e
              end
            end
          end
        end
      end

      # Merge contiguous blocks in the given fetch list, up to
      # `:merge_max` bytes.
      #
      # Returns `[offset, size, [offset1, offset2, ...]]` tuples.
      def merge_fetch_list(orig_fl)
        fl = orig_fl.dup
        r = []
        until fl.empty? do
          cur = fl.shift
          if r.last \
            && (r.last[0] + r.last[1]) == cur[0] \
            && (r.last[1] + cur[1]) <= @merge_max
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

      # Parse the header of the MAF file.
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

      # Parse all alignment blocks until EOF.
      #
      # Delegates to {#parse_blocks_parallel} if `:threads` is set
      # under JRuby.
      #
      # @return [Enumerator<Block>] enumerator of alignment blocks.
      # @api public
      def parse_blocks
        if RUBY_PLATFORM == 'java' && @opts.has_key?(:threads)
          parse_blocks_parallel
        else
          Enumerator.new do |y|
            until at_end
              y << parse_block()
            end
          end
        end
      end

      # Parse alignment blocks with a worker thread.
      #
      # @return [Enumerator<Block>] enumerator of alignment blocks.
      # @api private
      def parse_blocks_parallel
        queue = java.util.concurrent.LinkedBlockingQueue.new(128)
        worker = Thread.new do
          begin
            until at_end
              queue.put(parse_block())
            end
            queue.put(:eof)
          rescue
            $stderr.puts "worker exiting: #{$!.class}: #{$!}"
            $stderr.puts $!.backtrace.join("\n")
          end
        end
        Enumerator.new do |y|
          saw_eof = false
          n_final_poll = 0
          while true
            block = queue.poll(1, java.util.concurrent.TimeUnit::SECONDS)
            if block == :eof
              saw_eof = true
              break
            elsif block
              y << block
            else
              # timed out
              n_final_poll += 1 unless worker.alive?
            end
            break if n_final_poll > 1
          end
          unless saw_eof
            raise "worker exited unexpectedly!"
          end
        end
      end

      def each_block
        until at_end
          yield parse_block()
        end
      end

    end

    # Exposes parser internals for unit tests.
    class DummyParser
      include MAFParsing
    end

  end
  
end
