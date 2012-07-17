require 'strscan'
require 'java' if RUBY_PLATFORM == 'java'

# @api public
module Bio
  # @api public
  module MAF

    # @api public
    class ParseError < Exception; end

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
      def _parse_block
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
        filtered = false
        lines = payload.split("\n")
        until lines.empty?
          line = lines.shift
          first = line.getbyte(0)
          if first == S
            seq = parse_seq_line(line, sequence_filter)
            if seq
              seqs << seq
            else
              filtered = true
            end
          elsif first == E && parse_empty
            e_seq = parse_empty_line(line, sequence_filter)
            if e_seq
              seqs << e_seq
            else
              filtered = true
            end
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
        Block.new(block_vars,
                  seqs,
                  block_offset,
                  s.pos - block_start_pos,
                  filtered)
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
      attr_accessor :f, :s, :cr, :parser, :opts
      attr_accessor :chunk_start, :last_block_pos, :at_end

      def initialize(fd, chunk_size, parser)
        @f = fd
        @parser = parser
        @opts = parser.opts
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
        if block_given?
          start_chunk_read_if_needed(offset, len)
          # read chunks until we have the entire merged set of
          # blocks ready to parse
          # to avoid fragment joining
          append_chunks_to(len)
          # parse the blocks
          block_offsets.each do |expected_offset|
            block = _parse_block
            parse_error("expected a block at offset #{expected_offset} but could not parse one!") unless block
            parse_error("got block with offset #{block.offset}, expected #{expected_offset}!") unless block.offset == expected_offset
            yield block
          end
        else
          enum_for(:fetch_blocks, offset, len, block_offsets)
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
    #  * `:remove_gaps`: remove gaps left after filtering sequences
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

      def close
        f.close
      end

      # Create a {ParseContext} for random access, using the given
      # chunk size.
      #
      # @return [ParseContext]
      # @api private
      def context(chunk_size)
        # IO#dup calls dup(2) internally, but seems broken on JRuby...
        fd = File.open(file_spec)
        ParseContext.new(fd, chunk_size, self)
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
      # @yield [block] each block matched, in turn
      # @return [Enumerable<Block>] each matching {Block}, if no block given
      def fetch_blocks(fetch_list, &blk)
        if blk
          merged = merge_fetch_list(fetch_list)
          if RUBY_PLATFORM == 'java' && @opts.fetch(:threads, 1) > 1
            fun = lambda { |&b2| fetch_blocks_merged_parallel(merged, &b2) }
          else
            fun = lambda { |&b2| fetch_blocks_merged(merged, &b2) }
          end
          wrap_block_seq(fun, &blk)
        else
          enum_for(:fetch_blocks, fetch_list)
        end
      end

      # Fetch and parse the blocks given by the merged fetch list.
      #
      # @param [Array] fetch_list merged fetch list from {#merge_fetch_list}.
      # @return [Array<Block>] the requested alignment blocks
      def fetch_blocks_merged(fetch_list, &blk)
        start = Time.now
        total_size = fetch_list.collect { |e| e[1] }.reduce(:+)
        with_context(@random_access_chunk_size) do |ctx|
          fetch_list.each do |e|
            ctx.fetch_blocks(*e, &blk)
          end
        end
        elapsed = Time.now - start
        # TODO: debug log
        # rate = (total_size / 1048576.0) / elapsed
        # $stderr.printf("Fetched blocks in %.3fs, %.1f MB/s.\n",
        #                elapsed, rate)
      end

      # Fetch and parse the blocks given by the merged fetch list, in
      # parallel. Uses the number of threads specified by the
      # `:threads` parser option.
      #
      # @param [Array] fetch_list merged fetch list from {#merge_fetch_list}.
      # @return [Array<Block>] the requested alignment blocks
      def fetch_blocks_merged_parallel(fetch_list)
        total_size = fetch_list.collect { |e| e[1] }.reduce(:+)
        start = Time.now
        n_threads = @opts.fetch(:threads, 1)
        # TODO: break entries up into longer runs for more
        # sequential I/O
        jobs = java.util.concurrent.ConcurrentLinkedQueue.new(fetch_list)
        ct = CompletionTracker.new(fetch_list)
        completed = ct.queue
        threads = []
        n_threads.times { threads << make_worker(jobs, ct) }

        n_res = 0
        while n_res < fetch_list.size
          c = completed.poll(1, java.util.concurrent.TimeUnit::SECONDS)
          unless c
            raise "Worker failed!" if threads.find { |t| t.status.nil? }
            next
          end
          c.each do |block|
            yield block
          end
          n_res += 1
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

      # Create a worker thread for parallel parsing.
      #
      # @see #fetch_blocks_merged_parallel
      def make_worker(jobs, ct)
        Thread.new do
          begin
            with_context(@random_access_chunk_size) do |ctx|
              while true
                req = jobs.poll
                break unless req
                n_blocks = req[2].size
                blocks = ctx.fetch_blocks(*req).to_a
                if blocks.size != n_blocks
                  raise "expected #{n_blocks}, got #{blocks.size}: #{e.inspect}"
                end
                ct << blocks
              end
            end
          rescue Exception => e
            $stderr.puts "Worker failing: #{e.class}: #{e}"
            $stderr.puts e.backtrace.join("\n")
            raise e
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
      # @return [Enumerator<Block>] enumerator of {Block}s if no block given.
      # @yield [block] Passes each {Block} in turn to a block
      # @api public
      def each_block(&blk)
        if block_given?
          if RUBY_PLATFORM == 'java' && @opts.has_key?(:threads)
            fun = method(:parse_blocks_parallel)
          else
            fun = method(:each_block_seq)
          end
          wrap_block_seq(fun, &blk)
        else
          enum_for(:each_block)
        end
      end
      alias_method :parse_blocks, :each_block

      def each_block_seq
        until at_end
          block = _parse_block()
          yield block if block
        end
      end

      def parse_block
        b = nil
        wrap_block_seq(lambda { |&blk| blk.call(_parse_block()) }) do |block|
          b = block
        end
        b
      end

      WRAP_OPTS = [:as_bio_alignment, :join_blocks, :remove_gaps]

      def wrap_block_seq(fun, &blk)
        opts = WRAP_OPTS.find_all { |o| @opts[o] }
        opts << :sequence_filter if sequence_filter && (! sequence_filter.empty?)
        _wrap(opts, fun, &blk)
      end

      # options should be [:outer, ..., :inner]
      def _wrap(options, fun, &blk)
        first = options.shift
        case first
        when nil
          fun.call(&blk)
        when :sequence_filter
          conv_map(options,
                   fun,
                   lambda { |b| b if b.sequences.size > 1 },
                   &blk)
        when :join_blocks
          block_joiner(options, fun, &blk)
        when :as_bio_alignment
          conv_send(options,
                    fun,
                    :to_bio_alignment,
                    &blk)
        when :remove_gaps
          conv_map(options,
                   fun,
                   lambda { |b| b.remove_gaps! if b.filtered?; b },
                   &blk)
        else
          raise "unhandled wrapper mode: #{first}"
        end
      end

      def filter_seq_count(fun)
        fun.call() do |block|
          yield block if block.filtered? && block.sequences.size > 1
        end
      end

      def block_joiner(options, fun)
        prev = nil
        _wrap(options, fun) do |cur|
          if prev && (prev.filtered? || cur.filtered?) \
            && prev.joinable_with?(cur)
            prev = prev.join(cur)
          else
            yield prev if prev
            prev = cur
          end
        end
        yield prev if prev
      end

      def conv_map(options, search, fun)
        _wrap(options, search) do |block|
          v = fun.call(block)
          yield v if v
        end
      end

      def conv_send(options, search, sym)
        _wrap(options, search) do |block|
          v = block.send(sym)
          yield v if v
        end
      end

      # Parse alignment blocks with a worker thread.
      #
      # @block block handler
      # @api private
      def parse_blocks_parallel
        queue = java.util.concurrent.LinkedBlockingQueue.new(128)
        worker = Thread.new do
          begin
            until at_end
              block = _parse_block()
              queue.put(block) if block
            end
            queue.put(:eof)
          rescue
            $stderr.puts "worker exiting: #{$!.class}: #{$!}"
            $stderr.puts $!.backtrace.join("\n")
          end
        end
        saw_eof = false
        n_final_poll = 0
        while true
          block = queue.poll(1, java.util.concurrent.TimeUnit::SECONDS)
          if block == :eof
            saw_eof = true
            break
          elsif block
            yield block
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

    class CompletionTracker
      attr_reader :queue, :offsets, :delayed

      def initialize(fetch_list)
        @offsets = fetch_list.collect { |e| e[0] }
        @queue = java.util.concurrent.LinkedBlockingQueue.new(128)
        @delayed = {}
        @sem = Mutex.new
      end

      def next_expected
        offsets.first
      end

      def <<(blocks)
        @sem.synchronize do
          f_offset = blocks.first.offset
          if f_offset == next_expected
            offsets.shift
            queue.put(blocks)
            drain_delayed
          else
            # out of order
            delayed[f_offset] = blocks
          end
        end
      end

      def drain_delayed
        while e = delayed.delete(next_expected)
          offsets.shift
          queue.put(e)
        end
      end
    end

    # Exposes parser internals for unit tests.
    class DummyParser
      include MAFParsing
    end

  end
  
end
