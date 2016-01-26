require 'strscan'
require 'zlib'
require 'java' if RUBY_PLATFORM == 'java'
require 'bio-bgzf'
#121 145
# @api public
module Bio
  # @api public
  module MAF
    LOG = Bio::Log::LoggerPlus['bio-maf']
    
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
      # Currently always reads size_hint bytes.
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

    class BGZFChunkReader
      attr_reader :f, :r

      def initialize(f, _chunk_size)
        @f = f
        @r = Bio::BGZF::Reader.new(f)
      end

      def pos
        r.tell
      end

      def read_chunk
        r.read_block
      end

      def read_chunk_at(vo, _size)
        r.read_block_at(vo)
      end
    end

    class ThreadedChunkReaderWrapper

      attr_reader :cr, :pos

      def initialize(cr, buffer_size=64)
        @cr = cr
        @buffer = java.util.concurrent.LinkedBlockingQueue.new(buffer_size)
        @eof_reached = false
        @first_seq_read = false
		@restarted = false
      end
      
      # Spawn a read-ahead thread. Called from {#initialize}.
      def start_read_ahead
        LOG.debug { "Starting read-ahead thread." }
        @read_thread = Thread.new { read_ahead }
      end

      def f
        cr.f
      end
     
      # Read ahead into queue.
      def read_ahead
        # n = 0
        begin
          until f.eof?
            chunk = cr.read_chunk
            c_pos = cr.pos            
			if @restarted
				if first_chunk_read.eql?(chunk)
					@eof_reached = true
				end
			end
			@buffer.put([c_pos, chunk])
          end
          @buffer.put(:eof)
          @eof_reached = true
        rescue Exception
          @read_ahead_ex = $!
          LOG.error $!
          @buffer.put($!)
        end
      end

      def read_chunk
        if ! @first_seq_read
          # this is the first read_chunk call to read the header
          # not necessarily indicative of sequential access
          @first_seq_read = true
          chunk = cr.read_chunk
		  first_chunk_read = chunk
          @pos = cr.pos
          return chunk, first_chunk_read
        elsif @read_ahead_ex
          raise @read_ahead_ex
        elsif @eof_reached
		  restarted = true
		  cr.read_chunk_at(0, chunk_size)
		  @pos = cr.pos
		  return chunk, restarted
        else
          start_read_ahead if @read_thread.nil?
          e = @buffer.take
          case
          when e == :eof
            @eof_reached = nil
			return nil
          when e.is_a?(Exception)
            raise e
          else
            c_pos, chunk = e
            @pos = c_pos
            return chunk
          end
        end
      end

      def read_chunk_at(*args)
        cr.read_chunk_at(*args)
      end
    end

    # MAF parsing code useful for sequential and random-access parsing.
    module MAFParsing
      
      BLOCK_START = /^(?=a)/
      BLOCK_START_OR_EOS = /(?:^(?=a))|\z/
      EOL_OR_EOF = /\n|\z/
      JRUBY_P = (RUBY_PLATFORM == 'java')

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
            if opts[:strict]
              parse_error "unexpected line: '#{line}'"
            else
              LOG.warn "Ignoring invalid MAF line: '#{line}'"
            end
          end
        end
        b = Block.new(block_vars,
                      seqs,
                      block_offset,
                      s.pos - block_start_pos,
                      filtered)
        if opts[:retain_text]
          b.orig_text = s.string.slice(block_start_pos...(s.pos))
        end
        return b
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
        @cr = parser.base_reader.new(@f, chunk_size)
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
          LOG.debug { "fetching blocks from #{offset} (length #{len}): #{block_offsets.inspect}" }
          start_chunk_read_if_needed(offset, len)
          # read chunks until we have the entire merged set of
          # blocks ready to parse
          # to avoid fragment joining
          append_chunks_to(len)
          # parse the blocks
          block_offsets.each do |expected_offset|
            # skip ahead, in case there is a gap resulting from a
            # block that is not being parsed
            rel_offset = expected_offset - offset
            if s.pos < rel_offset
              s.pos = rel_offset
            end
            # now actually parse the block data
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
    #  * `:join_blocks`: join blocks where possible
    #  * `:upcase`: fold sequence data to upper case
    #  * `:chunk_size`: read MAF file in chunks of this many bytes
    #  * `:random_chunk_size`: as above, but for random access ({#fetch_blocks})
    #  * `:merge_max`: merge up to this many bytes of blocks for
    #    random access
    #  * `:threads`: number of threads to use for parallel
    #    parsing. Only useful under JRuby.
    #  * `:strict`: abort on un-parseable lines instead of continuing with
    #    a warning.
    # @api public
    
    class Parser
      include MAFParsing

      # @return [Header] header of the MAF file being parsed.
      attr_reader :header
      # @return [String] path of MAF file being parsed.
      attr_reader :file_spec
      # @return [IO] file handle for MAF file.
      attr_reader :f
      # May be gzip-compressed.
      # @return [IO] file handle for physical MAF file.
      # @api private
      attr_reader :phys_f
      # @return [StringScanner] scanner for parsing.
      attr_reader :s
      # @return [ChunkReader] ChunkReader.
      attr_reader :cr
      # @return [Class] ChunkReader class to use for random access
      # @see ParseContext
      attr_reader :base_reader
      # @return [Boolean] whether EOF has been reached.
      attr_reader :at_end
      # @return [Hash] parser options.
      attr_reader :opts
      # @return [Integer] starting offset of the current chunk.
      attr_reader :chunk_start
      # @return [Integer] offset of the last block start in this chunk.
      attr_reader :last_block_pos
      # @return [Symbol] compression method used for this file, or nil
      attr_reader :compression

      # @api private
      attr_accessor :parse_extended
      attr_accessor :parse_empty

      SEQ_CHUNK_SIZE = 131072
      RANDOM_CHUNK_SIZE = 4096
      MERGE_MAX = SEQ_CHUNK_SIZE

      DEFAULT_OPTS = {
        :chunk_size => SEQ_CHUNK_SIZE,
        :random_chunk_size => RANDOM_CHUNK_SIZE,
        :merge_max => MERGE_MAX,
        :parse_extended => false,
        :parse_empty => false,
        :readahead_thread => true,
        :seq_parse_thread => true
      }
      if JRUBY_P
        DEFAULT_OPTS[:threads] = java.lang.Runtime.runtime.availableProcessors
      end

      # Create a new parser instance.
      #
      # @param [String] file_spec path of file to parse.
      # @param [Hash] parse_opts parser options.
      # @api public
      def initialize(file_spec, parse_opts={})
        opts = DEFAULT_OPTS.merge(parse_opts)
        @opts = opts
        @random_access_chunk_size = opts[:random_chunk_size]
        @merge_max = opts[:merge_max]
        @parse_extended = opts[:parse_extended]
        @parse_empty = opts[:parse_empty]
        @chunk_start = 0
        if file_spec.respond_to? :flush
          # an IO object
          # guess what, Pathnames respond to :read...
          @f = file_spec
          @file_spec = @f.path if @f.respond_to?(:path)
          # TODO: test for gzip?
        else
          # a pathname (or Pathname)
          @file_spec = file_spec
          @phys_f = File.open(file_spec)
          if file_spec.to_s.end_with?(".maf.gz")
            @f = Zlib::GzipReader.new(@phys_f)
            @compression = :gzip
          else
            @f = @phys_f
          end
        end
        if @file_spec.to_s =~ /\.bgzf?$/
          @base_reader = BGZFChunkReader
          @compression = :bgzf
        else
          @base_reader = ChunkReader
        end
        @cr = base_reader.new(@f, opts[:chunk_size])
        if JRUBY_P && opts[:readahead_thread]
          LOG.debug "Using ThreadedChunkReaderWrapper."
          @cr = ThreadedChunkReaderWrapper.new(@cr)
        end
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
        if file_spec
          fd = File.open(file_spec)
        else
          fd = f.dup
        end
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
          if JRUBY_P && @opts.fetch(:threads, 1) > 1
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
        count = 0
        with_context(@random_access_chunk_size) do |ctx|
          fetch_list.each do |e|
            ctx.fetch_blocks(*e, &blk)
            count += 1
          end
        end
        elapsed = Time.now - start
        rate = (total_size / 1048576.0) / elapsed
        LOG.debug { sprintf("Fetched %d blocks in %.3fs, %.1f MB/s.",
                            count, elapsed, rate) }
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
        LOG.debug { sprintf("Fetched blocks from %d threads in %.1fs.",
                            n_threads,
                            elapsed) }
        mb = total_size / 1048576.0
        LOG.debug { sprintf("%.3f MB processed (%.1f MB/s).",
                            mb,
                            mb / elapsed) }
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
            LOG.error "Worker failing: #{e.class}: #{e}"
            LOG.error e
            raise e
          end
        end
      end

      # Merge contiguous blocks in the given fetch list, up to
      # `:merge_max` bytes.
      #
      # Returns `[offset, size, [offset1, offset2, ...]]` tuples.
      def merge_fetch_list(orig_fl)
        case compression
        when nil
          _merge_fetch_list(orig_fl)
        when :bgzf
          _merge_bgzf_fetch_list(orig_fl)
        end
      end

      def _merge_fetch_list(orig_fl)
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

      # Build a merged fetch list in a BGZF-aware way.  This will
      # group together all MAF blocks from a single BGZF block. These
      # MAF blocks may not be consecutive.
      def _merge_bgzf_fetch_list(orig_fl)
        block_e = orig_fl.chunk { |entry|
          Bio::BGZF::vo_block_offset(entry[0])
        }
        block_e.collect do |bgzf_block, fl|
          # text size to read from disk, from the start of the first
          # block to the end of the last block
          text_size = fl.last[0] + fl.last[1] - fl.first[0]
          offsets = fl.collect { |e| e[0] }
          [fl.first[0], text_size, offsets]
        end
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
        if ! s.skip_until(BLOCK_START)
          @at_end = true
        end
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
          if JRUBY_P && opts[:seq_parse_thread]
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

      WRAP_OPTS = [:as_bio_alignment, :join_blocks, :remove_gaps, :upcase]

      def wrap_block_seq(fun, &blk)
        opts = WRAP_OPTS.find_all { |o| @opts[o] }
        opts << :sequence_filter if sequence_filter && (! sequence_filter.empty?)
        LOG.debug { "wrapping #{fun} with #{opts.inspect}" }
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
        when :upcase
          conv_send(options,
                    fun,
                    :upcase!,
                    true,
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

      def conv_send(options, search, sym, always_yield_block=false)
        _wrap(options, search) do |block|
          v = block.send(sym)
          if always_yield_block
            yield block
          else
            yield v if v
          end
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
            LOG.debug "Starting parse worker."
            until at_end
              block = _parse_block()
              queue.put(block) if block
            end
            queue.put(:eof)
            LOG.debug { "Parse worker reached EOF." }
          rescue Exception
            LOG.error $!
            Thread.current[:exception] = $!
            raise
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
            unless worker.alive?
              LOG.debug "Worker has exited."
              n_final_poll += 1
            end
          end
          break if n_final_poll > 1
        end
        unless saw_eof
          raise "worker exited unexpectedly from #{worker[:exception]}!"
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

    def handle_logging_options(opts)
      opts.on("--logger filename", String,
              "Log to file (default STDERR)") do |name|
        Bio::Log::CLI.logger(name)
      end
      opts.on("--trace options", String,
              "Set log level",
              "(default INFO, see bio-logger)") do |s|
        Bio::Log::CLI.trace(s)
      end
      opts.on("-q", "--quiet", "Run quietly") do
        Bio::Log::CLI.trace('error')
      end
      opts.on("-v", "--verbose", "Run verbosely") do
        Bio::Log::CLI.trace('info')
      end
      opts.on("--debug", "Run with extra debugging output") do
        Bio::Log::CLI.trace('debug')
      end
    end
    module_function :handle_logging_options

  end
  
end
