require 'kyotocabinet'
require 'jruby/profiler' if RUBY_PLATFORM == 'java'

#require 'bio-ucsc-api'
require 'bio-genomic-interval'

module Bio

  module MAF

    # Binary record packing and unpacking.
    # @api private
    module KVHelpers

      KEY = Struct.new([[:marker,    :uint8],
                        [:seq_id,    :uint8],
                        [:bin,       :uint16],
                        [:seq_start, :uint32],
                        [:seq_end,   :uint32]])

      VAL = Struct.new([[:offset,      :uint64],
                        [:length,      :uint32],
                        [:text_size,   :uint32],
                        [:n_seq,       :uint8],
                        [:species_vec, :uint64]])

      KEY_FMT = KEY.fmt
      KEY_SCAN_FMT = KEY.extractor_fmt(:seq_id, :bin, :seq_start, :seq_end)
      CHROM_BIN_PREFIX_FMT = KEY.extractor_fmt(:marker, :seq_id, :bin)

      VAL_FMT = VAL.fmt
      VAL_IDX_OFFSET_FMT = VAL.extractor_fmt(:offset, :length)
      VAL_TEXT_SIZE_FMT = VAL.extractor_fmt(:text_size)
      VAL_N_SEQ_FMT = VAL.extractor_fmt(:n_seq)
      VAL_SPECIES_FMT = VAL.extractor_fmt(:species_vec)

      module_function

      def extract_species_vec(entry)
        entry[1].unpack(VAL_SPECIES_FMT)[0]
      end

      def extract_n_sequences(entry)
        entry[1].unpack(VAL_N_SEQ_FMT)[0]
      end

      def extract_index_offset(entry)
        entry[1].unpack(VAL_IDX_OFFSET_FMT)
      end

      def extract_text_size(entry)
        entry[1].unpack(VAL_TEXT_SIZE_FMT)[0]
      end

      def unpack_key(ks)
        ks.unpack(KEY_FMT)
      end

      def bin_start_prefix(chrom_id, bin)
        [0xFF, chrom_id, bin].pack(CHROM_BIN_PREFIX_FMT)
      end
    end

    class Access

      attr_accessor :parse_options

      def self.maf_dir(dir, options={})
        self.new(dir, options)
      end

      def initialize(dir, options={})
        @dir = dir
        @parse_options = options
        @maf_files = Dir.glob("#{dir}/*.maf")
        @indices = {}
        @maf_by_chrom = {}
        scan_indices!
      end

      def close
        @indices.values.each { |ki| ki.close }
      end

      def scan_indices!
        @maf_files.each do |maf|
          base = File.basename(maf, '.maf')
          index_f = "#{@dir}/#{base}.kct"
          if File.exists? index_f
            index = KyotoIndex.open(index_f)
            @indices[index.ref_seq] = index
            @maf_by_chrom[index.ref_seq] = maf
          end
        end
      end

      def find(intervals, &blk)
        if block_given?
          by_chrom = intervals.group_by { |i| i.chrom }
          by_chrom.keys.each do |chrom|
            unless @indices.has_key? chrom
              raise "No index available for chromosome #{chrom}!"
            end
          end
          by_chrom.each do |chrom, c_intervals|
            index = @indices[chrom]
            with_parser(chrom) do |parser|
              index.find(c_intervals, parser, &blk)
            end
          end
        else
          enum_for(:find, intervals)
        end
      end

      def with_parser(chrom)
        parser = Parser.new(@maf_by_chrom[chrom], @parse_options)
        begin
          yield parser
        ensure
          parser.close
        end
      end

    end

    class KyotoIndex
      include KVHelpers

      attr_reader :db, :species, :species_max_id
      attr_accessor :index_sequences, :ref_seq

      FORMAT_VERSION_KEY = 'bio-maf:index-format-version'
      FORMAT_VERSION = 2
      REF_SEQ_KEY = 'bio-maf:reference-sequence'
      MAX_SPECIES = 64

      ## Key-value store index format
      ##
      ## This format is designed for Kyoto Cabinet but should work on
      ## other key-value databases allowing binary data.
      ##
      ## Index metadata is stored as ASCII text, but index data is
      ## stored as packed binary values.
      ##
      ## Index metadata:
      ##
      ##   Sequence IDs:
      ##     sequence:<name> => <id>
      ##
      ##     Each indexed sequence has a corresponding entry of this
      ##     kind. The <name> parameter is the sequence or chromosome
      ##     name as found in the MAF file, e.g. mm8.chr7. The <id>
      ##     parameter is assigned when the sequence is indexed, and
      ##     can be from 0 to 255.
      ##
      ##   Species IDs:
      ##     species:<name> => <id>
      ##
      ##     Each indexed species has a corresponding entry of this
      ##     kind. The <name> parameter is the species part of the
      ##     sequence name as found in the MAF file, e.g. 'mm8' for
      ##     'mm8.chr7'. The <id> parameter is assigned when the
      ##     species is indexed, and can be from 0 to 255.
      ##
      ## Index data:
      ##
      ##   For each sequence upon which an index is built, one index
      ##   entry is generated per MAF alignment block. The key
      ##   identifies the sequence, the UCSC index bin, and the
      ##   zero-based start and end positions of the sequence. The
      ##   value gives the offset and size of the alignment block
      ##   within the MAF file.
      ##
      ##   All values are stored as big-endian, unsigned packed binary
      ##   data.
      ##
      ## Keys: (12 bytes) [CCS>L>L>]
      ##
      ##   0xFF (1 byte):
      ##      index entry prefix
      ##   Sequence chromosome ID (1 byte):
      ##      corresponds to sequence:<name> entries
      ##   UCSC bin (16 bits)
      ##   Sequence start, zero-based, inclusive (32 bits)
      ##   Sequence end, zero-based, exclusive (32 bits)
      ##
      ## Values (25 bytes) [Q>L>L>CQ>]
      ##
      ##   MAF file offset (64 bits)
      ##   MAF alignment block length (32 bits)
      ##   Block text size (32 bits)
      ##   Number of sequences in block (8 bits)
      ##   Species bit vector (64 bits)
      ##
      ## Example:
      ##
      ##  For a block with sequence 0, bin 1195, start 80082334, end
      ##       80082368, MAF offset 16, and MAF block length 1087:
      ##
      ##     |  |id| bin | seq_start | seq_end   |
      ## key: FF 00 04 AB 04 C5 F5 9E 04 C5 F5 C0
      ##
      ##     |         offset        |  length   |   ts   |ns|  species_vec  |
      ## val: 00 00 00 00 00 00 00 10 00 00 04 3F  [TODO]

      #### Public API

      # Open an existing index for reading.
      # @param [String] path path to existing Kyoto Cabinet index
      # @return [KyotoIndex]
      def self.open(path)
        return KyotoIndex.new(path)
      end

      # Build a new index from the MAF file being parsed by `parser`,
      # and store it in `path`.
      # @param [Parser] parser MAF parser for file to index
      # @param [String] path path to index file to create
      # @return [KyotoIndex]
      def self.build(parser, path)
        idx = self.new(path)
        idx.build_default(parser)
        return idx
      end

      # Find all alignment blocks in the genomic regions in the list
      # of Bio::GenomicInterval objects, and parse them with the given
      # parser.
      #
      # An optional Hash of filters may be passed in. The following
      # keys are used:
      #
      #  * `:with_all_species => ["sp1", "sp2", ...]`
      #
      #      Only match alignment blocks containing all given species.
      #
      #  * `:at_least_n_sequences => n`
      #
      #      Only match alignment blocks with at least N sequences.
      #
      #  * `:min_size => n`
      #
      #      Only match alignment blocks with text size at least N.
      #
      #  * `:max_size => n`
      #
      #      Only match alignment blocks with text size at most N.
      #
      # @param [Enumerable<Bio::GenomicInterval>] intervals genomic
      #  intervals to parse.
      # @param [Parser] parser MAF parser for file to fetch blocks
      #  from.
      # @param [Hash] filter Block filter expression.
      # @yield [block] each {Block} matched, in turn
      # @return [Enumerable<Block>] each matching {Block}, if no block given
      # @api public
      def find(intervals, parser, filter={}, &blk)
        # start = Time.now
        fl = fetch_list(intervals, filter)
        # $stderr.printf("Built fetch list of %d items in %.3fs.\n",
        #                fl.size,
        #                Time.now - start)
        parser.fetch_blocks(fl, &blk)
      end

      # Close the underlying Kyoto Cabinet database handle.
      def close
        db.close
      end

      #### KyotoIndex Internals
      # @api private

      def initialize(path, db_arg=nil)
        @species = {}
        @species_max_id = -1
        if db_arg || ((path.size > 1) and File.exist?(path))
          mode = KyotoCabinet::DB::OREADER
        else
          mode = KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE
        end
        @db = db_arg || KyotoCabinet::DB.new
        @path = path
        unless db_arg || db.open(path.to_s, mode)
          raise "Could not open DB file!"
        end
        if mode == KyotoCabinet::DB::OREADER
          self.ref_seq = db[REF_SEQ_KEY]
          load_index_sequences
          load_species
        end
      end

      # Reopen the same DB handle read-only. Only useful for unit tests.
      def reopen
        KyotoIndex.new(@path, @db)
      end

      def dump(stream=$stdout)
        stream.puts "KyotoIndex dump: #{@path}"
        stream.puts
        if db.count == 0
          stream.puts "Empty database!"
          return
        end
        db.cursor_process do |cur|
          stream.puts "== Metadata =="
          cur.jump('')
          while true
            k, v = cur.get(false)
            raise "unexpected end of records!" unless k
            break if k[0] == "\xff"
            stream.puts "#{k}: #{v}"
            unless cur.step
              raise "could not advance cursor!"
            end
          end
          stream.puts "== Index records =="
          while pair = cur.get(true)
            _, chr, bin, s_start, s_end = pair[0].unpack(KEY_FMT)
            offset, len, text_size, n_seq, species_vec = pair[1].unpack(VAL_FMT)
            stream.puts "#{chr} [bin #{bin}] #{s_start}:#{s_end}"
            stream.puts "  offset #{offset}, length #{len}"
            stream.puts "  text size: #{text_size}"
            stream.puts "  sequences in block: #{n_seq}"
            stream.printf("  species vector: %016x\n", species_vec)
          end
        end
      end

      ## Retrieval:
      ##  1. merge the intervals of interest
      ##  2. for each interval, compute the bins with #bin_all
      ##  3. for each bin to search, make a list of intervals of
      ##     interest
      ##  4. compute the spanning interval for that bin
      ##  5. start at the beginning of the bin
      ##  6. if a record intersects the spanning interval: 
      ##    A. #find an interval it intersects
      ##    B. if found, add to the fetch list
      ##  7. if a record starts past the end of the spanning interval,
      ##     we are done scanning this bin.
      ##
      ## Optimizations:
      ##  * once we reach the start of the spanning interval,
      ##    all records start in it until we see a record starting
      ##    past it.
      ##  * as record starts pass the start of intervals of interest,
      ##    pull those intervals off the list

      # Build a fetch list of alignment blocks to read, given an array
      # of Bio::GenomicInterval objects
      def fetch_list(intervals, filter_spec={})
        filter_spec ||= {}
        filters = Filters.build(filter_spec, self)
        chrom = intervals.first.chrom
        chrom_id = index_sequences[chrom]
        unless chrom_id
          raise "chromosome #{chrom} not indexed!"
        end
        if intervals.find { |i| i.chrom != chrom }
          raise "all intervals must be for the same chromosome!"
        end
        # for each bin, build a list of the intervals to look for there
        bin_intervals = Hash.new { |h, k| h[k] = [] }
        intervals.each do |i|
          i.bin_all.each do |bin|
            bin_intervals[bin] << (i.zero_start...i.zero_end)
          end
        end
        bin_intervals.values.each do |intervals|
          intervals.sort_by! {|i| i.begin}
        end
        matches = if RUBY_PLATFORM == 'java' && bin_intervals.size > 4
                    scan_bins_parallel(chrom_id, bin_intervals, filters)
                  else
                    scan_bins(chrom_id, bin_intervals, filters)
                  end
        matches.sort_by! { |e| e[0] } # sort by offset in file
      end # #fetch_list

      # Scan the index for blocks matching the given bins and intervals.
      def scan_bins(chrom_id, bin_intervals, filters)
        to_fetch = []
        db.cursor_process do |cur|
          bin_intervals.each do |bin, bin_intervals_raw|
            matches = scan_bin(cur, chrom_id, bin, bin_intervals_raw, filters)
            to_fetch.concat(matches)
          end 
        end
        to_fetch
      end

      def with_profiling
        if RUBY_PLATFORM == 'java' && ENV['profile']
          rv = nil
          pdata = JRuby::Profiler.profile do
            rv = yield
          end
          printer = JRuby::Profiler::FlatProfilePrinter.new(pdata)
          printer.printProfile(STDERR)
          return rv
        else
          yield
        end
      end

      def scan_bins_parallel(chrom_id, bin_intervals, filters)
        start = Time.now
        n_threads = ENV['profile'] ? 1 : java.lang.Runtime.runtime.availableProcessors
        jobs = java.util.concurrent.ConcurrentLinkedQueue.new(bin_intervals.to_a)
        completed = java.util.concurrent.LinkedBlockingQueue.new(128)
        threads = []
        n_threads.times do
          threads << make_scan_worker(jobs, completed) do |cur, req|
            bin, intervals = req
            scan_bin(cur, chrom_id, bin, intervals, filters)
          end
        end
        n_completed = 0
        to_fetch = []
        while (n_completed < bin_intervals.size)
          c = completed.poll(5, java.util.concurrent.TimeUnit::SECONDS)
          if c.nil?
            if threads.find { |t| t.alive? }
              next
            else
              raise "No threads alive, completed #{n_completed}/#{bin_intervals.size} jobs!"
            end
          end
          raise "worker failed: #{c}" if c.is_a? Exception
          to_fetch.concat(c)
          n_completed += 1
        end
        threads.each { |t| t.join }
        $stderr.printf("Matched %d index records with %d threads in %.3f seconds.\n",
                       to_fetch.size, n_threads, Time.now - start)
        to_fetch
      end

      def make_scan_worker(jobs, completed)
        Thread.new do
          with_profiling do
            db.cursor_process do |cur|
              while true
                req = jobs.poll
                break unless req
                begin
                  result = yield(cur, req)
                  completed.put(result)
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
      end

      def scan_bin(cur, chrom_id, bin, bin_intervals, filters)
        # bin_intervals is sorted by zero_start
        # compute the start and end of all intervals of interest
        spanning_start = bin_intervals.first.begin
        spanning_end = bin_intervals.map {|i| i.end}.max
        # scan from the start of the bin
        cur.jump(bin_start_prefix(chrom_id, bin))
        matches = []
        while pair = cur.get(true)
          c_chr, c_bin, c_start, c_end = pair[0].unpack(KEY_SCAN_FMT)
          if (c_chr != chrom_id) \
            || (c_bin != bin) \
            || c_start >= spanning_end
            # we've hit the next bin, or chromosome, or gone past
            # the spanning interval, so we're done with this bin
            break
          end
          if c_end >= spanning_start # possible overlap
            # any intervals that end before the start of the current
            # block are no longer relevant
            while bin_intervals.first.end < c_start
              bin_intervals.shift
            end
            bin_intervals.each do |i|
              i_start = i.begin
              break if i_start > c_end
              if ((c_start <= i_start && i_start < c_end) \
                  || i.include?(c_start)) \
                  && filters.match(pair)
                # match
                matches << extract_index_offset(pair)
                break
              end
            end
          end
        end
        matches
      end

      def overlaps?(gi, i_start, i_end)
        g_start = gi.begin

        (i_start <= g_start && g_start < i_end) \
         || gi.include?(i_start)
      end

      def build_default(parser)
        first_block = parser.parse_block
        self.ref_seq = first_block.sequences.first.source
        db[REF_SEQ_KEY] = ref_seq
        db[FORMAT_VERSION_KEY] = FORMAT_VERSION
        @index_sequences = { ref_seq => 0 }
        store_index_sequences!
        index_blocks([first_block])
        parser.enum_for(:each_block).each_slice(1000).each do |blocks|
          index_blocks(blocks)
        end
        db.synchronize(true)
      end

      def index_blocks(blocks)
        h = blocks.map { |b| entries_for(b) }.reduce(:merge!)
        db.set_bulk(h, false)
      end

      def load_index_sequences
        h = {}
        db.match_prefix("sequence:").each do |key|
          _, name = key.split(':', 2)
          id = db[key].to_i
          h[name] = id
        end
        @index_sequences = h
      end

      def store_index_sequences!
        index_sequences.each do |name, id|
          db.set("sequence:#{name}", id.to_s)
        end
      end

      def load_species
        db.match_prefix("species:").each do |key|
          _, name = key.split(':', 2)
          id = db[key].to_i
          @species[name] = id
        end
        @species_max_id = @species.values.sort.last || -1
      end

      def species_id_for_seq(seq)
        # NB can have multiple dots
        # example: otoGar1.scaffold_104707.1-93001
        parts = seq.split('.', 2)
        if parts.size == 2
          species_name = parts[0]
          if species.has_key? species_name
            return species[species_name]
          else
            species_id = @species_max_id + 1
            if species_id >= MAX_SPECIES
              raise "cannot index MAF file with more than #{MAX_SPECIES} species"
            end
            species[species_name] = species_id
            db["species:#{species_name}"] = species_id
            @species_max_id = species_id
            return species_id
          end
        else
          # not in species.sequence format, apparently
          return nil
        end
      end

      def build_block_value(block)
        bits = block.sequences.collect {|s| 1 << species_id_for_seq(s.source) }
        vec = bits.reduce(0, :|)
        return [block.offset,
                block.size,
                block.text_size,
                block.sequences.size,
                vec].pack(VAL_FMT)
      end

      def entries_for(block)
        unless block.ref_seq.source == @ref_seq
          raise "Inconsistent reference sequence: expected #{@ref_seq}, got #{block.ref_seq.source}"
        end
        h = {}
        val = build_block_value(block)
        block.sequences.each do |seq|
          seq_id = index_sequences[seq.source]
          next unless seq_id
          seq_end = seq.start + seq.size
          bin = Bio::Ucsc::UcscBin.bin_from_range(seq.start, seq_end)
          key = [255, seq_id, bin, seq.start, seq_end].pack(KEY_FMT)
          h[key] = val
        end
        return h
      end
    end # class KyotoIndex

    class Filter
      include KVHelpers

      def call(e)
        match(e)
      end
    end

    class AllSpeciesFilter < Filter
      attr_reader :bs
      def initialize(species, idx)
        ids = species.collect {|s| 1 << idx.species.fetch(s) }
        @mask = ids.reduce(0, :|)
      end

      def match(entry)
        vec = extract_species_vec(entry)
        (@mask & vec) == @mask
      end
    end

    class AtLeastNSequencesFilter < Filter
      attr_reader :n
      def initialize(n, idx)
        @n = n
      end

      def match(entry)
        extract_n_sequences(entry) >= @n
      end
    end

    class MaxSizeFilter < Filter
      def initialize(n, idx)
        @n = n
      end
      def match(entry)
        extract_text_size(entry) <= @n
      end
    end

    class MinSizeFilter < Filter
      def initialize(n, idx)
        @n = n
      end
      def match(entry)
        extract_text_size(entry) >= @n
      end
    end

    class Filters
      include KVHelpers

      FILTER_CLASSES = {
        :with_all_species => MAF::AllSpeciesFilter,
        :at_least_n_sequences => MAF::AtLeastNSequencesFilter,
        :min_size => MAF::MinSizeFilter,
        :max_size => MAF::MaxSizeFilter
      }

      def self.build(spec, idx)
        l = spec.collect do |key, val|
          if FILTER_CLASSES.has_key? key
            FILTER_CLASSES[key].new(val, idx)
          else
            raise "Unsupported filter key #{key}!"
          end
        end
        return Filters.new(l)
      end

      def initialize(l)
        @l = l
      end

      def match(entry)
        return ! @l.find { |f| ! f.call(entry) }
      end
    end

  end # module MAF
  
end
