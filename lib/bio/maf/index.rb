require 'dbi'
require 'kyotocabinet'

require 'bio-ucsc-api'
require 'bio-genomic-interval'

module Bio

  module MAF

    class KyotoIndex

      attr_reader :db
      attr_accessor :index_sequences

      KEY_FMT = "CCS>L>L>"
      KEY_SCAN_FMT = "xCS>L>L>"
      CHROM_BIN_PREFIX_FMT = "CCS>"
      VAL_FMT = "Q>Q>"

      def find(intervals, parser)
        parser.fetch_blocks(fetch_list(intervals))
      end

      ## keys:
      ##  0xFF<chrom><bin><start><end>
      ## values:
      ##  <offset>:<length>
      ##
      ## bin: 5 digits (could be 4, or 4 hex, or 16 bits)
      ## start, end: 10 digits (could do 8 hex, or 32 bits)
      ## offset: hex
      ## length: hex
      ##
      ## ex: 01195:80082334:80082367              (23 bytes)
      ## =    04AB:4C5F59E:4C5F5BF                (20 bytes)
      ## = [1195, 80082334, 80082367].pack("SLL") (10 bytes)
      ##
      ## Intervals: Ruby 1...9 intervals (three dots) are half-open
      ##
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
      def fetch_list(intervals)
        to_fetch = []
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
          i.bin_all.each { |bin| bin_intervals[bin] << i }
        end
        db.cursor_process do |cur|
          bin_intervals.each do |bin, bin_intervals_raw|
            bin_intervals = bin_intervals_raw.sort_by { |i| i.zero_start }
            # compute the start and end of all intervals of interest
            spanning_start = bin_intervals.first.zero_start
            spanning_end = bin_intervals.collect {|i| i.zero_end}.sort.last
            # scan from the start of the bin
            cur.jump(bin_start_prefix(chrom_id, bin))
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
                c_int = GenomicInterval.zero_based(chrom, c_start, c_end)
                if bin_intervals.find { |i| i.overlapped?(c_int) }
                  to_fetch << pair[1].unpack(VAL_FMT)
                end
              end
            end
          end # bin_intervals.each
        end # #cursor_process
        return to_fetch
      end # #fetch_list

      def unpack_key(ks)
        ks.unpack(KEY_FMT)
      end

      def bin_start_prefix(chrom_id, bin)
        [0xFF, chrom_id, bin].pack(CHROM_BIN_PREFIX_FMT)
      end

      def self.build(parser, path)
        idx = self.new(path)
        idx.build_default(parser)
        return idx
      end

      def initialize(path)
        if (path.size > 1) and File.exist?(path)
          mode = KyotoCabinet::DB::OREADER
        else
          mode = KyotoCabinet::DB::OWRITER | KyotoCabinet::DB::OCREATE
        end
        @db = KyotoCabinet::DB.new
        @path = path
        unless db.open(path, mode)
          raise "Could not open DB file!"
        end
      end

      def build_default(parser)
        first_block = parser.parse_block
        ref_seq = first_block.sequences.first.source
        @index_sequences = { ref_seq => 0 }
        index_block(first_block)
        parser.each_block { |b| index_block(b) }
      end

      def index_block(block)
        entries_for(block).each do |k, v|
          db.set(k, v)
        end
      end

      def entries_for(block)
        e = []
        block.sequences.each do |seq|
          seq_id = index_sequences[seq.source]
          next unless seq_id
          seq_end = seq.start + seq.size
          bin = Bio::Ucsc::UcscBin.bin_from_range(seq.start, seq_end)
          key = [255, seq_id, bin, seq.start, seq_end].pack(KEY_FMT)
          val = [block.offset, block.size].pack(VAL_FMT)
          e << [key, val]
        end
        return e
      end
    end

    class SQLiteIndex

      attr_accessor :sequence, :db, :table_name

      def self.build(parser, idx_path)
        if File.exist? idx_path
          raise "Cannot build index: #{idx_path} already exists!"
        end
        idx = self.new(idx_path)
        idx.create_schema
        idx.build_default(parser)
        return idx
      end

      def self.open(path)
        unless File.exist? path.to_s
          raise "Cannot open index: #{path} does not exist!"
        end
        idx = self.new(path)
        idx.load
        return idx
      end

      def initialize(path)
        # TODO: for JRuby we need DBI:jdbc:sqlite:<path>
        # and 'driver' => 'org.sqlite.JDBC'
        @path = path
        @db = DBI.connect("DBI:SQLite3:#{path.to_s}", "", "")
      end

      def find(ranges)
        fetch = fetch_list(ranges)
        # TODO: fetch blocks from MAF file
      end

      def fetch_list(ranges)
        # returns an array of
        # [pos, [start, end]]
        # where pos is the block's start offset
        # and start and end are the expected reference sequence start
        # and end positions
        fetch_h = {}
        ranges.each do |range|
          # XXX: check whether these are properly half-open
          bins = Bio::Ucsc::UcscBin.bin_all(range.begin, range.end)
          bin_list = bins.join(',')
          query = <<-EOF
          SELECT pos, start, end
          FROM #{table_name}
          WHERE bin IN (#{bin_list})
          AND (end BETWEEN #{range.begin} AND #{range.end}
               OR #{range.end} BETWEEN start AND end)
          EOF
          db.execute(query) do |stmt|
            stmt.fetch do |row|
              fetch_h[row[0]] = [row[1], row[2]]
            end
          end
        end
        return fetch_h.to_a.sort_by! { |e| e[0] }
      end

      def load
        unless count_tables("metadata") == 1
          raise "#{@path} is not a usable index database!"
        end
        row = db.select_one("select value from metadata where key = 'ref_seq'")
        unless row
          raise "#{@path} is not a usable index database: missing reference sequence name!"
        end
        @sequence = row[0]
        select_table_name!
      end

      def build_default(parser)
        first_block = parser.parse_block
        @sequence = first_block.sequences.first.source
        db.do("insert into metadata (key, value) values (?, ?)",
              "ref_seq", @sequence)
        select_table_name!
        create_index_table

        insert_index_row(index_tuple(first_block))
        parser.each_block do |b|
          insert_index_row(index_tuple(b))
        end
        build_covering_index!
      end

      def build_covering_index!
        db.do(<<-EOF)
CREATE INDEX idx_#{table_name}
ON #{table_name} (bin, start, end, pos)
        EOF
      end

      def insert_index_row(tuple)
        @db.do(<<-EOF, *tuple)
INSERT INTO #{table_name} (bin, start, end, pos)
VALUES (?, ?, ?, ?)
EOF
      end

      def select_table_name!
        tname_base = sequence.gsub(/[^a-zA-Z0-9]/, '_')
        @table_name = "seq_#{tname_base}"
      end

      def create_index_table
        db.do(<<-EOF)
CREATE TABLE #{table_name} (
    bin integer not null,
    start integer not null,
    end integer not null,
    pos integer not null
)
        EOF
      end

      def create_schema
        db.do(<<-EOF)
CREATE TABLE metadata (
    key varchar(32) primary key not null,
    value varchar(32) not null
)
        EOF
      end

      def count_tables(name_pat)
        if name_pat =~ /%/
          op = 'LIKE'
        else
          op = '='
        end
        query = <<-"EOF"
           SELECT COUNT(*)
           FROM sqlite_master
           WHERE name #{op} '#{name_pat}'
           AND type = 'table'
        EOF
        res = db.select_one(query)
        ## XXX: what? DBI is returning a string here?
        return res[0].to_i
      end

      def index_tuple(block)
        seq = block.sequences.find { |s| s.source == sequence }
        seq_end = seq.start + seq.size - 1
        bin = Bio::Ucsc::UcscBin.bin_from_range(seq.start, seq_end)
        return [bin, seq.start, seq_end, block.offset]
      end
    end
    
  end
  
end
