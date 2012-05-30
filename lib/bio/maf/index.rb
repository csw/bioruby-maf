require 'dbi'

require 'bio-ucsc-api'

module Bio

  module MAF

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
      end

      def initialize(path)
        # TODO: for JRuby we need DBI:jdbc:sqlite:<path>
        # and 'driver' => 'org.sqlite.JDBC'
        @path = path
        @db = DBI.connect("DBI:SQLite3:#{path.to_s}", "", "")
      end

      def load
        unless count_tables("metadata") == 1
          raise "#{@path} is not a usable index database!"
        end
        #row = db.select_one("select value from metadata where key = 'ref_seq'")
        #unless row
        #  raise "Could not fetch reference sequence name!"
        #end
        #@sequence = row[0]
      end

      def build_default(parser)
        first_block = parser.parse_block
        @sequence = first_block.sequences.first.source
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
        seq_end = seq.start + seq.size
        bin = Bio::Ucsc::UcscBin.bin_from_range(seq.start, seq_end)
        return [bin, seq.start, seq_end, block.pos]
      end
    end
    
  end
  
end
