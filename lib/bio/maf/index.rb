require 'dbi'

require 'bio-ucsc-api'

module Bio

  module MAF

    class SQLiteIndex

      attr_accessor :sequence, :db

      def self.build(parser, idx_path)
        if File.exist? idx_path
          raise "Cannot build index: #{idx_path} already exists!"
        end
        idx = self.new(idx_path)

      end

      def initialize(path)
        # TODO: for JRuby we need DBI:jdbc:sqlite:<path>
        # and 'driver' => 'org.sqlite.JDBC'
        @db = DBI.connect("DBI:SQLite3:#{path.to_s}", "", "")
      end

      def index_tuple(block)
        seq = block.sequences.find { |s| s.source == sequence }
        seq_end = seq.start + seq.size
        bin = Bio::Ucsc::UcscBin.bin_from_range(seq.start, seq_end)
        return [bin, seq.start, seq_end]
      end
    end
    
  end
  
end
