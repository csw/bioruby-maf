require 'dbi'

module Bio

  module MAF

    class SQLiteIndex

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

      
    end
    
  end
  
end
