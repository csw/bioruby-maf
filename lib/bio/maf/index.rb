require 'sqlite3'

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
        @db = SQLite3::Database.new(path.to_s)
      end
      
    end
    
  end
  
end
