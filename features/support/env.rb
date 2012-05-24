unless ENV.has_key? 'TRAVIS'
  begin
    require 'simplecov'
  rescue
    $stderr.puts "WARNING: could not require 'simplecov': $!"
  end
end

require 'pathname'
require 'tempfile'

$LOAD_PATH << File.expand_path('../../../lib', __FILE__)

require 'bio-maf'

$test_data = Pathname.new 'test/data'
