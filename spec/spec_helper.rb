unless ENV.has_key? 'TRAVIS'
  begin
    require 'simplecov'
  rescue
    $stderr.puts "WARNING: could not require 'simplecov': $!"
  end
end

require 'rspec'
require 'pathname'

require 'bio-maf'

TestData = Pathname.new(__FILE__) + '../../test/data'
