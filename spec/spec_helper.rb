unless ENV.has_key?('TRAVIS') || RUBY_PLATFORM == 'java'
  begin
    require 'simplecov'
  rescue LoadError
    $stderr.puts "WARNING: could not require 'simplecov': #{$!}"
  end
end

require 'rspec'
require 'pathname'

require 'bio-maf'

TestData = Pathname.new(__FILE__) + '../../test/data'
