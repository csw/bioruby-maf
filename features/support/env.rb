require 'bundler/setup'

unless ENV.has_key?('TRAVIS') || RUBY_PLATFORM == 'java'
  begin
    require 'simplecov'
  rescue LoadError
    $stderr.puts "WARNING: could not require 'simplecov': #{$!}"
  end
end

require 'pathname'
require 'tempfile'

lib_dir = File.expand_path('../../../lib', __FILE__)
$LOAD_PATH << lib_dir
ENV['RUBYLIB'] = lib_dir

require 'bio-maf'

$test_data = Pathname.new 'test/data'
