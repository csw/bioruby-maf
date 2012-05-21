require 'pathname'
require 'tempfile'

$LOAD_PATH << File.expand_path('../../../lib', __FILE__)

require 'bio-maf'

$test_data = Pathname.new 'test/data'
