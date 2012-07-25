# Please require your code below, respecting the naming conventions in the
# bioruby directory tree.
#
# For example, say you have a plugin named bio-plugin, the only uncommented
# line in this file would be 
#
#   require 'bio/bio-plugin/plugin'
#
# In this file only require other files. Avoid other source code.

require 'bio-logger'
log = Bio::Log::LoggerPlus.new('bio-maf')
log.outputters = Bio::Log::Outputter.stderr
log.level = Bio::Log::WARN

require 'bio/ucsc'
require 'bio/maf'
