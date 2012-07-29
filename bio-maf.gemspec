# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = "bio-maf"
  s.version = "0.3.2"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Clayton Wheeler"]
  s.date = "2012-07-26"
  s.description = "Multiple Alignment Format parser for BioRuby."
  s.email = "cswh@umich.edu"
  s.extra_rdoc_files = [
    "LICENSE.txt",
    "README.md"
                       ]
  s.files         = `git ls-files`.split("\n")
  s.test_files    = `git ls-files -- {test,spec,features}/*`.split("\n")
  s.executables   = `git ls-files -- bin/*`.split("\n").map {
    |f| File.basename(f)
  }

  s.homepage = "http://github.com/csw/bioruby-maf"
  s.licenses = ["MIT"]
  s.require_paths = ["lib"]
  s.rubygems_version = "1.8.24"
  s.summary = "MAF parser for BioRuby"

  s.specification_version = 3

  if RUBY_PLATFORM == 'java'
    s.platform = 'java'
  end

  s.add_runtime_dependency('bio-alignment', ["~> 0.0.7"])
  s.add_runtime_dependency('bio-bgzf', ["~> 0.1.1"])
  s.add_runtime_dependency('bio-genomic-interval', ["~> 0.1.2"])
  s.add_runtime_dependency('bio-logger', ["~> 1.0.1"])
  if RUBY_PLATFORM == 'java'
    s.add_runtime_dependency('kyotocabinet-java', ["~> 0.3.0"])
  else    
    s.add_runtime_dependency('kyotocabinet-ruby', ["~> 1.27.1"])
  end

end
