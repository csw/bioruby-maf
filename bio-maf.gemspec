# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = "bio-maf"
  s.version = "0.1.0"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["Clayton Wheeler"]
  s.date = "2012-06-29"
  s.description = "Multiple Alignment Format parser for BioRuby."
  s.email = "cswh@umich.edu"
  s.executables = ["maf_count", "maf_dump_blocks", "maf_extract_ranges_count", "maf_index", "maf_parse_bench", "maf_to_fasta", "maf_write", "random_ranges"]
  s.extra_rdoc_files = [
    "LICENSE.txt",
    "README.md"
  ]
  s.files = [
    ".document",
    ".simplecov",
    ".travis.yml",
    ".yardopts",
    "DEVELOPMENT.md",
    "Gemfile",
    "LICENSE.txt",
    "README.md",
    "Rakefile",
    "VERSION",
    "benchmarks/dispatch_bench",
    "benchmarks/iter_bench",
    "benchmarks/read_bench",
    "benchmarks/sort_bench",
    "benchmarks/split_bench",
    "bin/maf_count",
    "bin/maf_dump_blocks",
    "bin/maf_extract_ranges_count",
    "bin/maf_index",
    "bin/maf_parse_bench",
    "bin/maf_to_fasta",
    "bin/maf_write",
    "bin/random_ranges",
    "features/maf-indexing.feature",
    "features/maf-output.feature",
    "features/maf-parsing.feature",
    "features/maf-querying.feature",
    "features/maf-to-fasta.feature",
    "features/step_definitions/convert_steps.rb",
    "features/step_definitions/index_steps.rb",
    "features/step_definitions/output_steps.rb",
    "features/step_definitions/parse_steps.rb",
    "features/step_definitions/query_steps.rb",
    "features/step_definitions/ucsc_bin_steps.rb",
    "features/support/env.rb",
    "features/ucsc-bins.feature",
    "lib/bio-maf.rb",
    "lib/bio-maf/maf.rb",
    "lib/bio/maf.rb",
    "lib/bio/maf/index.rb",
    "lib/bio/maf/parser.rb",
    "lib/bio/maf/struct.rb",
    "lib/bio/maf/writer.rb",
    "lib/bio/ucsc.rb",
    "lib/bio/ucsc/genomic-interval-bin.rb",
    "lib/bio/ucsc/ucsc_bin.rb",
    "man/.gitignore",
    "man/maf_index.1",
    "man/maf_index.1.markdown",
    "man/maf_index.1.ronn",
    "man/maf_to_fasta.1",
    "man/maf_to_fasta.1.ronn",
    "spec/bio/maf/index_spec.rb",
    "spec/bio/maf/parser_spec.rb",
    "spec/bio/maf/struct_spec.rb",
    "spec/spec_helper.rb",
    "test/data/big-block.maf",
    "test/data/chr22_ieq.maf",
    "test/data/chrY-1block.maf",
    "test/data/empty",
    "test/data/empty.db",
    "test/data/mm8_chr7_tiny.kct",
    "test/data/mm8_chr7_tiny.maf",
    "test/data/mm8_mod_a.maf",
    "test/data/mm8_single.maf",
    "test/data/mm8_subset_a.maf",
    "test/data/t1-bad1.maf",
    "test/data/t1.fasta",
    "test/data/t1.maf",
    "test/data/t1a.maf",
    "test/helper.rb",
    "test/test_bio-maf.rb",
    "travis-ci/install_kc",
    "travis-ci/install_kc_java",
    "travis-ci/report_errors"
  ]
  s.homepage = "http://github.com/csw/bioruby-maf"
  s.licenses = ["MIT"]
  s.require_paths = ["lib"]
  s.rubygems_version = "1.8.24"
  s.summary = "MAF parser for BioRuby"

  s.specification_version = 3

  if RUBY_PLATFORM == 'java'
    s.platform = 'java'
  end

  s.add_runtime_dependency('bio-bigbio', [">= 0"])
  s.add_runtime_dependency('bio-genomic-interval', ["~> 0.1.2"])
  if RUBY_PLATFORM == 'java'
    s.add_runtime_dependency('kyotocabinet-java', ["~> 0.2.0"])
  else    
    s.add_runtime_dependency('kyotocabinet-ruby', ["~> 1.27.1"])
  end

end
