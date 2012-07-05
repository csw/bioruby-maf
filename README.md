# bio-maf

[![Build Status](https://secure.travis-ci.org/csw/bioruby-maf.png)](http://travis-ci.org/csw/bioruby-maf)

This is a plugin for [BioRuby](http://bioruby.open-bio.org/) adding
support for the
[Multiple Alignment Format](http://genome.ucsc.edu/FAQ/FAQformat#format5)
(MAF), used in bioinformatics to store whole-genome sets of multiple
sequence alignments.

Ultimately it will provide indexed and sequential access to MAF data,
as well as performing various manipulations on it and writing modified
MAF files. So far, it only supports simple sequential parsing.

For more information, see the
[project wiki](https://github.com/csw/bioruby-maf/wiki).

Developer documentation generated with YARD is available at
[rubydoc.info](http://rubydoc.info/github/csw/bioruby-maf/).

This is being developed by Clayton Wheeler as
[part of](http://www.bioruby.org/wiki/Google_Summer_of_Code) the
Google Summer of Code 2012, under the auspices of the Open
Bioinformatics Foundation. The development
[blog](http://csw.github.com/bioruby-maf/) may be of interest.

## Dependencies

[Kyoto Cabinet][] is a database library, required for building MAF
indexes. Install the core library in the appropriate way for your
platform, as documented [here][].

[Kyoto Cabinet]: http://fallabs.com/kyotocabinet/
[here]: https://github.com/csw/bioruby-maf/wiki/Kyoto-Cabinet

If you're using MRI, the [kyotocabinet-ruby][] gem will be used to
interact with Kyoto Cabinet. For best performance, however, you should
really consider using JRuby. On JRuby, the [kyotocabinet-java][] gem
will be used instead; this builds a Java library using JNI to call
into Kyoto Cabinet. Please file a [bug report][] if you encounter
problems building or using this gem, which is still fairly new.

[kyotocabinet-ruby]: https://rubygems.org/gems/kyotocabinet-ruby
[kyotocabinet-java]: https://github.com/csw/kyotocabinet-java
[bug report]: https://github.com/csw/kyotocabinet-java/issues


## Installation

`bio-maf` is now published as a Ruby [gem](https://rubygems.org/gems/bio-maf).

    $ gem install bio-maf

## Usage

### Create an index on a MAF file

Much of the functionality of this library relies on an index. You can
create one with [maf_index(1)][], like so:

[maf_index(1)]: http://csw.github.com/bioruby-maf/man/maf_index.1.html


    $ maf_index test/data/mm8_chr7_tiny.maf /tmp/mm8_chr7_tiny.kct
    
Or programmatically:

    require 'bio-maf'
    parser = Bio::MAF::Parser.new("test/data/mm8_chr7_tiny.maf")
    idx = Bio::MAF::KyotoIndex.build(parser, "/tmp/mm8_chr7_tiny.kct")

### Extract blocks from an indexed MAF file, by genomic interval

Refer to [`mm8_chr7_tiny.maf`](https://github.com/csw/bioruby-maf/blob/master/test/data/mm8_chr7_tiny.maf).


    require 'bio-maf'
    parser = Bio::MAF::Parser.new('test/data/mm8_chr7_tiny.maf')
    idx = Bio::MAF::KyotoIndex.open('test/data/mm8_chr7_tiny.kct')

    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 80082592, 80082766)]
    idx.find(q, parser).each do |block|
      ref_seq = block.sequences[0]
      puts "Matched block at #{ref_seq.start}, #{ref_seq.size} bases"
    end

    # => Matched block at 80082592, 121 bases
    # => Matched block at 80082713, 54 bases

### Filter species returned in alignment blocks

    require 'bio-maf'
    parser = Bio::MAF::Parser.new('test/data/mm8_chr7_tiny.maf')
    idx = Bio::MAF::KyotoIndex.open('test/data/mm8_chr7_tiny.kct')

    parser.sequence_filter = { :only_species => %w(hg18 mm8 rheMac2) }
    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 80082592, 80082766)]
    blocks = idx.find(q, parser)
    block = blocks.first
    puts "Block has #{block.sequences.size} sequences."

    # => Block has 3 sequences.

### Extract blocks matching certain conditions

See also the [Cucumber feature][] and [step definitions][] for this.

[Cucumber feature]: https://github.com/csw/bioruby-maf/blob/master/features/maf-querying.feature
[step definitions]: https://github.com/csw/bioruby-maf/blob/master/features/step_definitions/query_steps.rb

#### Match only blocks with all specified species

    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 80082471, 80082730)]
    filter = { :with_all_species => %w(panTro2 loxAfr1) }
    n_blocks = idx.find(q, parser, filter).count
    # => 1

#### Match only blocks with a certain number of sequences

    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 80082767, 80083008)]
    filter = { :at_least_n_sequences => 6 }
    n_blocks = idx.find(q, parser, filter).count
    # => 1

#### Match only blocks within a text size range

    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 0, 80100000)]
    filter = { :min_size => 72, :max_size => 160 }
    n_blocks = idx.find(q, parser, filter).count
    # => 3

### Process each block in a MAF file

    require 'bio-maf'
    p = Bio::MAF::Parser.new('test/data/mm8_chr7_tiny.maf')
    puts "MAF version: #{p.header.version}"
    # => MAF version: 1

    p.parse_blocks.each do |block|
      block.sequences.each do |seq|
        do_something(seq)
      end
    end

### Parse empty ('e') lines

Refer to [`chr22_ieq.maf`](https://github.com/csw/bioruby-maf/blob/master/test/data/chr22_ieq.maf).

    require 'bio-maf'
    p = Bio::MAF::Parser.new('test/data/chr22_ieq.maf',
                             :parse_empty => false)
    block = p.parse_block
    block.sequences.size
    # => 3

    p = Bio::MAF::Parser.new('test/data/chr22_ieq.maf',
                             :parse_empty => true)
    block = p.parse_block
    block.sequences.size
    # => 4
    block.sequences.find { |s| s.empty? }
    # => #<Bio::MAF::EmptySequence:0x007fe1f39882d0 
    #      @source="turTru1.scaffold_109008", @start=25049,
    #      @size=1601, @strand=:+, @src_size=50103, @text=nil,
    #      @status="I"> 


### Command line tools

Man pages for command line tools:

* [`maf_index(1)`](http://csw.github.com/bioruby-maf/man/maf_index.1.html)
* [`maf_to_fasta(1)`](http://csw.github.com/bioruby-maf/man/maf_to_fasta.1.html)
* [`maf_tile(1)`](http://csw.github.com/bioruby-maf/man/maf_tile.1.html)

### Other documentation

Also see the [API documentation][]. For more code examples see the
[RSpec][] and [Cucumber][] test files in the source tree.

[API documentation]: http://rubydoc.info/github/csw/bioruby-maf/
[RSpec]: https://github.com/csw/bioruby-maf/tree/master/spec/bio/maf
[Cucumber]: https://github.com/csw/bioruby-maf/tree/master/features 

Also, the scripts in the
[bin](https://github.com/csw/bioruby-maf/tree/master/bin) directory
provide good worked examples of how to use the existing parsing API.
        
## Project home page

For information on the source tree, documentation, examples, issues
and how to contribute, see

  <http://github.com/csw/bioruby-maf>

The BioRuby community is on IRC server: irc.freenode.org, channel: #bioruby.

## Cite

If you use this software, please cite one of
  
* [BioRuby: bioinformatics software for the Ruby programming language](http://dx.doi.org/10.1093/bioinformatics/btq475)
* [Biogem: an effective tool-based approach for scaling up open source software development in bioinformatics](http://dx.doi.org/10.1093/bioinformatics/bts080)

## Biogems.info

This Biogem will be published at [#bio-maf](http://biogems.info/index.html)

## Copyright

Copyright (c) 2012 Clayton Wheeler. See LICENSE.txt for further details.

