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

## Performance

This parser performs best under [JRuby][], particularly with Java
7. See the [Performance][] wiki page for more information. For best
results, use JRuby in 1.9 mode with the ObjectProxyCache disabled:

[JRuby]: http://jruby.org/
[Performance]: https://github.com/csw/bioruby-maf/wiki/Performance

    $ export JRUBY_OPTS='--1.9 -Xji.objectProxyCache=false'

Many parsing modes are multithreaded. Under JRuby, it will default to
using one parser thread per available core, but if desired this can be
configured with the `:threads` parser option.

Ruby 1.9.3 is fully supported, but does not perform as well,
especially since its concurrency features are not useful for this
workload.

## Usage

### Create an index on a MAF file

Much of the functionality of this library relies on an index. You can
create one with [maf_index(1)][], like so:

[maf_index(1)]: http://csw.github.com/bioruby-maf/man/maf_index.1.html


    $ maf_index test/data/mm8_chr7_tiny.maf /tmp/mm8_chr7_tiny.kct
    
Or programmatically:

```ruby
    require 'bio-maf'
    parser = Bio::MAF::Parser.new("test/data/mm8_chr7_tiny.maf")
    idx = Bio::MAF::KyotoIndex.build(parser, "/tmp/mm8_chr7_tiny.kct")
```

### Extract blocks from an indexed MAF file, by genomic interval

Refer to [`mm8_chr7_tiny.maf`](https://github.com/csw/bioruby-maf/blob/master/test/data/mm8_chr7_tiny.maf).

    require 'bio-maf'
    access = Bio::MAF::Access.maf_dir('test/data')

    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 80082592, 80082766)]
    access.find(q) do |block|
      ref_seq = block.sequences[0]
      puts "Matched block at #{ref_seq.start}, #{ref_seq.size} bases"
    end

    # => Matched block at 80082592, 121 bases
    # => Matched block at 80082713, 54 bases

Or, equivalently, one can work with a specific MAF file and index directly:

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

### Extract alignment blocks truncated to a given interval

Given a genomic interval of interest, one can also extract only the
subsets of blocks that intersect with that interval, using the
`#slice` method like so:

    require 'bio-maf'
    access = Bio::MAF::Access.maf_dir('test/data')
    int = Bio::GenomicInterval.zero_based('mm8.chr7', 80082350, 80082380)
    blocks = access.slice(int).to_a
    puts "Got #{blocks.size} blocks, first #{blocks.first.ref_seq.size} base pairs."
    # => Got 2 blocks, first 18 base pairs.

### Filter species returned in alignment blocks

    require 'bio-maf'
    access = Bio::MAF::Access.maf_dir('test/data')

    access.sequence_filter = { :only_species => %w(hg18 mm8 rheMac2) }
    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 80082592, 80082766)]
    blocks = access.find(q)
    block = blocks.first
    puts "Block has #{block.sequences.size} sequences."

    # => Block has 3 sequences.

### Extract blocks matching certain conditions

See also the [Cucumber feature][] and [step definitions][] for this.

[Cucumber feature]: https://github.com/csw/bioruby-maf/blob/master/features/maf-querying.feature
[step definitions]: https://github.com/csw/bioruby-maf/blob/master/features/step_definitions/query_steps.rb

#### Match only blocks with all specified species

    access = Bio::MAF::Access.maf_dir('test/data')
    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 80082471, 80082730)]
    access.block_filter = { :with_all_species => %w(panTro2 loxAfr1) }
    n_blocks = access.find(q).count
    # => 1

#### Match only blocks with a certain number of sequences

    access = Bio::MAF::Access.maf_dir('test/data')
    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 80082767, 80083008)]
    access.block_filter = { :at_least_n_sequences => 6 }
    n_blocks = access.find(q).count
    # => 1

#### Match only blocks within a text size range

    access = Bio::MAF::Access.maf_dir('test/data')
    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 0, 80100000)]
    access.block_filter = { :min_size => 72, :max_size => 160 }
    n_blocks = access.find(q).count
    # => 3

### Process each block in a MAF file

    require 'bio-maf'
    p = Bio::MAF::Parser.new('test/data/mm8_chr7_tiny.maf')
    puts "MAF version: #{p.header.version}"
    # => MAF version: 1

    p.each_block do |block|
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

Such options can also be set on a Bio::MAF::Access object:

    require 'bio-maf'
    access = Bio::MAF::Access.maf_dir('test/data')
    access.parse_options[:parse_empty] = true

### Remove gaps from parsed blocks

After filtering out species with
[`Parser#sequence_filter`](#filter-species-returned-in-alignment-blocks),
gaps may be left where there was an insertion present only in
sequences that were filtered out. Such gaps can be removed by setting
the `:remove_gaps` parser option:

    require 'bio-maf'
    access = Bio::MAF::Access.maf_dir('test/data')
    access.parse_options[:remove_gaps] = true

### Join blocks after filtering together

Similarly, filtering out species may remove a species which had caused
two adjacent alignment blocks to be split. By enabling the
`:join_blocks` parser option, such blocks can be joined together:

    require 'bio-maf'
    access = Bio::MAF::Access.maf_dir('test/data')
    access.parse_options[:join_blocks] = true

See the [Cucumber feature][] for more details.

[Cucumber feature]: https://github.com/csw/bioruby-maf/blob/master/features/block-joining.feature

### Extract bio-alignment representations of blocks

When the `:as_bio_alignment` parser option is given, blocks will be
returned as [Bio::BioAlignment::Alignment][] objects as used in the
[bio-alignment] Biogem. This offers a great deal of built-in
functionality for column-wise operations, alignment manipulation, and
more.

[Bio::BioAlignment::Alignment]: http://rdoc.info/gems/bio-alignment/Bio/BioAlignment/Alignment
[bio-alignment]: https://github.com/pjotrp/bioruby-alignment

    require 'bio-maf'
    access = Bio::MAF::Access.maf_dir('test/data')
    access.parse_options[:as_bio_alignment] = true
    q = [Bio::GenomicInterval.zero_based('mm8.chr7', 80082592, 80082766)]
    access.find(q) do |aln|
      col = aln.columns[3]
      puts "bases in column 3: #{col}"
    end

### Tile blocks together over an interval

Extracts alignment blocks overlapping the given genomic interval and
constructs a single alignment block covering the entire interval for
the specified species. Optionally, any gaps in coverage of the MAF
file's reference sequence can be filled in from a FASTA sequence
file. See the Cucumber [feature][] for examples of output, and also
the
[`maf_tile(1)`](http://csw.github.com/bioruby-maf/man/maf_tile.1.html)
man page.

[feature]: https://github.com/csw/bioruby-maf/blob/master/features/tiling.feature

    require 'bio-maf'
    access = Bio::MAF::Access.maf_dir('test/data')
    interval = Bio::GenomicInterval.zero_based('mm8.chr7',
                                               80082334,
                                               80082468)
    access.tile(interval) do |tiler|
      # reference is optional
      tiler.reference = 'reference.fa.gz'
      tiler.species = %w(mm8 rn4 hg18)
      # species_map is optional
      tiler.species_map = {
        'mm8' => 'mouse',
        'rn4' => 'rat',
        'hg18' => 'human'
      }
      tiler.write_fasta($stdout)
    end

### Command line tools

Man pages for command line tools:

* [`maf_index(1)`](http://csw.github.com/bioruby-maf/man/maf_index.1.html)
* [`maf_to_fasta(1)`](http://csw.github.com/bioruby-maf/man/maf_to_fasta.1.html)
* [`maf_tile(1)`](http://csw.github.com/bioruby-maf/man/maf_tile.1.html)

With [gem-man](https://github.com/defunkt/gem-man) installed, these
can be read with:

    $ gem man bio-maf

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

This Biogem is published at [biogems.info](http://biogems.info/index.html#bio-maf).

## Copyright

Copyright (c) 2012 Clayton Wheeler. See LICENSE.txt for further details.

