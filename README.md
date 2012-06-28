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

```sh
    gem install bio-maf
```

## Usage

```ruby
    require 'bio-maf'

    p = Bio::MAF::Parser.new(path)
    header = p.header
    p.each_block do |block|
      block.sequences.each do |seq|
        do_something(seq)
      end
    end
```

Man pages for command line tools:

* [maf_index](http://csw.github.com/bioruby-maf/man/maf_index.1.html)
* [maf_to_fasta](http://csw.github.com/bioruby-maf/man/maf_to_fasta.1.html)

The API doc is online. For more code examples see the
[RSpec](https://github.com/csw/bioruby-maf/tree/master/spec/bio/maf)
and
[Cucumber](https://github.com/csw/bioruby-maf/tree/master/features)
test files in the source tree.

Also, the scripts in the
[bin](https://github.com/csw/bioruby-maf/tree/master/bin) directory
provide good worked examples of how to use the existing parsing API.
        
## Project home page

Information on the source tree, documentation, examples, issues and
how to contribute, see

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

