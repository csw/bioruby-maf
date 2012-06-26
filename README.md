# bio-maf

<!--
[![Build Status](https://secure.travis-ci.org/csw/bioruby-maf.png)](http://travis-ci.org/csw/bioruby-maf)
-->

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
platform (with Homebrew on Mac OS X: `brew install kyoto-cabinet`).

[Kyoto Cabinet]: http://fallabs.com/kyotocabinet/

If you're using MRI, the `kyotocabinet-ruby` gem is included as a
dependency. A quick install for MRI of the current version of Kyoto
Cabinet (1.2.76) on Debian-like systems is

```sh
  sudo apt-get install liblzo2-dev lzma-dev
  ./configure --enable-zlib --enable-lzo --enable-lzma --prefix=/usr --disable-atomic
  make
  sudo make install
```

(you may need to comment out the lines in configure as suggested by
this
[patch](https://gitorious.org/kyotocabinet-debian/kyotocabinet/blobs/master/debian/patches/0001-disable-march-native.patch).

For best performance, however, you should really consider using
JRuby. Kyoto Cabinet support under JRuby requires that you build and
install the [kyotocabinet-java][] JNI library on your system. By
default, its `make install` target will install `libjkyotocabinet.so`
in `/usr/local/lib`, which is on the JNI load path
(`java.library.path`) on Mac OS X. Ubuntu users will need to copy this
library elsewhere: `/usr/lib/jni` is one location that should work.

[kyotocabinet-java]: http://fallabs.com/kyotocabinet/javapkg/

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

