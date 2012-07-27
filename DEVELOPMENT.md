# Development guide

Here are notes on less obvious aspects of the development process for
this library.

## Gem build / tagging / release

This now uses [rubygems-tasks][] for building and releasing gems.

[rubygems-tasks]: https://github.com/postmodern/rubygems-tasks

We build two gem platform variants: a 'default' one for MRI with no
platform set, and a JRuby one with `platform = 'java'`. These get
built as `bio-maf-X.Y.Z.gem` and `bio-maf-X.Y.Z-java.gem`. At least
for now, this is done by running `gem release` separately under JRuby
and MRI. SCM tagging and pushing is done under MRI only, but the gems
will be built and pushed to rubygems.org separately under each
platform.

The version is simply set by hand in `bio-maf.gemspec`. Don't forget
to increment it!

First, verify that you are on the `master` branch:

    $ git branch

Testing the build:

    $ rake build
    $ rake install

Release:

    $ rvm use 1.9.3@bioruby-maf
    $ rake release
    $ rvm use jruby-1.6.7.2@bioruby-maf
    $ rake release

## kyotocabinet-java

Running `bio-maf` on JRuby requires the [kyotocabinet-java][] gem, a
wrapper around the Kyoto Cabinet Java interface providing a Ruby API
compatible with the standard Kyoto Cabinet Ruby API.

[kyotocabinet-java]: https://github.com/csw/kyotocabinet-java

## Man pages

Man pages are developed with [ronn][] and live in `man/`; see
[maf_index.1.ronn][] for an example. The generated man pages,
e.g. `maf_index.1`, are added to Git for [gem-man][] support.

[ronn]: https://github.com/rtomayko/ronn
[gem-man]: https://github.com/defunkt/gem-man
[maf_index.1.ronn]: https://github.com/csw/bioruby-maf/blob/master/man/maf_index.1.ronn

HTML and roff versions are built with:

    $ rake man

The HTML versions are published through Octopress to Github Pages,
e.g. <http://csw.github.com/bioruby-maf/man/maf_index.1.html>. This is
a separate step, and necessarily dependent on the local filesystem
layout. Specifically, there must be an `octopress` directory at the
same level as `bioruby-maf`, containing a checked-out copy of
<https://github.com/csw/bioruby-maf-blog>. Then, to publish the man
pages, run:

    $ rake man:publish
    
After this, in that Octopress instance, run:

    $ rake deploy
