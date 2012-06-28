# Development guide

Here are notes on less obvious aspects of the development process for
this library.

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
